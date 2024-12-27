fit.iv.owl = function(itr_formula, iv_formula, trt_formula, mean_formula = NULL, trt, iv, outcome, dat, 
                      na.rm = T, minimize = F, center.outcome = F, 
                      kernel = "linear", cross = 10, scale = T, seed, ...){ 
  
  
  set.seed(seed)
  ## Error handling
  if (!inherits(itr_formula, "formula")) stop("itr_formula must be a formula object.")
  if (!inherits(iv_formula, "formula")) stop("iv_formula must be a formula object.")
  if (!inherits(trt_formula, "formula")) stop("trt_formula must be a formula object.")
  if (!trt %in% names(dat)) stop("The treatment variable is not in the dataset.")
  if (!outcome %in% names(dat)) stop("The outcome variable is not in the dataset.")
  
  
  ## Extract variables from formulas
  itr_vars <- all.vars(itr_formula)
  iv_vars <- all.vars(iv_formula) #variables in the 
  trt_vars <- all.vars(trt_formula)
  
  missing_itr_vars <- setdiff(itr_vars, names(dat))
  missing_iv_vars <- setdiff(iv_vars, names(dat))
  missing_trt_vars <- setdiff(trt_vars, names(dat))
  
  if(length(missing_itr_vars) > 0){
    stop(paste("The following variables in itr_formula are missing from the dataset:", 
               paste(missing_itr_vars, collapse = ", ")))
  }
  
  if(length(missing_iv_vars) > 0){
    stop(paste("The following variables in prop_formula are missing from the dataset:", 
               paste(missing_iv_vars, collapse = ", ")))
  }
  
  if(length(missing_trt_vars) > 0){
    stop(paste("The following variables in trt_formula are missing from the dataset:", 
               paste(missing_trt_vars, collapse = ", ")))
  }
  
  ## Data preparing
  names(dat)[names(dat) ==  iv ] = "Z"
  names(dat)[names(dat) ==  trt] = "A"
  names(dat)[names(dat) == outcome] = "Y" 
  
  
  # spliting data into training and testing set
  dat$train = sample(c(T,F), nrow(dat), prob = c(0.7,0.3), replace = T )
  
  # remove missing data
  if (na.rm) {
    columns_to_check <- unique(c(itr_vars, iv_vars, trt_vars, "A", "Y","Z"))
    dat <- dat[complete.cases(dat[, columns_to_check]), ]
    if (nrow(dat) == 0) {
      stop("All rows were removed due to missing values.")
    }
  }
  
  
  # OWL maximizes Y. If we want to minimize Y, set minimize = T
  if(minimize){
    dat$Y = dat$Y*(-1)
  }
  
  # center the outcome to avoid 0 value impacting the weights
  if(center.outcome){
    dat$Y = dat$Y - mean(dat$Y)
  }
  
  ### Modeling nuisance parameters
  ##  IV model
  
 
  #convert to a formula
  iv.logistic.formula = paste( "as.factor(Z)" , deparse(iv_formula, width.cutoff = 500), collapse = "")  
 
  
  iv.logistic = glm(iv.logistic.formula, data = dat, family = binomial(link = "logit"))
  
  pred.Z = predict(iv.logistic, type = "response", newdata = dat)
  dat$Z.pred = ifelse(dat$Z == 1, pred.Z, 1 - pred.Z)
  
  ## Trt model
  
  #convert to a formula
 
  trt.logistic.formula = paste( "as.factor(A)" , deparse(trt_formula, width.cutoff = 500), collapse = "") 
  
  trt.logistic = glm(trt.logistic.formula, data = dat, family = binomial(link = "logit"))
  
  pred.A = predict(trt.logistic, type = "response", newdata = dat)
  dat$A.pred = ifelse(dat$A == 1, pred.A, 1 - pred.A)
 
  dat_mod = dat 
  #estimating E[A|Z = 1, ...]
  dat_mod$Z = 1
  dat$A.Z1 = predict(trt.logistic, type = "response", newdata = dat_mod)
  #estimating E[A|Z = -1, ...]
  dat_mod$Z = -1
  dat$A.Zn1 = predict(trt.logistic, type = "response", newdata = dat_mod)
  
  #calculating the weight
  dat$delta.pred = dat$A.Z1 - dat$A.Zn1 #calculate the delta term
  dat$weight = with(dat, Y*A*Z/(Z.pred*delta.pred))
  
  dat$fitted.owl = NA
  
  ## calculate the weight
  dat$lab = sign(dat$weight)*dat$A
  
  formula_str <- paste("as.factor(lab)", deparse(itr_formula, width.cutoff = 500), collapse = "")
  formula.owl <- as.formula(formula_str)
  
  ## Model fitting
  mod.owl = WeightSVM::wsvm(formula.owl, data = dat[dat$train,], weight = abs(dat$weight[dat$train]), kernel = "linear", 
                            cross = 10, scale = T) #OWL method
  
  # Predicting 
  fitted.owl = predict(mod.owl, newdata = dat)
  
  
  
  fitted.owl = as.numeric(as.character(fitted.owl)) #opt treatment regimes
  fitted.A = rep(-1, length(fitted.owl)) #always -1 strategy. A is -1 
  fitted.B = rep(1, length(fitted.owl)) #always 1 strategy. B is 1 
  
  
  #return(table(fitted.owl, true.value[!dat$train]))
  
  owl.coef = svm.coef(mod.owl, type = "line")
  names(owl.coef) = c("Intercept", itr_vars)
  
  return(table(fitted.owl))
  
  dat.bt$fitted.maob = rep(-1, nrow(dat.bt))
  dat.bt$fitted.dra = rep(1, nrow(dat.bt))
  # dat$fz.pred = fz.pred 
  # dat$delta.pred = delta.pred
  
  dat.bt$g.p.opt = dat.bt$outcome*dat.bt$A*I(dat.bt$fitted.owl == dat.bt$A)
  dat.bt$g.p.MAOB = dat.bt$outcome*dat.bt$A*I(dat.bt$fitted.maob == dat.bt$A)
  dat.bt$g.p.DRA = dat.bt$outcome*dat.bt$A*I(dat.bt$fitted.dra == dat.bt$A)
  
  ### parametric estimator of gamma
  
  ############ opt ################
  ### estimator of gamma prime
  
  g.p.opt.mod = lm(g.p.opt ~ Z+  updrs_1 +  updrs_2 + updrs_3 + time_from_diagnosis.2nd + levodopa_dose.2nd +
                     gender + age.2nd  + bmi  + Time_since_first_Levo , data = dat.bt)
  temp.dat.opt = dat.bt[,c("Z", "updrs_1", "updrs_2" , "updrs_3" , "time_from_diagnosis.2nd" , "levodopa_dose.2nd" , "gender" , "age.2nd"  ,"bmi", "Time_since_first_Levo" )]
  temp.dat.opt$Z = -1
  dat.bt$g.p.opt.pred = predict(g.p.opt.mod, newdata = temp.dat.opt  )   
  
  #### estimator of gamma
  
  #dat$gamma.pseudo.opt = (dat$g.p.opt - dat$g.p.opt.pred)*dat$Z/dat$fz.pred 
  #dat$gamma.weight = (dat$A - pred.A1.Zn1)*dat$Z/(2*dat$fz.pred)
  dat.bt$gamma.pseudo.opt = (dat.bt$g.p.opt - dat.bt$g.p.opt.pred)*2/(dat.bt$A - dat.bt$pred.A1.Zn1) 
  dat.bt$gamma.weight = ((dat.bt$A - dat.bt$pred.A1.Zn1)*dat.bt$Z/(2*dat.bt$fz.pred))^2
  
  gamma.opt.mod = lm(gamma.pseudo.opt ~  updrs_1 +  updrs_2 + updrs_3 + time_from_diagnosis.2nd + levodopa_dose.2nd +
                       gender + age.2nd  + bmi + Time_since_first_Levo  , data = dat.bt, weights = gamma.weight)
  dat.bt$gamma.opt.pred = predict(gamma.opt.mod, newdata = dat.bt)
  
  #### maob ############################
  ### estimator of gamma prime
  
  g.p.maob.mod = lm(g.p.MAOB ~ Z+  updrs_1 +  updrs_2 + updrs_3 + time_from_diagnosis.2nd + levodopa_dose.2nd +
                      gender + age.2nd  + bmi + Time_since_first_Levo  , data = dat.bt)
  
  dat.bt$g.p.maob.pred = predict(g.p.maob.mod, newdata = temp.dat.opt)   
  
  dat.bt$gamma.pseudo.maob = (dat.bt$g.p.MAOB - dat.bt$g.p.maob.pred)*2/(dat.bt$A - dat.bt$pred.A1.Zn1) 
  
  ### estimator of gamma
  gamma.maob.mod = lm(gamma.pseudo.maob ~ updrs_1 +  updrs_2 + updrs_3 + time_from_diagnosis.2nd + levodopa_dose.2nd +
                        gender + age.2nd  + bmi + Time_since_first_Levo , data = dat.bt, weights = gamma.weight)
  dat.bt$gamma.maob.pred = predict(gamma.maob.mod, newdata = dat.bt)
  
  ########### dra #######################################
  ## gamma prime 
  g.p.dra.mod = lm(g.p.DRA ~ Z+  updrs_1 +  updrs_2 + updrs_3 + time_from_diagnosis.2nd + levodopa_dose.2nd +
                     gender + age.2nd  + bmi  + Time_since_first_Levo , data = dat.bt)
  
  dat.bt$g.p.dra.pred = predict(g.p.dra.mod, newdata = temp.dat.opt)   
  
  dat.bt$gamma.pseudo.dra = (dat.bt$g.p.DRA - dat.bt$g.p.dra.pred)*2/(dat.bt$A - dat.bt$pred.A1.Zn1)
  ## gamma 
  gamma.dra.mod = lm(gamma.pseudo.dra ~ updrs_1 +  updrs_2 + updrs_3 + time_from_diagnosis.2nd + levodopa_dose.2nd +
                       gender + age.2nd  + bmi + Time_since_first_Levo , data = dat.bt, weights = gamma.weight)
  dat.bt$gamma.dra.pred = predict(gamma.dra.mod, newdata = dat.bt)
  
  
  #calculate the value function
  value.opt = with(dat.bt[!dat.bt$train,], mean((outcome*I(fitted.owl == A)*Z*A*delta.marg)/( pmax(fz.pred*delta.pred,lb) )))
  value.maob = with(dat.bt[!dat.bt$train,], mean((outcome*I(fitted.maob == A)*Z*A*delta.marg)/ pmax(fz.pred*delta.pred,lb) ))
  value.dra = with(dat.bt[!dat.bt$train,], mean((outcome*I(fitted.dra == A)*Z*A*delta.marg)/pmax(fz.pred*delta.pred,lb)))
  value.behavior = mean(dat.bt$outcome[!dat.bt$train]) 
  
  
  ## calculate the mr value function
  ##### opt
  
  component1.opt = with(dat.bt[!dat.bt$train,], (outcome*I(fitted.owl == A)*Z*A*delta.marg)/pmax(fz.pred*delta.pred,lb))
  component2.opt = with(dat.bt[!dat.bt$train,],  (Z*g.p.opt.pred*delta.marg)/pmax(fz.pred*delta.pred,lb))
  component3.opt = dat.bt$gamma.opt.pred[!dat.bt$train]
  component4.opt = with(dat.bt[!dat.bt$train,], Z*(A - pred.A)*gamma.opt.pred*delta.marg/(2*pmax(fz.pred*delta.pred,lb) ))
  
  value.opt.mr = mean(component1.opt - component2.opt + component3.opt - component4.opt) 
  
  ####### maob
  component1.maob = with(dat.bt[!dat.bt$train,], (outcome*I(fitted.maob == A)*Z*A*delta.marg)/pmax(fz.pred*delta.pred,lb))
  component2.maob = with(dat.bt[!dat.bt$train,],  (Z*g.p.maob.pred*delta.marg)/pmax(fz.pred*delta.pred,lb))
  component3.maob = dat.bt$gamma.maob.pred[!dat.bt$train]
  component4.maob = with(dat.bt[!dat.bt$train,], Z*(A - pred.A)*gamma.maob.pred*delta.marg/(2*pmax(fz.pred*delta.pred,lb) ))
  
  value.maob.mr = mean(component1.maob - component2.maob + component3.maob - component4.maob)
  
  ####### dra
  component1.dra = with(dat.bt[!dat.bt$train,], (outcome*I(fitted.dra == A)*Z*A*delta.marg)/pmax(fz.pred*delta.pred,lb))
  component2.dra = with(dat.bt[!dat.bt$train,],  (Z*g.p.dra.pred*delta.marg)/pmax(fz.pred*delta.pred,lb))
  component3.dra = dat.bt$gamma.dra.pred[!dat.bt$train]
  component4.dra = with(dat.bt[!dat.bt$train,], Z*(A - pred.A)*gamma.dra.pred*delta.marg/(2*pmax(fz.pred*delta.pred,lb) ))
  
  value.dra.mr = mean(component1.dra -  component2.dra + component3.dra - component4.dra)
  
  #variance of the mr
  est.var.opt = var((component1.opt - component2.opt + component3.opt - component4.opt))*(1/length(dat.bt$outcome[!dat.bt$train]))
  est.var.dra = var((component1.dra -  component2.dra + component3.dra - component4.dra))*(1/length(dat.bt$outcome[!dat.bt$train]))
  est.var.maob = var((component1.maob - component2.maob + component3.maob - component4.maob))*(1/length(dat.bt$outcome[!dat.bt$train]))
  
}

#unit test

fit.iv.owl(itr_formula = ~ updrs_1 + updrs_2 + updrs_3 + time_from_diagnosis.2nd + levodopa_dose.2nd + gender + age.2nd + bmi + dys + fluc + nausea + sleep + Time_since_first_Levo,
           iv_formula =  ~ updrs_1 + updrs_2 + updrs_3 + time_from_diagnosis + levodopa_dose + gender + age + bmi +Time_since_first_Levo, 
           trt_formula = ~ Z + updrs_1 + updrs_2 + updrs_3 + time_from_diagnosis.2nd + levodopa_dose.2nd + gender + age.2nd + bmi + Time_since_first_Levo,
           dat = iv.dat2a, trt ="A", iv = "Z", outcome = "outcome_updrs3", seed = 3431, minimize = T, center.outcome = T )

fit.iv.owl(itr_formula = ~ updrs_1 ,
           iv_formula =  ~ updrs_1 + updrs_2 + updrs_3 + time_from_diagnosis + age + Time_since_first_Levo , 
           trt_formula = ~ Z ,
           dat = iv.dat2a, trt ="A", iv = "Z", outcome = "outcome_updrs3", seed = 3431 )
