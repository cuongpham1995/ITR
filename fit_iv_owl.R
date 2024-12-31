fit.iv.owl = function(itr_formula, iv_formula, trt_formula, mean_formula = NULL, method, trt, iv, outcome, dat, 
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
  
  trt_formula2 <- update(trt_formula, paste("~ . +", "Z")) #adding the instrument variable
  trt.logistic.formula = paste( "as.factor(A)" , deparse(trt_formula2, width.cutoff = 500), collapse = "") 
  
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
  
  if(method == "iv owl"){ 
    dat$weight = with(dat, Y*A*Z/(Z.pred*delta.pred))
  }else if(method == "owl"){
    
    #if the owl method is used, we do not use the instrument variable 
    trt.logistic.formula = paste( "as.factor(A)" , deparse(trt_formula, width.cutoff = 500), collapse = "") 
  
    trt.logistic = glm(trt.logistic.formula, data = dat, family = binomial(link = "logit"))
    
    pred.A = predict(trt.logistic, type = "response", newdata = dat)
    dat$A.pred.2 = ifelse(dat$A == 1, pred.A, 1 - pred.A)
    
    dat$weight = with(dat, Y/A.pred.2)
  }else {stop("Invalid method. Enter: iv owl or owl")}
  
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
  
  
  
  dat$fitted.owl = as.numeric(as.character(fitted.owl)) #opt treatment regimes
  dat$fitted.A = rep(-1, length(fitted.owl)) #always -1 strategy. A is -1 
  dat$fitted.B = rep(1, length(fitted.owl)) #always 1 strategy. B is 1 
  
  

  
  owl.coef = svm.coef(mod.owl, type = "line")
  names(owl.coef) = c("Intercept", itr_vars)
  

  
  #pseudo outcome
  dat$g.p.opt = dat$Y*dat$A*I(dat$fitted.owl == dat$A)
  dat$g.p.A = dat$Y*dat$A*I(dat$fitted.A == dat$A)
  dat$g.p.B = dat$Y*dat$A*I(dat$fitted.B == dat$A)
  
  trt.logistic.formula = paste( "as.factor(A)" , deparse(trt_formula, width.cutoff = 500), collapse = "")
  
  ### parametric estimator of gamma
  
  ############ opt ################
  ### estimator of gamma prime
  mean_formula2 <- update(mean_formula, paste("~ . +", "Z"))
  gamma.prime.opt.formula = paste("g.p.opt", deparse(mean_formula2, width.cutoff = 500), collapse = "")
  gamma.opt.formula = paste("gamma.pseudo.opt", deparse(mean_formula, width.cutoff = 500), collapse = "")
  
  g.p.opt.mod = lm(gamma.prime.opt.formula, data = dat)

  dat_mod = dat
  dat_mod$Z = - 1
  dat$g.p.opt.pred = predict(g.p.opt.mod, newdata = dat_mod)   
  
  
  #### estimator of gamma
  
  dat$gamma.pseudo.opt = (dat$g.p.opt - dat$g.p.opt.pred)*2/(dat$A - dat$A.Zn1) 
  dat$gamma.weight = ((dat$A - dat$A.Zn1)*dat$Z/(2*dat$Z.pred))^2
  
  gamma.opt.mod = lm( gamma.opt.formula, data = dat, weights = gamma.weight)
  dat$gamma.opt.pred = predict(gamma.opt.mod, newdata = dat)
  
  
  #### maob ############################
  ### estimator of gamma prime
  gamma.prime.A.formula = paste("g.p.A", deparse(mean_formula2, width.cutoff = 500), collapse = "")
  gamma.A.formula = paste("gamma.pseudo.A", deparse(mean_formula, width.cutoff = 500), collapse = "")
  
  g.p.A.mod = lm(gamma.prime.A.formula, data = dat)
  
  dat$g.p.A.pred = predict(g.p.A.mod, newdata = dat_mod)   
  
  dat$gamma.pseudo.A = (dat$g.p.A - dat$g.p.A.pred)*2/(dat$A - dat$A.Zn1) 
  
  ### estimator of gamma
  gamma.A.mod = lm(gamma.A.formula, data = dat, weights = gamma.weight)
  dat$gamma.A.pred = predict(gamma.A.mod, newdata = dat)
  
  
  ########### dra #######################################
  gamma.prime.B.formula = paste("g.p.B", deparse(mean_formula2, width.cutoff = 500), collapse = "")
  gamma.B.formula = paste("gamma.pseudo.B", deparse(mean_formula, width.cutoff = 500), collapse = "")
  
  g.p.B.mod = lm(gamma.prime.B.formula, data = dat)
  
  dat$g.p.B.pred = predict(g.p.B.mod, newdata = dat_mod)   
  
  dat$gamma.pseudo.B = (dat$g.p.B - dat$g.p.B.pred)*2/(dat$A - dat$A.Zn1) 
  
  ### estimator of gamma
  gamma.B.mod = lm(gamma.B.formula, data = dat, weights = gamma.weight)
  dat$gamma.B.pred = predict(gamma.B.mod, newdata = dat)
  
  
  #calculate the value function
  value.opt = with(dat[!dat$train,], mean((Y*I(fitted.owl == A)*Z*A)/( Z.pred*delta.pred )))
  value.A = with(dat[!dat$train,], mean((Y*I(fitted.A == A)*Z*A)/( Z.pred*delta.pred )))
  value.B = with(dat[!dat$train,], mean((Y*I(fitted.B == A)*Z*A)/( Z.pred*delta.pred )))
  value.behavior = mean(dat$Y[!dat$train]) 
  
  
  ## calculate the mr value function
  ##### opt
  
  component1.opt = with(dat[!dat$train,], (Y*I(fitted.owl == A)*Z*A)/ (Z.pred*delta.pred))
  component2.opt = with(dat[!dat$train,], (Z*g.p.opt.pred)/(Z.pred*delta.pred))
  component3.opt = dat$gamma.opt.pred[!dat$train]
  component4.opt = with(dat[!dat$train,], Z*(A - A.pred)*gamma.opt.pred/(2*(Z.pred*delta.pred)))
  
  value.opt.mr = mean(component1.opt - component2.opt + component3.opt - component4.opt) 
  
  #### trt A
  component1.A = with(dat[!dat$train,], (Y*I(fitted.A == A)*Z*A)/ (Z.pred*delta.pred))
  component2.A = with(dat[!dat$train,], (Z*g.p.A.pred)/(Z.pred*delta.pred))
  component3.A = dat$gamma.A.pred[!dat$train]
  component4.A = with(dat[!dat$train,], Z*(A - A.pred)*gamma.A.pred/(2*(Z.pred*delta.pred)))
  
  value.A.mr = mean(component1.A - component2.A + component3.A - component4.A) 
  #### trt B
  component1.B = with(dat[!dat$train,], (Y*I(fitted.B == A)*Z*A)/ (Z.pred*delta.pred))
  component2.B = with(dat[!dat$train,], (Z*g.p.B.pred)/(Z.pred*delta.pred))
  component3.B = dat$gamma.B.pred[!dat$train]
  component4.B = with(dat[!dat$train,], Z*(A - A.pred)*gamma.B.pred/(2*(Z.pred*delta.pred)))
  
  value.B.mr = mean(component1.B - component2.B + component3.B - component4.B) 
  
 
  #variance of the mr
  var.opt.mr = var((component1.opt - component2.opt + component3.opt - component4.opt))/sum( 1 - dat$train )
  var.B.mr= var((component1.B -  component2.B + component3.B - component4.B))/sum( 1 - dat$train )
  var.A.mr = var((component1.A - component2.A + component3.A - component4.A))/sum( 1 - dat$train )
  
  
  #return the results 
  value.function =data.frame(est.ipw = c(value.opt, value.A, value.B), 
                             est.mr = c(value.opt.mr, value.A.mr, value.B.mr))
  colnames(value.function) = c("IPW estimator", "MR estimator")
  rownames(value.function) = c("Opt Trt.", "Trt A", "Trt B")
  
  var.value.function = c(var.opt.mr, var.A.mr, var.B.mr)
  names(var.value.function) = c("Var Opt.", "Var Trt A.", "Var Trt B.")
  
  
  return(list(model = mod.owl, 
              value.funcion.IPW =c(value.opt, value.A, value.B),
              value.funcion.MR = c(value.opt.mr, value.A.mr, value.B.mr),
              variance = var.value.function, 
              coef.owl = owl.coef) )
}

#unit test

#res = fit.iv.owl(itr_formula = ~ updrs_1 + updrs_2 + updrs_3 + time_from_diagnosis.2nd + levodopa_dose.2nd + gender + age.2nd + bmi + dys + fluc + nausea + sleep + Time_since_first_Levo,
#           iv_formula =  ~ updrs_1 + updrs_2 + updrs_3 + time_from_diagnosis + levodopa_dose + gender + age + bmi +Time_since_first_Levo, 
#           trt_formula = ~ updrs_1 + updrs_2 + updrs_3 + time_from_diagnosis.2nd + levodopa_dose.2nd + gender + age.2nd + bmi + Time_since_first_Levo,
#           mean_formula = ~ updrs_1 +  updrs_2 + updrs_3 + time_from_diagnosis.2nd + levodopa_dose.2nd +
#             gender + age.2nd  + bmi  + Time_since_first_Levo, method = "iv owl",
#           dat = iv.dat2a, trt ="A", iv = "Z", outcome = "outcome_updrs3", seed = 3431, minimize = T, center.outcome = T )
#res$value.funcion
