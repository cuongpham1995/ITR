no.obs = function(x){
  return(length(x[!is.na(x)]))
}

expit <- function(x) {
  1 / (1 + exp(-x))
}

svm.coef = function(mod, type = c("line", "slope")){
  type = match.arg(type)
  #the linear decision is in the form b0 = b1*x1 + b2*x2 + ... + bn*xn
  coef = t(mod$coefs[,1]) %*% mod$SV
  intercept = mod$rho[1]*(-1)
  beta = c(intercept, coef)
  
  if(type == "line"){
    return(beta)
  }else{
    #the linear decision is in the form  x2 = b0/b2  - b1/b2 *x1 - b3/b2 *x3 - ... - bn/b2 *xn
    slop = -beta/beta[3]
    slop[1] = -slop[1]
    return(slop[-3])
  }
  
}


#this is the case where we assume there are no unmeasured confounders
fit.owl.bt.ls1.postcept = function(seed, dat, covariates, kernel.type, degree.poly,  outcome.col, minimize = T, complier = F, center.outcome = T, scale.covariates = T, lb = 0, logistic){
  tryCatch({ 
    set.seed(seed)
    require(WeightSVM)
    
    # dat = iv.dat2a
    #  kernel.type = "linear"
    #  outcome.col = "outcome_updrs1"
    #  complier = F 
    #  minimize = T 
    #  lb = 0
    #  remove.outlier = T
    #  center.outcome = T
    #  scale.covariates = T
    #  covariates = c("updrs_1", "updrs_2" , "updrs_3" , "time_from_diagnosis.2nd" , "levodopa_dose.2nd","gender" , "age.2nd" ,"bmi" , "dys")
    # 
    if(any(covariates %in% c("dys", "fluc", "nausea", "sleep", "ortho")) ){
      dat = dat %>% filter(!is.na(dys))
    }
    
    #prepare data
    if(minimize){
      dat$outcome = dat[,outcome.col]*-1
    }else{dat$outcome = dat[,outcome.col]}
    
    
    if(center.outcome){
      dat$outcome = dat$outcome - mean(dat$outcome)
    }
    ############calculate nuisance parameters
    #propensity model
    
    if(scale.covariates){
      dat[,c("updrs_1", "updrs_2", "updrs_3", "time_from_diagnosis", "time_from_diagnosis.2nd", "levodopa_dose","levodopa_dose.2nd",
             "age", "age.2nd", "bmi", "Time_since_first_Levo")] = apply( dat[,c("updrs_1", "updrs_2", "updrs_3", "time_from_diagnosis", "time_from_diagnosis.2nd", "levodopa_dose","levodopa_dose.2nd", "age", "age.2nd", "bmi", "Time_since_first_Levo")],2,scale)
    }
    
    dat.bt = dat
    dat.bt$train = sample(c(T,F), nrow(dat.bt), prob = c(0.7,0.3), replace = T )
    
    logistic.A = glm(as.factor(A) ~ Z +  updrs_1 + updrs_2 + updrs_3 + time_from_diagnosis.2nd + levodopa_dose.2nd + gender + age.2nd + bmi + Time_since_first_Levo , data = dat.bt, family = binomial(link = "logit"))
    
    pred.A1 = predict(logistic.A, type = "response", newdata = dat.bt)
    dat.bt$pred.A = ifelse(dat.bt$A == 1, pred.A1, 1 - pred.A1)
    
    logistic.A.v2 = glm(as.factor(A) ~  updrs_1 + updrs_2 + updrs_3 + time_from_diagnosis.2nd + levodopa_dose.2nd + gender + age.2nd + bmi + Time_since_first_Levo , data = dat.bt, family = binomial(link = "logit"))
    
    pred.A1.v2 = predict(logistic.A.v2, type = "response", newdata = dat.bt)
    dat.bt$pred.A.v2 = ifelse(dat.bt$A == 1, pred.A1.v2, 1 - pred.A1.v2)
    
    dat.Z1 = dat.bt %>% mutate(Z = 1)
    dat.bt$pred.A1.Z1 = predict(logistic.A, newdata = dat.Z1, type = "response")
    
    dat.Zn1 = dat.bt %>% mutate(Z = -1)
    dat.bt$pred.A1.Zn1 = predict(logistic.A, newdata = dat.Zn1, type = "response")
    
    
    logistic.Z = glm(as.factor(Z) ~ updrs_1 + updrs_2 + updrs_3 + time_from_diagnosis + levodopa_dose + gender + age + bmi + Time_since_first_Levo , data = dat.bt, family = binomial(link = "logit"))
    
    pred.z = predict(logistic.Z, newdata = dat.bt, type = "response") 
    dat.bt$fz.pred = ifelse(dat.bt$Z == 1, pred.z , 1 - pred.z) 
    
    if(complier){
      dat.bt$delta.pred = rep(sum(dat.bt$A == 1 & dat.bt$Z == 1)/sum(dat.bt$Z == 1) - sum(dat.bt$A == 1 & dat.bt$Z == -1)/sum(dat.bt$Z == -1), nrow(dat.bt))
      weight = dat.bt$outcome[dat.bt$train]*dat.bt$A[dat.bt$train]*dat.bt$Z[dat.bt$train]/(dat.bt$fz.pred[dat.bt$train]) #OWL weight
    }else{
      dat.bt$delta.pred = dat.bt$pred.A1.Z1 - dat.bt$pred.A1.Zn1
      dat.bt$weight = with(dat.bt, outcome/( pmax(pred.A.v2,lb) ) )
    }
    
    #uncomment to remove the outliers by weight.  
    # outlier.cov = dat.bt %>% filter(weight <= quantile(weight, prob = 0.05) | weight >= quantile(weight, prob = 0.95) )
    # outlier.id = outlier.cov$patno %>% unique()
    # outlier.weight = outlier.cov$weight
    # outlier.cov = apply(outlier.cov[, c(covariates, "outcome")],2,mean)
    # 
    
    #dat.bt = dat.bt %>% filter(weight > quantile(weight, prob = 0.05) & weight < quantile(weight, prob = 0.95) )
    # non.outlier.weight = dat.bt$weight
    dat.bt$fitted.owl = NA
    
    #   included.cov = apply(dat.bt[,c(covariates, "outcome")],2,mean)
    
    ## calculate the weight
    
    dat.bt$lab = sign(dat.bt$weight)*dat.bt$A
    
    formula_str <- paste("as.factor(lab)", "~", paste(covariates, collapse = " + "))
    formula.owl <- as.formula(formula_str)
    
    mod.owl = wsvm(formula.owl
                   , data = dat.bt[dat.bt$train,], degree = degree.poly , weight = abs(dat.bt$weight[dat.bt$train]), kernel = kernel.type, 
                   cross = 10, scale = T) #OWL method
    
  
    ## MODEL FITTING
    
    fitted.owl =  predict(mod.owl, newdata = dat.bt[, covariates])# optTx(fitOWL, dat[dat$train == 0,])$optimalTx
    fitted.owl = as.numeric(as.character(fitted.owl))
    
    dat.bt$fitted.owl = fitted.owl
    
    owl.coef = svm.coef(mod.owl, type = "line")
    names(owl.coef) = c("Intercept", covariates)
    
    
    MAOB.same = apply(dat.bt[dat.bt$A == -1 & dat.bt$A == dat.bt$fitted.owl & dat.bt$train, covariates],2,mean, na.rm =T)
    MAOB.diff = apply(dat.bt[dat.bt$A == -1 & dat.bt$A != dat.bt$fitted.owl & dat.bt$train, covariates],2,mean, na.rm = T)
    
    DRA.same = apply(dat.bt[dat.bt$A == 1 & dat.bt$A == dat.bt$fitted.owl & dat.bt$train, covariates],2,mean, na.rm = T)
    DRA.diff = apply(dat.bt[dat.bt$A == 1 & dat.bt$A != dat.bt$fitted.owl & dat.bt$train, covariates],2,mean, na.rm = T)
    
    
    dat.bt$fitted.maob = rep(-1, nrow(dat.bt))
    dat.bt$fitted.dra = rep(1, nrow(dat.bt))
    
    
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
                         gender + age.2nd  + bmi  + Time_since_first_Levo , data = dat.bt, weights = gamma.weight)
    dat.bt$gamma.opt.pred = predict(gamma.opt.mod, newdata = dat.bt)
    
    #### maob ############################
    ### estimator of gamma prime
    
    g.p.maob.mod = lm(g.p.MAOB ~ Z+  updrs_1 +  updrs_2 + updrs_3 + time_from_diagnosis.2nd + levodopa_dose.2nd +
                        gender + age.2nd  + bmi + Time_since_first_Levo  , data = dat.bt)
    
    dat.bt$g.p.maob.pred = predict(g.p.maob.mod, newdata = temp.dat.opt)   
    
    dat.bt$gamma.pseudo.maob = (dat.bt$g.p.MAOB - dat.bt$g.p.maob.pred)*2/(dat.bt$A - dat.bt$pred.A1.Zn1) 
    
    ### estimator of gamma
    gamma.maob.mod = lm(gamma.pseudo.maob ~ updrs_1 +  updrs_2 + updrs_3 + time_from_diagnosis.2nd + levodopa_dose.2nd +
                          gender + age.2nd  + bmi + Time_since_first_Levo  , data = dat.bt, weights = gamma.weight)
    dat.bt$gamma.maob.pred = predict(gamma.maob.mod, newdata = dat.bt)
    
    ########### dra #######################################
    ## gamma prime 
    g.p.dra.mod = lm(g.p.DRA ~ Z+  updrs_1 +  updrs_2 + updrs_3 + time_from_diagnosis.2nd + levodopa_dose.2nd +
                       gender + age.2nd  + bmi + Time_since_first_Levo  , data = dat.bt)
    
    dat.bt$g.p.dra.pred = predict(g.p.dra.mod, newdata = temp.dat.opt)   
    
    dat.bt$gamma.pseudo.dra = (dat.bt$g.p.DRA - dat.bt$g.p.dra.pred)*2/(dat.bt$A - dat.bt$pred.A1.Zn1)
    ## gamma 
    gamma.dra.mod = lm(gamma.pseudo.dra ~ updrs_1 +  updrs_2 + updrs_3 + time_from_diagnosis.2nd + levodopa_dose.2nd +
                         gender + age.2nd  + bmi , data = dat.bt, weights = gamma.weight)
    dat.bt$gamma.dra.pred = predict(gamma.dra.mod, newdata = dat.bt)
    
    
    #calculate the value function
    value.opt = with(dat.bt[!dat.bt$train,], mean((outcome*I(fitted.owl == A)*Z*A)/( pmax(fz.pred*delta.pred,lb) )))
    value.maob = with(dat.bt[!dat.bt$train,], mean((outcome*I(fitted.maob == A)*Z*A)/ pmax(fz.pred*delta.pred,lb) ))
    value.dra = with(dat.bt[!dat.bt$train,], mean((outcome*I(fitted.dra == A)*Z*A)/pmax(fz.pred*delta.pred,lb)))
    value.behavior = mean(dat.bt$outcome[!dat.bt$train]) 
    
    
    ## calculate the mr value function
    ##### opt
    
    component1.opt = with(dat.bt[!dat.bt$train,], (outcome*I(fitted.owl == A)*Z*A)/pmax(fz.pred*delta.pred,lb))
    component2.opt = with(dat.bt[!dat.bt$train,],  (Z*g.p.opt.pred)/pmax(fz.pred*delta.pred,lb))
    component3.opt = dat.bt$gamma.opt.pred[!dat.bt$train]
    component4.opt = with(dat.bt[!dat.bt$train,], Z*(A - pred.A)*gamma.opt.pred/(2*pmax(fz.pred*delta.pred,lb) ))
    
    value.opt.mr = mean(component1.opt - component2.opt + component3.opt - component4.opt) 
    
    ####### maob
    component1.maob = with(dat.bt[!dat.bt$train,], (outcome*I(fitted.maob == A)*Z*A)/pmax(fz.pred*delta.pred,lb))
    component2.maob = with(dat.bt[!dat.bt$train,],  (Z*g.p.maob.pred)/pmax(fz.pred*delta.pred,lb))
    component3.maob = dat.bt$gamma.maob.pred[!dat.bt$train]
    component4.maob = with(dat.bt[!dat.bt$train,], Z*(A - pred.A)*gamma.maob.pred/(2*pmax(fz.pred*delta.pred,lb) ))
    
    value.maob.mr = mean(component1.maob - component2.maob + component3.maob - component4.maob)
    
    ####### dra
    component1.dra = with(dat.bt[!dat.bt$train,], (outcome*I(fitted.dra == A)*Z*A)/pmax(fz.pred*delta.pred,lb))
    component2.dra = with(dat.bt[!dat.bt$train,],  (Z*g.p.dra.pred)/pmax(fz.pred*delta.pred,lb))
    component3.dra = dat.bt$gamma.dra.pred[!dat.bt$train]
    component4.dra = with(dat.bt[!dat.bt$train,], Z*(A - pred.A)*gamma.dra.pred/(2*pmax(fz.pred*delta.pred,lb) ))
    
    value.dra.mr = mean(component1.dra -  component2.dra + component3.dra - component4.dra)
    
    #variance of the mr
    est.var.opt = var((component1.opt - component2.opt + component3.opt - component4.opt))*(1/length(dat.bt$outcome[!dat.bt$train]))
    est.var.dra = var((component1.dra -  component2.dra + component3.dra - component4.dra))*(1/length(dat.bt$outcome[!dat.bt$train]))
    est.var.maob = var((component1.maob - component2.maob + component3.maob - component4.maob))*(1/length(dat.bt$outcome[!dat.bt$train]))
    
    
    
    ########################################################################
    #################### Return results ####################################   
    
    value.ipw = c( value.opt, value.maob, value.dra, value.behavior)
    names(value.ipw) = c("Optimal Value Function", "MAOB Value Function", "DRA Value Function", "Behavior Policy")
    value.mr = c(value.opt.mr, value.maob.mr, value.dra.mr)
    names(value.mr) = c("Optimal Value Function MR", "MAOB Value Function MR", "DRA Value Function MR")
    var.value.mr = c( est.var.opt, est.var.maob, est.var.dra )
    names(var.value.mr) = c("Var Optimal Value Function MR", "Var MAOB Value Function MR", "Var DRA Value Function MR")
    
    prop.dra = mean(fitted.owl == 1)
    
    return(list(prop.dra = mean(prop.dra), 
                value.func.ipw = value.ipw, 
                value.func.mr = value.mr, 
                var.value.func.mr = var.value.mr,  
                coef.beta = owl.coef, 
                MAOB.same = MAOB.same, 
                MAOB.diff = MAOB.diff, 
                DRA.same = DRA.same, 
                DRA.diff = DRA.diff, 
                weight.prop = with(dat.bt[!dat.bt$train,], pmax(fz.pred*delta.pred,lb))))
  }, error = function(error_message){
    message(error_message)
    return(NA)
  })
  
} 



fit.owl.mr.sim = function(seed, dat, covariates, kernel.type, degree.poly,  outcome.col, minimize = T, complier = F, center.outcome = T, scale.covariates = T, lb = 0, logistic ){
  tryCatch({ 
    set.seed(seed)
    require(WeightSVM) #to fit weighted svm
    
    #dat = iv.dat2a %>% filter(Time_baseline_2nd_line <= 90)
    #covariates = sel.cov
    #outcome.col = "outcome_updrs3"
    #minimize = T
    #center.outcome = T
    #scale.covariates = F
    #degree.poly = 1
    #kernel.type = "linear"
    
    if(any(covariates %in% c("dys", "fluc", "nausea", "sleep", "ortho")) ){
      dat = dat %>% filter(!is.na(dys))
    }
    
    #prepare data
    if(minimize){
      dat$outcome = dat[,outcome.col]*-1
    }else{dat$outcome = dat[,outcome.col]}
    
    
    if(center.outcome){
      dat$outcome = dat$outcome - mean(dat$outcome)
    }
    ############calculate nuisance parameters
    #propensity model
    
    if(scale.covariates){
      dat[,c("updrs_1", "updrs_2", "updrs_3", "time_from_diagnosis", "time_from_diagnosis.2nd", "levodopa_dose","levodopa_dose.2nd",
             "age", "age.2nd", "bmi", "Time_since_first_Levo")] = apply( dat[,c("updrs_1", "updrs_2", "updrs_3", "time_from_diagnosis", "time_from_diagnosis.2nd", "levodopa_dose","levodopa_dose.2nd", "age", "age.2nd", "bmi", "Time_since_first_Levo")],2,scale)
    }
    
    dat.bt = dat
    dat.bt$train = sample(c(T,F), nrow(dat.bt), prob = c(0.7,0.3), replace = T )
    
    
    #prepare data
    
    
    ############calculate nuisance parameters
    #propensity model
    logistic.A = glm(as.factor(A) ~ updrs_1 + updrs_2 +
                       updrs_3 + time_from_diagnosis.2nd + levodopa_dose.2nd + gender + age.2nd + bmi + Time_since_first_Levo , data = dat.bt, family = binomial(link = "logit"))
    
    pred.a = predict(logistic.A, newdata = dat.bt, type = "response") 
    dat.bt$fa.pred = ifelse(dat.bt$A == 1, pred.a , 1 - pred.a) 
    
    #mean outcome model
    outcome.mod = lm(outcome ~ A + updrs_1 + updrs_2 + updrs_3 + time_from_diagnosis.2nd + 
                        levodopa_dose.2nd + gender + age.2nd + bmi + Time_since_first_Levo, data = dat.bt)
    
    dat.bt.A1 = dat.bt %>% mutate(A = 1)
    dat.bt$outcome.A1 = predict(outcome.mod, newdata = dat.bt.A1)
    dat.bt.An1 = dat.bt %>% mutate(A = -1)
    dat.bt$outcome.An1 = predict(outcome.mod, newdata =  dat.bt.An1)
    dat.bt$outcome.pred = predict(outcome.mod)
    
    # weight = (dat$Y[dat$train == 1] - train.outcome)/fz.train + dat$A[dat$train == 1]*(train.outcome.A1 - train.outcome.An1) #MR OWL weight
  
     dat.bt$weight = with(dat.bt, outcome/fa.pred) #OWL weight
    
     dat.bt$lab = with(dat.bt,sign(weight)*A)
    
    
    #models building
     formula_str <- paste("as.factor(lab)", "~", paste(covariates, collapse = " + "))
     formula.owl <- as.formula(formula_str)
    

    mod.owl = wsvm(formula.owl
                   , data = dat.bt[dat.bt$train,], degree = degree.poly , weight = abs(dat.bt$weight[dat.bt$train]), kernel = kernel.type, 
                   cross = 10, scale = T) #OWL method
    
    ## MODEL FITTING
    
    fitted.owl =  predict(mod.owl, newdata = dat.bt[, covariates])# optTx(fitOWL, dat[dat$train == 0,])$optimalTx
    fitted.owl = as.numeric(as.character(fitted.owl))
    
    dat.bt$fitted.owl = fitted.owl
    dat.bt$fitted.maob = rep(-1, nrow(dat.bt))
    dat.bt$fitted.dra = rep(1, nrow(dat.bt))
    
    owl.coef = svm.coef(mod.owl, type = "line")
    names(owl.coef) = c("Intercept", covariates)
    
    
    MAOB.same = apply(dat.bt[dat.bt$A == -1 & dat.bt$A == dat.bt$fitted.owl & dat.bt$train, covariates],2,mean, na.rm =T)
    MAOB.diff = apply(dat.bt[dat.bt$A == -1 & dat.bt$A != dat.bt$fitted.owl & dat.bt$train, covariates],2,mean, na.rm = T)
    
    DRA.same = apply(dat.bt[dat.bt$A == 1 & dat.bt$A == dat.bt$fitted.owl & dat.bt$train, covariates],2,mean, na.rm = T)
    DRA.diff = apply(dat.bt[dat.bt$A == 1 & dat.bt$A != dat.bt$fitted.owl & dat.bt$train, covariates],2,mean, na.rm = T)
    
    dat.bt$fitted.maob = rep(-1, nrow(dat.bt))
    dat.bt$fitted.dra = rep(1, nrow(dat.bt))
    
    #estimated value function
    value.opt = mean( ((dat.bt$outcome*I(dat.bt$fitted.owl == dat.bt$A))/dat.bt$fa.pred)[!dat.bt$train])
    value.maob =  mean(((dat.bt$outcome*I(dat.bt$fitted.maob == dat.bt$A))/dat.bt$fa.pred)[!dat.bt$train])
    value.dra = mean( ((dat.bt$outcome*I(dat.bt$fitted.dra == dat.bt$A))/dat.bt$fa.pred)[!dat.bt$train])
    value.behavior =  mean(dat.bt$outcome[!dat.bt$train]) 
    
    value.opt.train = mean( ((dat.bt$outcome*I(dat.bt$fitted.owl == dat.bt$A))/dat.bt$fa.pred)[dat.bt$train])
    value.maob.train =  mean(((dat.bt$outcome*I(dat.bt$fitted.maob == dat.bt$A))/dat.bt$fa.pred)[dat.bt$train])
    value.dra.train = mean( ((dat.bt$outcome*I(dat.bt$fitted.dra == dat.bt$A))/dat.bt$fa.pred)[dat.bt$train])
    value.behavior.train =  mean(dat.bt$outcome[dat.bt$train]) 
  
    
    #print(c(value.opt.train, value.maob.train, value.dra.train))
    
    #value function mr
    value.opt.mr = value.opt - mean( with(dat.bt[!dat.bt$train,], (outcome.pred*I(fitted.owl == A ))/fa.pred 
                                     - I(fitted.owl == 1)*outcome.A1 - I(fitted.owl == -1)*outcome.An1)) 
    
    value.maob.mr = value.maob - mean( with(dat.bt[!dat.bt$train,], (outcome.pred*I(fitted.maob == A ))/fa.pred 
                                            - I(fitted.maob == 1)*outcome.A1 - I(fitted.maob == -1)*outcome.An1)) 
    
    value.dra.mr = value.dra -  mean( with(dat.bt[!dat.bt$train,], (outcome.pred*I(fitted.dra == A ))/fa.pred 
                                           - I(fitted.dra == 1)*outcome.A1 - I(fitted.dra == -1)*outcome.An1)) 
    
    #variance of the mr
    var.value.opt.mr = var(with(dat.bt[!dat.bt$train,], ((outcome*I(fitted.owl == A))/fa.pred) -
                              (outcome.pred*I(fitted.owl == A ))/fa.pred 
                            + I(fitted.owl == 1)*outcome.A1 + I(fitted.owl == -1)*outcome.An1))
    
    var.value.maob.mr = var(with(dat.bt[!dat.bt$train,], ((outcome*I(fitted.maob == A))/fa.pred) -
                                   (outcome.pred*I(fitted.maob == A ))/fa.pred 
                                 + I(fitted.maob == 1)*outcome.A1 + I(fitted.maob == -1)*outcome.An1))
    
    var.value.dra.mr = var(with(dat.bt[!dat.bt$train,], ((outcome*I(fitted.dra == A))/fa.pred) -
                                  (outcome.pred*I(fitted.dra == A ))/fa.pred 
                                + I(fitted.dra == 1)*outcome.A1 + I(fitted.dra == -1)*outcome.An1))
    

   
    value.ipw = c( value.opt, value.maob, value.dra, value.behavior)
    names(value.ipw) = c("Optimal Value Function", "MAOB Value Function", "DRA Value Function", "Behavior Policy")
    value.mr = c(value.opt.mr, value.maob.mr, value.dra.mr)
    names(value.mr) = c("Optimal Value Function MR", "MAOB Value Function MR", "DRA Value Function MR")
    var.value.mr = c( var.value.opt.mr, var.value.maob.mr, var.value.dra.mr )
    names(var.value.mr) = c("Var Optimal Value Function MR", "Var MAOB Value Function MR", "Var DRA Value Function MR")
    
    prop.dra = mean(fitted.owl == 1)

    return(list(prop.dra = mean(prop.dra), 
                value.func.ipw = value.ipw, 
                value.func.mr = value.mr, 
                var.value.func.mr = var.value.mr,  
                coef.beta = owl.coef, 
                MAOB.same = MAOB.same, 
                MAOB.diff = MAOB.diff, 
                DRA.same = DRA.same, 
                DRA.diff = DRA.diff, 
                weight.prop = with(dat.bt[!dat.bt$train,], fa.pred)))
  }, error = function(error_message){
    message(error_message)
    return(NA)
  })
  
}



fit.iv.owl.bt.ls1.postcept = function(seed, dat, covariates, kernel.type, degree.poly,  outcome.col, minimize = T, complier = F, center.outcome = T, scale.covariates = T, lb = 0, logistic){
  tryCatch({ 
    set.seed(seed)
    require(WeightSVM)
    require(glmnet)
    
       # dat = iv.dat2a
       #  kernel.type = "linear"
       #  outcome.col = "outcome_updrs3"
       #  complier = F
       #  minimize = T 
       #  lb = 0
       #  logistic = F
       #  center.outcome = T
       #  scale.covariates = T
       #  covariates = sel.cov # c("updrs_1", "updrs_2" , "updrs_3" , "time_from_diagnosis.2nd" , "levodopa_dose.2nd","gender" , "age.2nd" ,"bmi" , "dys")
    # 
    
    
    if(any(covariates %in% c("dys", "fluc", "nausea", "sleep", "ortho")) ){
      dat = dat %>% filter(!is.na(dys))
    }
    
    #prepare data
    if(minimize){
      dat$outcome = dat[,outcome.col]*-1
    }else{dat$outcome = dat[,outcome.col]}
    
    
    if(center.outcome){
      dat$outcome = dat$outcome - mean(dat$outcome)
    }
    ############calculate nuisance parameters
    #propensity model
    
    if(scale.covariates){
      dat[,c("updrs_1", "updrs_2", "updrs_3", "time_from_diagnosis", "time_from_diagnosis.2nd", "levodopa_dose","levodopa_dose.2nd",
             "age", "age.2nd", "bmi", "Time_since_first_Levo" )] = apply( dat[,c("updrs_1", "updrs_2", "updrs_3", "time_from_diagnosis", "time_from_diagnosis.2nd", "levodopa_dose","levodopa_dose.2nd", "age", "age.2nd", "bmi", "Time_since_first_Levo")],2,scale)
    }
    
    dat.bt = dat
    dat.bt$train = sample(c(T,F), nrow(dat.bt), prob = c(0.7,0.3), replace = T )
    
    if(logistic){
    logistic.A = glm(as.factor(A) ~ Z +  updrs_1 + updrs_2 + updrs_3 + time_from_diagnosis.2nd + levodopa_dose.2nd + gender + age.2nd + bmi + Time_since_first_Levo  , data = dat.bt, family = binomial(link = "logit"))
    
    pred.A1 = predict(logistic.A, type = "response", newdata = dat.bt)
    dat.bt$pred.A = ifelse(dat.bt$A == 1, pred.A1, 1 - pred.A1)
    
    dat.Z1 = dat.bt %>% mutate(Z = 1)
    dat.bt$pred.A1.Z1 = predict(logistic.A, newdata = dat.Z1, type = "response")
    
    dat.Zn1 = dat.bt %>% mutate(Z = -1)
    dat.bt$pred.A1.Zn1 = predict(logistic.A, newdata = dat.Zn1, type = "response")
    
    logistic.Z = glm(as.factor(Z) ~ updrs_1 + updrs_2 + updrs_3 + time_from_diagnosis + levodopa_dose + gender + age + bmi +Time_since_first_Levo  , data = dat.bt, family = binomial(link = "logit"))
    
    pred.z = predict(logistic.Z, newdata = dat.bt, type = "response") 
    dat.bt$fz.pred = ifelse(dat.bt$Z == 1, pred.z , 1 - pred.z) 
    }else{
    
    ############ lasso model ##############################
    dat.bt$folds = sample(c(1:5), size = nrow(dat.bt), replace = T)
    for(i in 1:5){
     
    # Prepare data for LASSO
    x_A <- model.matrix(A ~ Z + updrs_1 + updrs_2 + updrs_3 + time_from_diagnosis.2nd + levodopa_dose.2nd + gender + age.2nd + bmi + Time_since_first_Levo, data = dat.bt)[,-1] # Remove intercept
    y_A <- as.factor(dat.bt$A)
    
    # Fit LASSO model for A
    lasso.A <- cv.glmnet(x_A[dat.bt$folds != i,], y_A[dat.bt$folds != i], family = "binomial", alpha = 1)  # alpha=1 for LASSO
    
    # Predict probabilities for A
    pred.A1 <- predict(lasso.A, newx = x_A[dat.bt$folds == i,], type = "response", s = "lambda.min")  # Use lambda.min from cross-validation
    dat.bt$pred.A[dat.bt$folds == i] <- ifelse(dat.bt$A[dat.bt$folds == i] == 1, pred.A1, 1 - pred.A1)
    
    # Predict with Z = 1
    dat.Z1 <- dat.bt[dat.bt$folds == i,] %>% mutate(Z = 1)
    x_A_Z1 <- model.matrix(A ~ Z + updrs_1 + updrs_2 + updrs_3 + time_from_diagnosis.2nd + levodopa_dose.2nd + gender + age.2nd + bmi + Time_since_first_Levo , data = dat.Z1)[,-1]
    dat.bt$pred.A1.Z1[dat.bt$folds == i] <- predict(lasso.A, newx = x_A_Z1, type = "response", s = "lambda.min")
    
    # Predict with Z = -1
    dat.Zn1 <- dat.bt[dat.bt$folds == i,] %>% mutate(Z = -1)
    x_A_Zn1 <- model.matrix(A ~ Z + updrs_1 + updrs_2 + updrs_3 + time_from_diagnosis.2nd + levodopa_dose.2nd + gender + age.2nd + bmi + Time_since_first_Levo , data = dat.Zn1)[,-1]
    dat.bt$pred.A1.Zn1[dat.bt$folds == i] <- predict(lasso.A, newx = x_A_Zn1, type = "response", s = "lambda.min")
    
    # LASSO model for Z
    x_Z <- model.matrix(Z ~ updrs_1 + updrs_2 + updrs_3 + time_from_diagnosis + levodopa_dose + gender + age + bmi + Time_since_first_Levo , data = dat.bt)[,-1]
    y_Z <- as.factor(dat.bt$Z)
    
    # Fit LASSO model for Z
    lasso.Z <- cv.glmnet(x_Z[dat.bt$folds != i,], y_Z[dat.bt$folds != i], family = "binomial", alpha = 1)
    
    # Predict probabilities for Z
    pred.z <- predict(lasso.Z, newx = x_Z[dat.bt$folds == i,], type = "response", s = "lambda.min")
    dat.bt$fz.pred[dat.bt$folds == i] = ifelse(dat.bt$Z[dat.bt$folds == i] == 1, pred.z, 1 - pred.z)
    
    }
    }
    #######################################################
    
    #try stablized weight
    dat.bt$delta.marg = 1 # sum(dat.bt$A == 1 & dat.bt$Z == 1)/sum(dat.bt$Z == 1) - sum(dat.bt$A == 1 & dat.bt$Z == -1)/sum(dat.bt$Z == -1) #ifelse(dat.bt$A == 1, mean(dat.bt$A == 1), 1 - mean(dat.bt$A == 1))
    
    if(complier){
      dat.bt$delta.pred = rep(sum(dat.bt$A == 1 & dat.bt$Z == 1)/sum(dat.bt$Z == 1) - sum(dat.bt$A == 1 & dat.bt$Z == -1)/sum(dat.bt$Z == -1), nrow(dat.bt))
      weight = dat.bt$outcome[dat.bt$train]*dat.bt$A[dat.bt$train]*dat.bt$Z[dat.bt$train]/(dat.bt$fz.pred[dat.bt$train]) #OWL weight
    }else{
      dat.bt$delta.pred = dat.bt$pred.A1.Z1 - dat.bt$pred.A1.Zn1
      dat.bt$weight = with(dat.bt, outcome*A*Z*delta.marg/( pmax(fz.pred*delta.pred,lb) ) )
    }
    
    #uncomment to remove the outliers by weight.  
    # outlier.cov = dat.bt %>% filter(weight <= quantile(weight, prob = 0.05) | weight >= quantile(weight, prob = 0.95) )
    # outlier.id = outlier.cov$patno %>% unique()
    # outlier.weight = outlier.cov$weight
    # outlier.cov = apply(outlier.cov[, c(covariates, "outcome")],2,mean)
    # 
    
    #dat.bt = dat.bt %>% filter(weight > quantile(weight, prob = 0.05) & weight < quantile(weight, prob = 0.95) )
    # non.outlier.weight = dat.bt$weight
    dat.bt$fitted.owl = NA
    
    #   included.cov = apply(dat.bt[,c(covariates, "outcome")],2,mean)
    
    ## calculate the weight
    
    dat.bt$lab = sign(dat.bt$weight)*dat.bt$A
    
    formula_str <- paste("as.factor(lab)", "~", paste(covariates, collapse = " + "))
    formula.owl <- as.formula(formula_str)
    
    mod.owl = wsvm(formula.owl
                   , data = dat.bt[dat.bt$train,], weight = abs(dat.bt$weight[dat.bt$train]), kernel = kernel.type, 
                   cross = 10, scale = T) #OWL method
    
     #mod.owl = best.tune_wsvm(formula.owl
    #                         , data = dat.bt[dat.bt$train,] , weight = abs(dat.bt$weight[dat.bt$train]), kernel = kernel.type, 
    #                         ranges = list(cost = c(0.01,0.1,1,5,10)), tunecontrol = tune.control(sampling = "cross") ) #OWL method
    
    
    ## MODEL FITTING
    
    fitted.owl =  predict(mod.owl, newdata = dat.bt[, covariates])# optTx(fitOWL, dat[dat$train == 0,])$optimalTx
    fitted.owl = as.numeric(as.character(fitted.owl))
    
    dat.bt$fitted.owl = fitted.owl
    
    owl.coef = svm.coef(mod.owl, type = "line")
    names(owl.coef) = c("Intercept", covariates)
    
    
    MAOB.same = apply(dat.bt[dat.bt$A == -1 & dat.bt$A == dat.bt$fitted.owl & dat.bt$train, covariates],2,mean, na.rm =T)
    MAOB.diff = apply(dat.bt[dat.bt$A == -1 & dat.bt$A != dat.bt$fitted.owl & dat.bt$train, covariates],2,mean, na.rm = T)
    
    DRA.same = apply(dat.bt[dat.bt$A == 1 & dat.bt$A == dat.bt$fitted.owl & dat.bt$train, covariates],2,mean, na.rm = T)
    DRA.diff = apply(dat.bt[dat.bt$A == 1 & dat.bt$A != dat.bt$fitted.owl & dat.bt$train, covariates],2,mean, na.rm = T)
    
    
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
    
    
    
    ########################################################################
    #################### Return results ####################################   
    
    value.ipw = c( value.opt, value.maob, value.dra, value.behavior)
    names(value.ipw) = c("Optimal Value Function", "MAOB Value Function", "DRA Value Function", "Behavior Policy")
    value.mr = c(value.opt.mr, value.maob.mr, value.dra.mr)
    names(value.mr) = c("Optimal Value Function MR", "MAOB Value Function MR", "DRA Value Function MR")
    var.value.mr = c( est.var.opt, est.var.maob, est.var.dra )
    names(var.value.mr) = c("Var Optimal Value Function MR", "Var MAOB Value Function MR", "Var DRA Value Function MR")
    
    prop.dra = mean(fitted.owl == 1)
    
    return(list(prop.dra = mean(prop.dra), 
                value.func.ipw = value.ipw,
                value.func.mr = value.mr, 
                var.value.func.mr = var.value.mr,  
                coef.beta = owl.coef, 
                MAOB.same = MAOB.same,
                MAOB.diff = MAOB.diff,
                DRA.same = DRA.same,
                DRA.diff = DRA.diff, 
                weight.prop = with(dat.bt[!dat.bt$train,], pmax(fz.pred*delta.pred,lb))))
  }, error = function(error_message){
    message(error_message)
    return(NA)
  })
  
} 

calculate.statistics = function(all_results, FUN){
  all_results = all_results[sapply(all_results, FUN = function(x) length(x) == length(all_results[[1]]))]
  average_return = lapply(2: (length(all_results[[1]]) - 1), function(j){
    component_values <- lapply(all_results, function(result) result[[j]])
    mat_comb = do.call(rbind, component_values)
    return(apply(mat_comb,2, FUN, na.rm = T))
    
  } )
  return(average_return)
}

proportion.mean.sd = function(all_results){
  
  mean.prop =  lapply(all_results, function(result) result[[1]]) %>% unlist() %>% mean() %>% round(.,2) 
  sd.prop = lapply(all_results, function(result) result[[1]]) %>% unlist() %>% sd() %>% round(.,2)
  mean.sd.prop = paste( mean.prop, " (", sd.prop, ")", sep = "" )
  return( mean.sd.prop )
}

read.updrs.list = function(all_results){
  all_results = all_results[sapply(all_results, FUN = function(x) length(x) == length(all_results[[1]]))] #remove missing data
  box.df = lapply(6:9, function(j){ #getting the distribution of MAOB same, MAOB diff, DRA same, and DRA diff
    component_values <- lapply(all_results, function(result) result[[j]])
    mat_comb = do.call(rbind, component_values)
    return(mat_comb)
    
  } )
  
  return(box.df)
}



wide.to.long = function(ind, type,box.df){
  
  temp.dat = as.data.frame(box.df[[ind]])
  temp.dat = temp.dat %>% pivot_longer(col = everything(), names_to = "Covariates", values_to = "Value")
  temp.dat$type = type
  
  return(temp.dat)
}

plot.box.func = function(all_results, title){
  
  box.df = read.updrs.list(all_results = all_results )  
  
  MAOB.same = wide.to.long(1, type = "MAOB same", box.df = box.df)
  MAOB.diff = wide.to.long(2, type = "MAOB diff", box.df = box.df)
  DRA.same = wide.to.long(3, type = "DRA same", box.df = box.df)
  DRA.diff = wide.to.long(4, type = "DRA diff", box.df = box.df) 
  
  plot.box.df = rbind(MAOB.same, MAOB.diff, DRA.same, DRA.diff)
  
  p1 = ggplot(data = plot.box.df, aes(x = Covariates, y = Value, col = as.factor(type))) + geom_boxplot() + theme_light() + labs(col = "Type", title = title) + theme(legend.position = "bottom")
  
  return(p1)
}



sim.updrs.func = function( FUN, dat, outcome.col, kernel, covariates, plot.title, lb, logistic){
  require(progress)
  
  total_iterations = 500
  pb <- progress_bar$new(
    format = "[:bar] :current/:total :percent eta: :eta",
    total = total_iterations
  )
  
  
  all_results_updrs = list()
  #all_results_updrs1.train = all_results_updrs2.train = all_results_updrs3.train = list()
  for(i in 1:total_iterations) {
   
    result_updrs <- FUN(seed = i, dat = dat, covariates = covariates, degree.poly = 3, kernel.type = kernel, outcome.col = outcome.col, center.outcome = T, scale.covariates = F, minimize = T, complier = F, lb = lb, logistic = logistic)
    
    all_results_updrs[[i]] <- result_updrs
    
    pb$tick()
  }
  
 # length(all_results_updrs[[1]])
  
  all_results_updrs = all_results_updrs[sapply(all_results_updrs, FUN = function(x) length(x) == 10 )]
  IPW.est <- lapply(all_results_updrs, function(result) result[[2]])
  IPW.est = do.call(rbind, IPW.est)
  
  DR.est <- lapply(all_results_updrs, function(result) result[[3]])
  DR.est = do.call(rbind, DR.est)
  
  opt.dec.coeff <- lapply(all_results_updrs, function(result) result[[5]])
  opt.dec.coeff = do.call(rbind, opt.dec.coeff)
  
  mean.updrs = calculate.statistics(all_results = all_results_updrs, FUN = mean)
  sd.updrs = calculate.statistics(all_results = all_results_updrs, FUN = sd)
  
  Mean.IPW.updrs = mean.updrs[[1]] %>% round(2)
  Mean.MR.updrs = mean.updrs[[2]] %>% round(2)
  sd.IPW.updrs = sd.updrs[[1]]  %>% round(2)
  sd.MR.updrs = sd.updrs[[2]]  %>% round(2)
  coeff.rule.mean = mean.updrs[[4]] %>% round(2)
  coeff.rule.sd = sd.updrs[[4]] %>% round(2)
  
  Mean.sd.IPW.updrs = paste0(Mean.IPW.updrs, " (", sd.IPW.updrs,")", sep = "")  
  Mean.sd.MR.updrs = paste0(Mean.MR.updrs, " (", sd.MR.updrs,")", sep = "")
  Mean.sd.coeff.rule = paste0(coeff.rule.mean, " (", coeff.rule.sd,")", sep = "")
  
  decision.plt = plot.box.func(all_results_updrs, title = plot.title)
  DRA.prop = proportion.mean.sd(all_results = all_results_updrs)
  
  return(list(IPW = Mean.sd.IPW.updrs, MR = Mean.sd.MR.updrs, DRA.prop = DRA.prop, plot = decision.plt, 
              decision.coef = Mean.sd.coeff.rule, DR.est = DR.est, IPW.est = IPW.est, opt.dec.coeff = opt.dec.coeff))
}


#remove the outliers and refitting the nuisance parameters based on the new data without the outlier
remove.outlier.func = function(dat, covariates, outcome.col, minimize = T, complier = F, center.outcome = T, scale.covariates = T, lb = 0, alpha.lo, alpha.hi){
  tryCatch({ 
    
    
    #### Identify and remove outlier
    save.dat = dat     
    if(any(covariates %in% c("dys", "fluc", "nausea", "sleep", "ortho")) ){
      dat = dat %>% filter(!is.na(dys))
    }
    
    #prepare data
    if(minimize){
      dat$outcome = dat[,outcome.col]*-1
    }else{dat$outcome = dat[,outcome.col]}
    
    
    if(center.outcome){
      dat$outcome = dat$outcome - mean(dat$outcome)
    }
    ############calculate nuisance parameters
    #propensity model
    
    if(scale.covariates){
      dat[,c("updrs_1", "updrs_2", "updrs_3", "time_from_diagnosis", "time_from_diagnosis.2nd", "levodopa_dose","levodopa_dose.2nd",
             "age", "age.2nd", "bmi")] = apply( dat[,c("updrs_1", "updrs_2", "updrs_3", "time_from_diagnosis", "time_from_diagnosis.2nd", "levodopa_dose","levodopa_dose.2nd", "age", "age.2nd", "bmi")],2,scale)
    }
    
    
    
    dat.bt = dat
    #dat.bt$train =  #dat.bt$study == "ls1"
    
    logistic.A = glm(as.factor(A) ~ Z +  updrs_1 + updrs_2 + updrs_3 + time_from_diagnosis.2nd + levodopa_dose.2nd + gender + age.2nd + bmi , data = dat.bt, family = binomial(link = "logit"))
    
    pred.A1 = predict(logistic.A, type = "response", newdata = dat.bt)
    dat.bt$pred.A = ifelse(dat.bt$A == 1, pred.A1, 1 - pred.A1)
    
    dat.Z1 = dat.bt %>% mutate(Z = 1)
    dat.bt$pred.A1.Z1 = predict(logistic.A, newdata = dat.Z1, type = "response")
    
    dat.Zn1 = dat.bt %>% mutate(Z = -1)
    dat.bt$pred.A1.Zn1 = predict(logistic.A, newdata = dat.Zn1, type = "response")
    
    
    logistic.Z = glm(as.factor(Z) ~ updrs_1 + updrs_2 + updrs_3 + time_from_diagnosis + levodopa_dose + gender + age + bmi , data = dat.bt, family = binomial(link = "logit"))
    
    pred.z = predict(logistic.Z, newdata = dat.bt, type = "response") 
    dat.bt$fz.pred = ifelse(dat.bt$Z == 1, pred.z , 1 - pred.z) 
    
    if(complier){
      dat.bt$delta.pred = rep(sum(dat.bt$A == 1 & dat.bt$Z == 1)/sum(dat.bt$Z == 1) - sum(dat.bt$A == 1 & dat.bt$Z == -1)/sum(dat.bt$Z == -1), nrow(dat.bt))
      weight = dat.bt$outcome[dat.bt$train]*dat.bt$A[dat.bt$train]*dat.bt$Z[dat.bt$train]/(dat.bt$fz.pred[dat.bt$train]) #OWL weight
    }else{
      dat.bt$delta.pred = dat.bt$pred.A1.Z1 - dat.bt$pred.A1.Zn1
      dat.bt$weight = with(dat.bt, outcome*A*Z/( pmax(fz.pred*delta.pred,lb) ) )
    }
    
    outlier.cov = dat.bt %>% filter(weight <= quantile(weight, prob = alpha.lo) | weight > quantile(weight, prob = alpha.hi) )
    outlier.id = outlier.cov$patno %>% unique()
    outlier.weight = outlier.cov$weight
    outlier.cov = outlier.cov[, c("patno", covariates, "outcome")]
    
    
    dat.bt = dat.bt %>% filter(weight > quantile(weight, prob = alpha.lo) & weight <= quantile(weight, prob = alpha.hi) )
    non.outlier.weight = dat.bt$weight
    
    dat = save.dat
    dat = dat %>% filter(!(patno %in% outlier.id) )
    
    return( list(dat = dat, outlier.id = outlier.id, outlier.weight = outlier.weight, non.outlier.weight = non.outlier.weight, outlier.dat = outlier.cov, non.outlier.cov = dat.bt[,c("patno", covariates, "outcome")]))
    
  }, error = function(error_message){
    message(error_message)
    return(NA)
  })
  
}


## plotting the decision covariates 

plt.decision = function(seed, dat, covariates, target.covariate, kernel.type, degree.poly,  outcome.col, minimize = T, complier = F, center.outcome = T, lb = 0){
  tryCatch({ 
    set.seed(seed)
    require(WeightSVM)
    
    if (!all(target.covariate %in% covariates) ) {
      stop("The target covariate must be one of the covariates")
    }
    
    if(any(covariates %in% c("dys", "fluc", "nausea", "sleep", "ortho")) ){
      dat = dat %>% filter(!is.na(dys))
    }
    
    #prepare data
    if(minimize){
      dat$outcome = dat[,outcome.col]*-1
    }else{dat$outcome = dat[,outcome.col]}
    
    
    if(center.outcome){
      dat$outcome = dat$outcome - mean(dat$outcome)
    }
    ############calculate nuisance parameters
    #propensity model
    
    dat.bt = dat
    dat.bt$train = sample(c(T,F), nrow(dat.bt), prob = c(0.7,0.3), replace = T )
    
    logistic.A = glm(as.factor(A) ~ Z +  updrs_1 + updrs_2 + updrs_3 + 
                       time_from_diagnosis.2nd + levodopa_dose.2nd + gender + age.2nd + bmi + Time_since_first_Levo , data = dat.bt, family = binomial(link = "logit"))
    
    pred.A1 = predict(logistic.A, type = "response", newdata = dat.bt)
    dat.bt$pred.A = ifelse(dat.bt$A == 1, pred.A1, 1 - pred.A1)
    
    dat.Z1 = dat.bt %>% mutate(Z = 1)
    dat.bt$pred.A1.Z1 = predict(logistic.A, newdata = dat.Z1, type = "response")
    
    dat.Zn1 = dat.bt %>% mutate(Z = -1)
    dat.bt$pred.A1.Zn1 = predict(logistic.A, newdata = dat.Zn1, type = "response")
    
    
    logistic.Z = glm(as.factor(Z) ~ updrs_1 + updrs_2 + updrs_3 + time_from_diagnosis
                     + levodopa_dose + gender + age + bmi + Time_since_first_Levo , data = dat.bt, family = binomial(link = "logit"))
    
    pred.z = predict(logistic.Z, newdata = dat.bt, type = "response") 
    dat.bt$fz.pred = ifelse(dat.bt$Z == 1, pred.z , 1 - pred.z) 
    
    if(complier){
      dat.bt$delta.pred = rep(sum(dat.bt$A == 1 & dat.bt$Z == 1)/sum(dat.bt$Z == 1) - sum(dat.bt$A == 1 & dat.bt$Z == -1)/sum(dat.bt$Z == -1), nrow(dat.bt))
      weight = dat.bt$outcome[dat.bt$train]*dat.bt$A[dat.bt$train]*dat.bt$Z[dat.bt$train]/(dat.bt$fz.pred[dat.bt$train]) #OWL weight
    }else{
      dat.bt$delta.pred = dat.bt$pred.A1.Z1 - dat.bt$pred.A1.Zn1
      dat.bt$weight = with(dat.bt, outcome*A*Z/( pmax(fz.pred*delta.pred,lb) ) )
    }
    
    #uncomment to remove the outliers by weight.  
    # outlier.cov = dat.bt %>% filter(weight <= quantile(weight, prob = 0.05) | weight >= quantile(weight, prob = 0.95) )
    # outlier.id = outlier.cov$patno %>% unique()
    # outlier.weight = outlier.cov$weight
    # outlier.cov = apply(outlier.cov[, c(covariates, "outcome")],2,mean)
    # 
    
    #dat.bt = dat.bt %>% filter(weight > quantile(weight, prob = 0.05) & weight < quantile(weight, prob = 0.95) )
    # non.outlier.weight = dat.bt$weight
    dat.bt$fitted.owl = NA
    
    #   included.cov = apply(dat.bt[,c(covariates, "outcome")],2,mean)
    
    ## calculate the weight
    
    dat.bt$lab = sign(dat.bt$weight)*dat.bt$A
    
    formula_str <- paste("as.factor(lab)", "~", paste(covariates, collapse = " + "))
    formula.owl <- as.formula(formula_str)
    
    mod.owl = wsvm(formula.owl
                   , data = dat.bt[dat.bt$train,], degree = degree.poly , weight = abs(dat.bt$weight[dat.bt$train]), kernel = kernel.type, 
                   cross = 10, scale = T) #OWL method
    
    #mod.owl = best.tune_wsvm(formula.owl
    #                         , data = dat.bt[dat.bt$train,] , weight = abs(dat.bt$weight[dat.bt$train]), kernel = kernel.type, 
    #                         ranges = list(cost = c(0.01,0.1,1,5,10)), tunecontrol = tune.control(sampling = "cross"), scale = T ) #OWL method
    
    
    ## MODEL FITTING
    
    grid_list <- lapply(covariates, function(x) {
      # Check if the covariate is numeric
    
      if (x %in% target.covariate ) {
        # Generate a sequence for numeric covariates
        return(seq(min(dat[[x]], na.rm = TRUE), 
                   max(dat[[x]], na.rm = TRUE), length.out = 50))
      } else {
        # For non-numeric covariates (or if you prefer a fixed value), return the mean
        if(length(unique(dat[[x]])) == 2){
          return( sort(unique(dat[[x]])) )
        }else{
          return(mean(dat[[x]], na.rm = TRUE))
        }
        
      }
    })
    
    grid.data = do.call(expand.grid, grid_list)
    colnames(grid.data) = covariates
    
    fitted.owl =  predict(mod.owl, newdata = grid.data)# optTx(fitOWL, dat[dat$train == 0,])$optimalTx
    fitted.owl = as.numeric(as.character(fitted.owl))
    grid.data$trt = (fitted.owl == 1)*1
    
    return(grid.data)
    
  }, error = function(error_message){
    message(error_message)
    return(NA)
  })
  
} 


plt.decision.owl = function(seed, dat, covariates, target.covariate, kernel.type, degree.poly,  outcome.col, minimize = T, complier = F, center.outcome = T, lb = 0){
  tryCatch({ 
    set.seed(seed)
    require(WeightSVM)
    
    if (!all(target.covariate %in% covariates) ) {
      stop("The target covariate must be one of the covariates")
    }
    
    if(any(covariates %in% c("dys", "fluc", "nausea", "sleep", "ortho")) ){
      dat = dat %>% filter(!is.na(dys))
    }
    
    #prepare data
    if(minimize){
      dat$outcome = dat[,outcome.col]*-1
    }else{dat$outcome = dat[,outcome.col]}
    
    
    if(center.outcome){
      dat$outcome = dat$outcome - mean(dat$outcome)
    }
    dat.bt = dat
    dat.bt$train = sample(c(T,F), nrow(dat.bt), prob = c(0.7,0.3), replace = T )
    
    logistic.A = glm(as.factor(A) ~ Z +  updrs_1 + updrs_2 + updrs_3 + 
                       time_from_diagnosis.2nd + levodopa_dose.2nd + gender
                       + age.2nd + bmi + Time_since_first_Levo , data = dat.bt, family = binomial(link = "logit"))
    
    pred.A1 = predict(logistic.A, type = "response", newdata = dat.bt)
    dat.bt$pred.A = ifelse(dat.bt$A == 1, pred.A1, 1 - pred.A1)
    
    logistic.A.v2 = glm(as.factor(A) ~  updrs_1 + updrs_2 + updrs_3 + 
                          time_from_diagnosis.2nd + levodopa_dose.2nd + gender + age.2nd + bmi + Time_since_first_Levo , data = dat.bt, family = binomial(link = "logit"))
    
    pred.A1.v2 = predict(logistic.A.v2, type = "response", newdata = dat.bt)
    dat.bt$pred.A.v2 = ifelse(dat.bt$A == 1, pred.A1.v2, 1 - pred.A1.v2)
    
    dat.Z1 = dat.bt %>% mutate(Z = 1)
    dat.bt$pred.A1.Z1 = predict(logistic.A, newdata = dat.Z1, type = "response")
    
    dat.Zn1 = dat.bt %>% mutate(Z = -1)
    dat.bt$pred.A1.Zn1 = predict(logistic.A, newdata = dat.Zn1, type = "response")
    
    
    logistic.Z = glm(as.factor(Z) ~ updrs_1 + updrs_2 + updrs_3 + time_from_diagnosis + 
                       levodopa_dose + gender + age + bmi + Time_since_first_Levo , data = dat.bt, family = binomial(link = "logit"))
    
    pred.z = predict(logistic.Z, newdata = dat.bt, type = "response") 
    dat.bt$fz.pred = ifelse(dat.bt$Z == 1, pred.z , 1 - pred.z) 
    
    if(complier){
      dat.bt$delta.pred = rep(sum(dat.bt$A == 1 & dat.bt$Z == 1)/sum(dat.bt$Z == 1) - sum(dat.bt$A == 1 & dat.bt$Z == -1)/sum(dat.bt$Z == -1), nrow(dat.bt))
      weight = dat.bt$outcome[dat.bt$train]*dat.bt$A[dat.bt$train]*dat.bt$Z[dat.bt$train]/(dat.bt$fz.pred[dat.bt$train]) #OWL weight
    }else{
      dat.bt$delta.pred = dat.bt$pred.A1.Z1 - dat.bt$pred.A1.Zn1
      dat.bt$weight = with(dat.bt, outcome/( pmax(pred.A.v2,lb) ) )
    }
    
    #uncomment to remove the outliers by weight.  
    # outlier.cov = dat.bt %>% filter(weight <= quantile(weight, prob = 0.05) | weight >= quantile(weight, prob = 0.95) )
    # outlier.id = outlier.cov$patno %>% unique()
    # outlier.weight = outlier.cov$weight
    # outlier.cov = apply(outlier.cov[, c(covariates, "outcome")],2,mean)
    # 
    
    #dat.bt = dat.bt %>% filter(weight > quantile(weight, prob = 0.05) & weight < quantile(weight, prob = 0.95) )
    # non.outlier.weight = dat.bt$weight
    dat.bt$fitted.owl = NA
    
    #   included.cov = apply(dat.bt[,c(covariates, "outcome")],2,mean)
    
    ## calculate the weight
    
    dat.bt$lab = sign(dat.bt$weight)*dat.bt$A
    
    formula_str <- paste("as.factor(lab)", "~", paste(covariates, collapse = " + "))
    formula.owl <- as.formula(formula_str)
    
    mod.owl = wsvm(formula.owl
                   , data = dat.bt[dat.bt$train,], degree = degree.poly , weight = abs(dat.bt$weight[dat.bt$train]), kernel = kernel.type, 
                   cross = 10, scale = T) #OWL method
    
    ## MODEL FITTING
    
    grid_list <- lapply(covariates, function(x) {
      # Check if the covariate is numeric
      
      if (x %in% target.covariate ) {
        # Generate a sequence for numeric covariates
        return(seq(min(dat[[x]], na.rm = TRUE), 
                   max(dat[[x]], na.rm = TRUE), length.out = 50))
      } else {
        # For non-numeric covariates (or if you prefer a fixed value), return the mean
        if(length(unique(dat[[x]])) == 2){
          return( sort(unique(dat[[x]])) )
        }else{
          return(mean(dat[[x]], na.rm = TRUE))
        }
        
      }
    })
    
    grid.data = do.call(expand.grid, grid_list)
    colnames(grid.data) = covariates
    
    fitted.owl =  predict(mod.owl, newdata = grid.data)# optTx(fitOWL, dat[dat$train == 0,])$optimalTx
    fitted.owl = as.numeric(as.character(fitted.owl))
    grid.data$trt = (fitted.owl == 1)*1
    
    return(grid.data)
    
  }, error = function(error_message){
    message(error_message)
    return(NA)
  })
  
} 



calc_phi <- function(table) {
  # Check if the table is 2x2
  if (!all(dim(table) == c(2, 2))) {
    stop("The input must be a 2x2 contingency table.")
  }
  
  # Extract the values from the table
  a <- table[1, 1]
  b <- table[1, 2]
  c <- table[2, 1]
  d <- table[2, 2]
  
  # Calculate the Phi coefficient
  phi_value <- (a * d - b * c) / sqrt((a + b) * (a + c) * (b + d) * (c + d))
  
  return(phi_value)
}

flip.decision.rule = function(s){
  
  # Split the string into two parts
  parts <- strsplit(s, " ")[[1]]
  
  # Convert the first part to a negative number
  negative_value <- as.character(-as.numeric(parts[1]))
  
  # Combine them back
  new_string <- paste(negative_value, parts[2])
  
  # Output the result
  new_string
}







