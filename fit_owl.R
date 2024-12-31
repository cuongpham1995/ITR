
res2 = fit.iv.owl.bt.ls1.postcept(seed = 3431, dat = iv.dat2a, covariates = sel.cov, 
                           kernel.type = "linear", degree.poly = 2, outcome.col = outcome,
                           minimize = T, complier = F, center.outcome = T, scale.covariates = F, lb = 0, logistic = T)

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
    
    #return(table(fitted.owl))
    
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
    
    return(c(value.opt.mr, value.maob.mr, value.dra.mr))
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

