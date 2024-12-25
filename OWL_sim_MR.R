#in this codes, we will try to see if the MR for OWL actually works 
#this is the scenario where the DRA is more dominant than MAOB

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

dat.gen.owl =function(seed, nsample){
  set.seed(seed)
  require(mvtnorm)
  require(dplyr)
  
  X <- matrix(runif(n = 3 * nsample, min = -1, max = 1), ncol = 3)
  sigma.mat = diag(nsample)*1
  A = as.matrix(2*rbinom(nsample, size = 1, prob = 0.5) - 1)
  
  #Y = rmvnorm(n =1, mean = 1 + 2*X[,1] + X[,2] + 0.5*X[,3] + 0.442*(1 - X[,1] - X[,2])*A[,1], sigma = sigma.mat) %>% t()
  Y = rmvnorm(n =1, mean = 1 + (0.5 - 0.5*X[,1] - X[,2])*A[,1], sigma = sigma.mat) %>% t()
  
  #potential values
  Y1 =as.matrix( 1 +  (0.5 - 0.5*X[,1] - X[,2])*1)
  Yn1 = as.matrix( 1+  (0.5 - 0.5*X[,1] - X[,2])*(-1)) 
  
  #table(Y1 > Yn1)
  
  dat = cbind(Y, A, X, Y1, Yn1) %>% as.data.frame() 
  colnames(dat) = c("Y", "A", paste("X", seq(1:ncol(X)), sep = "."), "Y1", "Yn1")
  dat$opt.trt = (dat$Y1 > dat$Yn1)*2 - 1 
  return(dat)
}


fit.owl.mr.sim = function(seed, nsample, kernel.type = "linear", method = "ipw" ){
  tryCatch({ 
    set.seed(seed)
    require(WeightSVM) #to fit weighted svm
    require(glmnet) #for lasso model
    require(dplyr) 
    require(hal9001)
    #require(survival)
    dat = dat.gen.owl(seed = seed, nsample = nsample)
    
    #prepare data
    dat$train = rbinom(nrow(dat), size = 1, prob = 0.8)
    
    train_data <- dat[dat$train == 1, c("Y", "A", "X.1", "X.2", "X.3", "Y1", "Yn1")]
    test_data  <- dat[dat$train == 0, c("Y", "A", "X.1", "X.2", "X.3", "Y1", "Yn1")]
    
    ############calculate nuisance parameters
    #propensity model
    
    prop.mod.cv =  fit_hal(X = train_data[,-c(1:2)], Y = train_data$A, family = "binomial", smoothness_orders = 0,
                           fit_control = list(cv_select = T, n_folds = 10, use_min = T), reduce_basis = NULL)
    
    
    pred.a.train = predict(prop.mod.cv, new_data = train_data[,-c(1:2)]) 
    fz.train = ifelse(dat$A[dat$train == 1] == 1, pred.a.train, 1 - pred.a.train)
    
    pred.a.test = predict(prop.mod.cv, new_data = test_data[,-c(1:2)]) #predict(cv.prop, s = cv.prop$lambda.min, newx = cov.x.test, type ="response")
    fz.test = ifelse(dat$A[dat$train == 0] == 1, pred.a.test, 1 - pred.a.test)
    
    #mean outcome model
    
    # Prepare your data
    train_data$A = as.numeric(train_data$A)
    test_data$A = as.numeric(test_data$A)
    
    #fit HAL model for the mean outcome model
    mean.mod.cv = fit_hal(X = train_data[,-1], Y = train_data$Y, family = "gaussian", smoothness_orders = 0,
                          fit_control = list(cv_select = T, n_folds = 10, use_min = T), reduce_basis = NULL)
    
    # Predict on the train data
    train.outcome <- predict(mean.mod.cv, new_data = train_data[,-1])
    train_data$A <- 1
    train.outcome.A1 <- predict(mean.mod.cv, new_data = train_data[, -1])
    train_data$A <- -1
    train.outcome.An1 <- predict(mean.mod.cv, new_data = train_data[, -1])
    
    # Predict on the test data
    test.outcome <- predict(mean.mod.cv, new_data = test_data[,-1])
    
    test_data$A <- 1
    test.outcome.A1 <- predict(mean.mod.cv, new_data = test_data[, -1])
    test_data$A <- -1
    test.outcome.An1 <- predict(mean.mod.cv, new_data = test_data[, -1])
    
    if(method == "mr"){
      weight = (dat$Y[dat$train == 1] - train.outcome)/fz.train + dat$A[dat$train == 1]*(train.outcome.A1 - train.outcome.An1) 
    }else if(method == "ipw"){
      weight = dat$Y[dat$train ==1]/fz.train #OWL weight
    }
    
    mean(dat$A[dat$train == 1]*(train.outcome.A1 - train.outcome.An1) )
    
    lab = sign(weight)*dat$A[dat$train == 1]
    
    
    #models building
    
    mod.owl = wsvm(as.factor(lab) ~ X.1 + X.2 + X.3
                   , data = dat[dat$train == 1,], weight = abs(weight), kernel = kernel.type, 
                   cross = 10, scale = F) #OWL method
    
    
    ## MODEL FITTING
    
    fitted.owl =  predict(mod.owl, newdata = dat[dat$train == 0, c("X.1", "X.2", "X.3")])# optTx(fitOWL, dat[dat$train == 0,])$optimalTx
    fitted.owl = as.numeric(as.character(fitted.owl))
    
    fitted.maob = rep(-1, length(fitted.owl))
    fitted.dra = rep(1, length(fitted.owl))
    
    fitted.owl.train = predict(mod.owl, newdata = dat[dat$train == 1, c("X.1", "X.2", "X.3")])
    fitted.owl.train = as.numeric(as.character(fitted.owl.train))
    
    fitted.maob.train = rep(-1, length(fitted.owl.train))
    fitted.dra.train = rep(1, length(fitted.owl.train))
    
    #correct rate
    #dat$opt.trt = ifelse(dat$Y1 > dat$Yn1,1,-1) 
    #tab.correct.rate = table( dat$opt.trt[dat$train == 0], fitted.owl)
    correct.rate = mean(fitted.owl == dat$opt.trt[dat$train == 0])
    
    #true value function
    true.value.opt = mean(ifelse(fitted.owl == 1, dat$Y1[dat$train == 0], dat$Yn1[dat$train == 0]))
    true.value.maob = mean(dat$Yn1[dat$train == 0])
    true.value.dra = mean(dat$Y1[dat$train == 0])
    
    #print(c( true.value.opt,true.value.maob,true.value.dra))
    
    #estimated value function
    value.opt = mean((dat$Y[dat$train == 0]*I(fitted.owl == dat$A[dat$train == 0]))/fz.test)
    value.maob =  mean((dat$Y[dat$train == 0]*I(fitted.maob == dat$A[dat$train == 0]) )/fz.test )
    value.dra = mean((dat$Y[dat$train == 0]*I(fitted.dra == dat$A[dat$train == 0]))/fz.test)
    value.behavior =  mean(dat$Y[dat$train == 0]) 
    
    #print(c(value.opt, value.maob, value.dra, value.behavior))
    
    value.opt.train = mean((dat$Y[dat$train == 1]*I(fitted.owl.train == dat$A[dat$train == 1]))/fz.train)
    value.maob.train =  mean((dat$Y[dat$train == 1]*I(fitted.maob.train == dat$A[dat$train == 1]) )/fz.train )
    value.dra.train = mean((dat$Y[dat$train == 1]*I(fitted.dra.train == dat$A[dat$train == 1]))/fz.train)
    value.behavior.train =  mean(dat$Y[dat$train == 1]) 
    
    #print(c(value.opt.train, value.maob.train, value.dra.train))
    
    #value function mr
    value.opt.mr = value.opt - mean( (test.outcome*I(fitted.owl == dat$A[dat$train == 0]))/fz.test 
                                     - I( fitted.owl == 1)*test.outcome.A1 - I( fitted.owl == -1)*test.outcome.An1) 
    
    value.maob.mr = value.maob - mean( (test.outcome*I(fitted.maob == dat$A[dat$train == 0]))/fz.test 
                                       - I(fitted.maob == 1)*test.outcome.A1 - I(fitted.maob == -1)*test.outcome.An1) 
    
    value.dra.mr = value.dra - mean( (test.outcome*I(fitted.dra == dat$A[dat$train == 0]))/fz.test 
                                     - I(fitted.dra == 1)*test.outcome.A1 - I(fitted.dra == -1)*test.outcome.An1) 
    
    #variance of the mr
    var.value.opt.mr = var((dat$Y[dat$train == 0]*I(fitted.owl == dat$A[dat$train == 0]))/fz.test
                           - (test.outcome*I(fitted.owl == dat$A[dat$train == 0]))/fz.test 
                           + I( fitted.owl == 1)*test.outcome.A1 + I( fitted.owl == -1)*test.outcome.An1
    )*(1/length(dat$Y[dat$train == 0]))
    
    var.value.maob.mr = var((dat$Y[dat$train == 0]*I(fitted.maob == dat$A[dat$train == 0]))/fz.test
                            - (test.outcome*I(fitted.maob == dat$A[dat$train == 0]))/fz.test 
                            + I( fitted.maob == 1)*test.outcome.A1 + I( fitted.maob == -1)*test.outcome.An1
    )*(1/length(dat$Y[dat$train == 0]))
    
    var.value.dra.mr = var((dat$Y[dat$train == 0]*I(fitted.dra == dat$A[dat$train == 0]))/fz.test
                           - (test.outcome*I(fitted.dra == dat$A[dat$train == 0]))/fz.test 
                           + I( fitted.dra == 1)*test.outcome.A1 + I( fitted.dra == -1)*test.outcome.An1
    )*(1/length(dat$Y[dat$train == 0]))
    
    #owl.coef = regimeCoef(fitOWL)
    owl.coef = svm.coef(mod.owl)
    value.mr = c(value.opt.mr, value.maob.mr, value.dra.mr) 
    value.ipw = c(value.opt, value.maob, value.dra, value.behavior)
    value.ipw.train = c(value.opt.train, value.maob.train, value.dra.train, value.behavior.train)
    prop.dra = c(mean(fitted.owl == 1), mean(fitted.owl.train == 1) )
    #names(res)[1:4] =  c("prop.fit.high.dose", "value.opt","value.maob","value.dra")
    
    
    #names(res)[5:length(res)] = names(owl.coef)
    #return(list(prop.fitted = mean(fitted.owl == 1), svm.coef = owl.coef, value.function = value.owl, value.obs = value.obs, 
    #           value.dopa = value.dopa, value.levo = value.levo))
    return(list(prop.dra = prop.dra, value.func.ipw = value.ipw, value.func.mr = value.mr, value.func.ipw.train = value.ipw.train, 
                coef.beta = owl.coef, var.est.mr = c(var.value.opt.mr, var.value.maob.mr, var.value.dra.mr), 
                #true value
                true.vf = c(true.value.opt, true.value.maob, true.value.dra), correct_rate = correct.rate))
  }, error = function(error_message){
    message(error_message)
    return(NA)
  })
  
}

#temp = fit.owl.mr.sim(seed = 92231, nsample = 500, kernel.type = "linear", method = "mr" )

##### Simulation ####
SLURM_ID <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
n = 500
method.use = "mr"
#calculate the true value
res = fit.owl.mr.sim(seed = SLURM_ID, nsample = n, kernel.type = "linear", method = method.use )
max_length <- max(sapply(res, length))
res = lapply(res, function(x) { if (length(x) < max_length) c(x, rep(NA, max_length - length(x))) else x })

output_correct <- paste0("HAL_IPW_",n,"_",method.use)
if (!file.exists(output_correct)){dir.create(output_correct)}

############## Save files ###############################
write.csv(res, file = paste(output_correct,"/result", SLURM_ID,".csv", sep = ""))


