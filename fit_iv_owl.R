fit.iv.owl = function(itr_formula, iv_formula, trt_formula, mean_formula = NULL,method, trt, iv, outcome, dat, 
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
  trt_formula <- all.vars(trt_formula)
  
  missing_itr_vars <- setdiff(itr_vars, names(dat))
  missing_iv_vars <- setdiff(iv_vars, names(dat))
  missing_trt_vars <- setdiff(trt_vars, names(dat))
  
  if(length(missing_itr_vars) > 0){
    stop(paste("The following variables in itr_formula are missing from the dataset:", 
               paste(missing_itr_vars, collapse = ", ")))
  }
  
  if(length(missing_prop_vars) > 0){
    stop(paste("The following variables in prop_formula are missing from the dataset:", 
               paste(missing_prop_vars, collapse = ", ")))
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
    columns_to_check <- unique(c(itr_vars, prop_vars, "A", "Y"))
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
  iv.logistic.formula = paste( "as.factor(Z)" , deparse(iv_formula), collapse = "")  
  
  iv.logistic = glm(iv.logistic.formula, data = dat, family = binomial(link = "logit"))
  
  pred.Z = predict(iv.logistic, type = "response", newdata = dat)
  dat$Z.pred = ifelse(dat$Z == 1, pred.Z, 1 - pred.Z)
  
  ## Trt model
  #convert to a formula
  trt.logistic.formula = paste( "as.factor(A)" , deparse(trt_formula), collapse = "") 
  
  trt.logistic = glm(trt.logistic.formula, data = dat, family = binomial(link = "logit"))
  
  pred.A = predict(trt.logistic, type = "response", newdata = dat)
  dat$A.pred = ifelse(dat$A == 1, pred.A, 1 - pred.A)
  
  ## Mean outcome model
  
  #convert to a formula
  if(is.null(mean_formula)){
    lm.formula = paste("Y", deparse(prop_formula), collapse = "" )
  }else{ lm.formula = paste("Y", deparse(mean_formula), collapse = "" )  }
  
  lm.mod = lm(lm.formula, data = dat)
  
  # Add predictions to the dataset
  dat$Y.pred <- predict(lm.mod)
  
  #estimating E[Y|A = 1, ...]
  dat_mod = dat
  dat_mod$A = 1
  dat$Y.A1 = predict(lm.mod, newdata = dat_mod)
  #estimating E[Y|A = -1, ...]
  dat_mod$A = -1
  dat$Y.An1 = predict(lm.mod, newdata = dat_mod)
  
  
  }