#' Fit an Outcome-Weighted Learning (OWL) Model
#'
#' This function fits an OWL model using a weighted SVM (wSVM) to optimize individual treatment rules.
#'
#' @param itr_formula A formula specifying the independent variables for treatment rule optimization.
#' @param prop_formula A formula specifying the independent variables for propensity score modeling.
#' @param mean_formula A formula specifying the independent variables for mean outcome model.
#' @param method Choose either "mr" for multiply robust method or "ipw" for inverse probability weighting method
#' @param trt A character string specifying the name of the treatment variable in the dataset.
#' @param dat A data frame containing the data for the model.
#' @param kernel A character string specifying the kernel type for the SVM (default is "linear").
#' @param scale Logical, whether to scale the data for the SVM (default is TRUE).
#' @param cross An integer specifying the number of folds for cross-validation in the SVM (default is 10).
#' @return A numeric vector of predicted treatment decisions.
#' @examples
#' # Example usage
#' sim_dat = dat.gen.owl(23, nsample = 250)
#'
#' model.owl = fit.owl(itr_formula = ~ X.1 + X.2 + X.3, prop_formula = ~ X.1 + X.2, mean_formula = ~ A*(X.1 + X.2) - X.1 - X.2, method = "mr", trt = "A", outcome = "outcome", dat = sim_dat,
#'                    center.outcome = F, seed = 98734)
#'
#' pred.val = predict(model.owl$model, newdata = sim_dat)
#' 
#' dat <- data.frame(Y = rnorm(100), A = rbinom(100, 1, 0.5), X1 = rnorm(100), X2 = rnorm(100))
#' itr_formula <- ~ X1 + X2
#' prop_formula <- ~ X1 + X2
#' fit.owl(itr_formula, prop_formula, trt = "A", dat = dat)

fit.owl = function(itr_formula, prop_formula, mean_formula = NULL,method, trt, outcome, dat, 
                   na.rm = T, minimize = F, center.outcome = F, 
                   kernel = "linear", cross = 10, scale = T, seed, ...){
    set.seed(seed)
    ## Error handling
    if (!inherits(itr_formula, "formula")) stop("itr_formula must be a formula object.")
    if (!inherits(prop_formula, "formula")) stop("prop_formula must be a formula object.")
    if (!inherits(mean_formula, "formula") & !is.null(mean_formula)) stop("mean_formula must be a formula object.")
    if (!trt %in% names(dat)) stop("The treatment variable is not in the dataset.")
    if (!outcome %in% names(dat)) stop("The outcome variable is not in the dataset.")
    
    
    ## Extract variables from formulas
    itr_vars <- all.vars(itr_formula)
    prop_vars <- all.vars(prop_formula)
   
    
    missing_itr_vars <- setdiff(itr_vars, names(dat))
    missing_prop_vars <- setdiff(prop_vars, names(dat))
   
    if(length(missing_itr_vars) > 0){
       stop(paste("The following variables in itr_formula are missing from the dataset:", 
               paste(missing_itr_vars, collapse = ", ")))
      }
  
    if(length(missing_prop_vars) > 0){
       stop(paste("The following variables in prop_formula are missing from the dataset:", 
               paste(missing_prop_vars, collapse = ", ")))
    }
    
    # dealing with the case that the mean formula is not null
    if(!is.null(mean_formula)){
      mean_vars <- all.vars(mean_formula)
      
      missing_mean_vars <- setdiff(mean_vars, names(dat))
      if(length(missing_mean_vars) > 0){
        stop(paste("The following variables in mean_formula are missing from the dataset:", 
                   paste(missing_mean_vars, collapse = ", ")))
      }
    }
  
    ## Data preparing
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
  
    ## Propensity score model
    
    #convert to a formula
    logistic.formula = paste( "as.factor(A)" , deparse(prop_formula), collapse = "")  
    
    logistic.A = glm(logistic.formula, data = dat, family = binomial(link = "logit"))
    
    pred.A = predict(logistic.A, type = "response", newdata = dat)
    dat$pred = ifelse(dat$A == 1, pred.A, 1 - pred.A)
    
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
    
    ## calculate the weight
    
    dat$weight = if(method == "ipw"){
                    with(dat, Y/pred)
    }else if(method == "mr"){
      
      with(dat, (dat$Y - Y.pred)/pred + A*(Y.A1 - Y.An1)) 
    }else{
      stop("Invalid method. Choose either 'ipw' or 'mr'.")
    }
    

    dat$fitted.owl = NA
    
    dat$lab = sign(dat$weight)*dat$A
    
    
    formula_str <- paste("as.factor(lab)", deparse(itr_formula), collapse = "")
    formula.owl <- as.formula(formula_str)
    
    ## Model fitting
    mod.owl = WeightSVM::wsvm(formula.owl, data = dat[dat$train,], weight = abs(dat$weight[dat$train]), kernel = "linear", 
                   cross = 10, scale = T) #OWL method
    
    # Predicting 
    fitted.owl = predict(mod.owl, newdata = dat[!dat$train,])
    
    
    
    fitted.owl = as.numeric(as.character(fitted.owl)) #opt treatment regimes
    fitted.A = rep(-1, length(fitted.owl)) #always -1 strategy. A is -1 
    fitted.B = rep(1, length(fitted.owl)) #always 1 strategy. B is 1 
    
    
    #return(table(fitted.owl, true.value[!dat$train]))
    
    owl.coef = svm.coef(mod.owl, type = "line")
    
    
    ## Evaluating the value function

    # estimating value function with IPW
    value.opt = mean((dat$Y[!dat$train]*I(fitted.owl == dat$A[!dat$train]))/dat$pred[!dat$train])
    value.A =  mean((dat$Y[!dat$train]*I(fitted.A == dat$A[!dat$train]) )/dat$pred[!dat$train] )
    value.B = mean((dat$Y[!dat$train]*I(fitted.B == dat$A[!dat$train]))/dat$pred[!dat$train])
    value.behavior =  mean(dat$Y[!dat$train]) 
    
    # estimating value function with MR
    #value function mr
    value.opt.mr = value.opt - mean((dat$Y.pred[!dat$train]*I(fitted.owl == dat$A[!dat$train]))/dat$pred[!dat$train] 
                                     - I( fitted.owl == 1)*dat$Y.A1[!dat$train] - I( fitted.owl == -1)*dat$Y.An1[!dat$train]) 
    
    value.A.mr = value.A - mean((dat$Y.pred[!dat$train]*I(fitted.A == dat$A[!dat$train]))/dat$pred[!dat$train] 
                                       - I(fitted.A == 1)*dat$Y.A1[!dat$train] - I(fitted.A == -1)*dat$Y.An1[!dat$train]) 
    
    value.B.mr = value.B - mean((dat$Y.pred[!dat$train]*I(fitted.B == dat$A[!dat$train]))/dat$pred[!dat$train] 
                                     - I(fitted.B == 1)*dat$Y.A1[!dat$train] - I(fitted.B == -1)*dat$Y.An1[!dat$train])
    
    #estimated variance
    var.opt.mr = var((dat$Y[!dat$train]*I(fitted.owl == dat$A[!dat$train]))/dat$pred[!dat$train] - 
                       ((dat$Y.pred[!dat$train]*I(fitted.owl == dat$A[!dat$train]))/dat$pred[!dat$train] - 
                         I( fitted.owl == 1)*dat$Y.A1[!dat$train] - I( fitted.owl == -1)*dat$Y.An1[!dat$train]))/sum( 1 - dat$train ) 
    
    var.A.mr = var((dat$Y[!dat$train]*I(fitted.A == dat$A[!dat$train]))/dat$pred[!dat$train] - 
                       ((dat$Y.pred[!dat$train]*I(fitted.A == dat$A[!dat$train]))/dat$pred[!dat$train] - 
                          I( fitted.A == 1)*dat$Y.A1[!dat$train] - I( fitted.A == -1)*dat$Y.An1[!dat$train]))/sum( 1 - dat$train ) 
    
    var.B.mr = var((dat$Y[!dat$train]*I(fitted.B == dat$A[!dat$train]))/dat$pred[!dat$train] - 
                       ((dat$Y.pred[!dat$train]*I(fitted.B == dat$A[!dat$train]))/dat$pred[!dat$train] - 
                          I( fitted.B == 1)*dat$Y.A1[!dat$train] - I( fitted.B == -1)*dat$Y.An1[!dat$train]))/sum( 1 - dat$train ) 
    
    #return the results 
    value.function =data.frame(est.ipw = c(value.opt, value.A, value.B), 
                               est.mr = c(value.opt.mr, value.A.mr, value.B.mr))
    colnames(value.function) = c("IPW estimator", "MR estimator")
    rownames(value.function) = c("Opt Trt.", "Trt A", "Trt B")
    
    var.value.function = c(var.opt.mr, var.A.mr, var.B.mr)
    names(var.value.function) = c("Var Opt.", "Var Trt A.", "Var Trt B.")
    
    return(list(model = mod.owl, 
                value.funcion = value.function, 
                variance = var.value.function ) )
}




