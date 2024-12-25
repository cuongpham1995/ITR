#' Fit an Outcome-Weighted Learning (OWL) Model
#'
#' This function fits an OWL model using a weighted SVM (wSVM) to optimize individual treatment rules.
#'
#' @param itr_formula A formula specifying the independent variables for treatment rule optimization.
#' @param prop_formula A formula specifying the independent variables for propensity score modeling.
#' @param trt A character string specifying the name of the treatment variable in the dataset.
#' @param dat A data frame containing the data for the model.
#' @param kernel A character string specifying the kernel type for the SVM (default is "linear").
#' @param scale Logical, whether to scale the data for the SVM (default is TRUE).
#' @param cross An integer specifying the number of folds for cross-validation in the SVM (default is 10).
#' @return A numeric vector of predicted treatment decisions.
#' @examples
#' # Example usage
#' dat <- data.frame(Y = rnorm(100), A = rbinom(100, 1, 0.5), X1 = rnorm(100), X2 = rnorm(100))
#' itr_formula <- ~ X1 + X2
#' prop_formula <- ~ X1 + X2
#' fit.owl(itr_formula, prop_formula, trt = "A", dat = dat)

fit.owl = function(itr_formula, prop_formula, trt, outcome, dat, 
                   na.rm = T, minimize = F, center.outcome = F, 
                   kernel = "linear", cross = 10, scale = T, true.value, ...){
     
    ## Error handling
    if (!inherits(itr_formula, "formula")) stop("itr_formula must be a formula object.")
    if (!inherits(prop_formula, "formula")) stop("prop_formula must be a formula object.")
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
  
    ## Data preparing
    names(dat)[names(dat) ==  trt] = "A"
    names(dat)[names(dat) == outcome] = "Y" 
    
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
    
    #convert it to a formula
    logistic.formula = paste( "as.factor(A)" , deparse(prop_formula), collapse = "")  
    
    logistic.A = glm(logistic.formula, data = dat, family = binomial(link = "logit"))
    
    pred.A = predict(logistic.A, type = "response", newdata = dat)
    dat$pred = ifelse(dat$A == 1, pred.A, 1 - pred.A)

    
    dat$weight = with(dat, Y/pred)
    dat$fitted.owl = NA
    
    
    ## calculate the weight
    
    dat$lab = sign(dat$weight)*dat$A
    
    
    formula_str <- paste("as.factor(lab)", deparse(itr_formula), collapse = "")
    formula.owl <- as.formula(formula_str)
    
    
    mod.owl = WeightSVM::wsvm(formula.owl, data = dat, weight = abs(dat$weight), kernel = "linear", 
                   cross = 10, scale = T) #OWL method
    
    
    ## MODEL FITTING
    fitted.owl = predict(mod.owl, newdata = dat)
    fitted.owl = as.numeric(as.character(fitted.owl))
  
    return( table(fitted.owl, true.value) )

    
}

## unit test
sim_dat = dat.gen.owl(23, nsample = 250)
colnames(sim_dat)[colnames(sim_dat) == "Y"] = "outcome"

sim_dat$X.1[12] = NA

fit.owl(itr_formula = ~ X.2 + X.3, prop_formula = ~ X.2 + X.3, trt = "A", outcome = "outcome", dat = sim_dat,
        true.value = sim_dat$opt.trt, center.outcome = T)


setdiff(c("X.1", "X.4"), names(sim_dat))
