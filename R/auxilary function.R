#' Extract Coefficients and Slope from an SVM Model
#'
#' This function extracts the coefficients and intercept from a fitted SVM model 
#' and optionally calculates the slope for a linear decision boundary.
#'
#' @param mod A fitted SVM model object (e.g., from the `e1071` package).
#' @param type A character string indicating the type of output. 
#'   `"line"` returns the coefficients and intercept as a vector. 
#'   `"slope"` returns the slope of the decision boundary, normalized by the coefficient of `x2`.
#'   Default is `"line"`.
#' @return A numeric vector containing either:
#'   - The coefficients and intercept (`type = "line"`)
#'   - The slope and intercept in terms of `x2` (`type = "slope"`)
#' @examples
#' # Example with a dummy SVM model
#' # Note: Replace with a real fitted SVM model
#' # svm_model <- e1071::svm(Y ~ ., data = dat, kernel = "linear")
#' # coef <- svm.coef(svm_model, type = "line")
#' # slope <- svm.coef(svm_model, type = "slope")
#' @export

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



#' Generate Synthetic Data for Outcome-Weighted Learning (OWL)
#'
#' This function generates synthetic data for testing outcome-weighted learning (OWL) methods.
#' The data includes features, treatment assignments, and potential outcomes.
#'
#' @param seed An integer seed for random number generation to ensure reproducibility.
#' @param nsample The number of samples to generate.
#' @return A data frame containing:
#'   - `Y`: The observed outcome.
#'   - `A`: The treatment assignment (-1 or 1).
#'   - `X.1`, `X.2`, ..., `X.n`: The covariates/features.
#'   - `Y1`: The potential outcome if `A = 1`.
#'   - `Yn1`: The potential outcome if `A = -1`.
#'   - `opt.trt`: The optimal treatment based on the potential outcomes.
#' @details
#' The synthetic data assumes a model where the outcome depends on the covariates `X` 
#' and treatment assignment `A` through a specified mean structure. The covariance of 
#' the outcome is controlled by a diagonal covariance matrix.
#' @examples
#' # Generate a dataset with 100 samples
#' dat <- dat.gen.owl(seed = 123, nsample = 100)
#' head(dat)
#' @export

dat.gen.owl =function(seed, nsample){
  set.seed(seed)
  require(mvtnorm)
  require(dplyr)
  
  X <- matrix(runif(n = 3 * nsample, min = -1, max = 1), ncol = 3)
  sigma.mat = diag(nsample)*1
  A = as.matrix(2*rbinom(nsample, size = 1, prob = 0.5) - 1)
  
  #Y = rmvnorm(n =1, mean = 1 + 2*X[,1] + X[,2] + 0.5*X[,3] + 0.442*(1 - X[,1] - X[,2])*A[,1], sigma = sigma.mat) %>% t()
  Y = rmvnorm(n =1, mean = 1 +  (0.5 - 0.5*X[,1] - X[,2])*A[,1], sigma = sigma.mat) %>% t()
  
  #potential values
  Y1 =as.matrix( 1 +  (0.5 - 0.5*X[,1] - X[,2])*1)
  Yn1 = as.matrix( 1+    (0.5 - 0.5*X[,1] - X[,2])*(-1)) 
  
  #table(Y1 > Yn1)
  
  dat = cbind(Y, A, X, Y1, Yn1) %>% as.data.frame() 
  colnames(dat) = c("Y", "A", paste("X", seq(1:ncol(X)), sep = "."), "Y1", "Yn1")
  dat$opt.trt = (dat$Y1 > dat$Yn1)*2 - 1 
  return(dat)
}




