library(truncnorm)


#' Sample from the proposal density in Idea V. Zero mean is assumed.
#' The `truncnorm` package uses accept-reject sampling and seems to be able to 
#'   sample from tail truncation although I haven't verified its accuracy in 
#'   tail sampling. 
#' Input:
#'   N - number of samples to draw
#'   veccCondMeanVarObj - contains information of the conditional mean 
#'     coefficient, the conditional variance, and the NN array of the Vecchia 
#'     approximation
#'   a - lower bound vector for TMVN
#'   b - upper bound vector for TMVN
#'   alpha - parameter of the proposal density
#'   beta - parameter of the proposal density
#' 
sample_5 <- function(N, veccCondMeanVarObj, a, b, alpha = rep(1, length(x)), 
                     beta = rep(0, length(x))){
  n <- length(a)
  X <- matrix(NA, n, N)
  for(i in 1 : n){
    ind <- veccCondMeanVarObj$nn[i, -1]
    sd <- alpha[i] * sqrt(veccCondMeanVarObj$cond_var) # scalar
    mu <- apply(veccCondMeanVarObj$cond_mean_coeff[i] * X[ind, ], 2, sum, 
                na.rm = T) + beta[i] * sd
    X[i, ] <- rtruncnorm(n = 1, a = a[i], b = b[i], mean = mu, sd = sd)
  }
  return(X)
}







