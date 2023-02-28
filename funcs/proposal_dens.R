source("utils.R")

#' Log-density function of the proposal density in Idea V. Zero mean is assumed.
#' Notice that currently, it seems that there isn't any R package that can 
#'   compute the log-pdf of truncated normal when the truncation is at the very 
#'   tail, e.g., 10. I tried the `truncnorm` package. Otherwise, we can use 
#'   existing packages to simplify the code here. But it's good to verify 
#'   correctness with the `truncnorm` package
#' Input:
#'   veccCondMeanVarObj - contains information of the conditional mean 
#'     coefficient, the conditional variance, and the NN array of the Vecchia 
#'     approximation
#'   x - an n X N matrix, where n is MVN dim and N is number of samples
#'   a - lower bound vector for TMVN
#'   b - upper bound vector for TMVN
#'   alpha - parameter of the proposal density
#'   beta - parameter of the proposal density
#' Return:
#'   a vector of length N, representing the log-pdfs
#' 
log_pdf_5 <- function(veccCondMeanVarObj, x, a, b, alpha = rep(1, length(x)), 
                       beta = rep(0, length(x))){
  logpdf <- 0
  n <- nrow(x)
  for(i in 1 : n){
    ind <- veccCondMeanVarObj$nn[i, -1]
    sd <- alpha[i] * sqrt(veccCondMeanVarObj$cond_var[i])
    mu <- apply(veccCondMeanVarObj$cond_mean_coeff[i, ] * x[ind, , drop = F], 
                2, sum, na.rm = T) + beta[i] * sd
    logpdf_i_nom <- dnorm(x = x[i, ], mean = mu, sd = sd, log = T)
    logpdf_i_denom <- lnNpr(a, b, mu, sd)
    logpdf <- logpdf + logpdf_i_nom - logpdf_i_denom
  }
  return(logpdf)
}