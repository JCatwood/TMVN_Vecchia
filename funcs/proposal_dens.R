#' Log-density function of the proposal density in Idea V
#' Input:
#'   veccCondMeanVarObj - contains information of the conditional mean 
#'     coefficient, the conditional variance, and the NN array of the Vecchia 
#'     approximation
#'   x - a sample vector 
#'   a - lower bound vector for TMVN
#'   b - upper bound vector for TMVN
#'   alpha - parameter of the proposal density
#'   beta - parameter of the proposal density
#' 
log_pdf_5 <- function(veccCondMeanVarObj, x, a, b, alpha = rep(1, length(x)), 
                       beta = rep(0, length(x))){
  
}