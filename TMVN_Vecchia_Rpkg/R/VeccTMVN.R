#' VeccTMVN
#' 
#' Exponential tilting method made scalable through the Vecchia approximation.
#' References are Botev (2017) and the TruncatedNormal R package
#' 
#' @importFrom Rcpp evalCpp
#' @useDynLib mypackage, .registration=TRUE
#' @name VeccTMVN

#' @examples
#' library(GpGp)
#' library(VeccTMVN)
#' set.seed(123)
#' n1 <- 10
#' n2 <- 10
#' n <- n1 * n2
#' locs <- as.matrix(expand.grid((1:n1) / n1, (1:n2) / n2))
#' covparms <- c(2, 0.3, 0)
#' cov_mat <- GpGp::matern15_isotropic(covparms, locs)
#' a <- rep(-Inf, n)
#' b <- -runif(n) * 2
#' m <- 30
#' est_Vecc <- VeccTMVN::pmvn(
#'   a, b, 0, locs,
#'   covName = "matern15_isotropic",
#'   covParms = covparms, m = m, verbose = F
#' )
#'
## usethis namespace: start
## usethis namespace: end
