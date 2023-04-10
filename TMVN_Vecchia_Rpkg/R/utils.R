#' Find ordered nearest neighbors based on a correlation Matrix. Assuming the
#' absolute value of the correlation is monotonically decreasing with distance.
#' Returns an n X (m + 1) matrix similar to `GpGp::find_ordered_nn`.
#'
#' @param corrMat the correlation matrix
#' @param m the number of nearest neighbors
#' @return an n X (m + 1) matrix
#' @example
#' library(GpGp)
#' library(VeccTMVN)
#' set.seed(123)
#' d <- 3
#' n <- 100
#' locs <- matrix(runif(d * n), n, d)
#' covparms <- c(2, 0.01, 0)
#' cov_mat <- GpGp::matern15_isotropic(covparms, locs)
#' m <- 10
#' NNarray_test <- GpGp::find_ordered_nn(locs, m = m)
#' NNarray <- find_nn_corr(cov_mat, m)
#' cat("Number of mismatch is", sum(NNarray != NNarray_test, na.rm = T))
#'
find_nn_corr <- function(corrMat, m) {
  NN <- find_nn_corr_internal(corrMat, m) + 1
  n <- nrow(NN)
  if (any(NN[, 1] != 1:n)) {
    warning("Input corrMat is not positive definite at machine precision\n")
  }
  return(NN)
}
