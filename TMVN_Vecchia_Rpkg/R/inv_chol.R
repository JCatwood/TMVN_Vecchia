#' Get the inverse upper Cholesky factor under the Vecchia approximation
#'
#' @param covMat the covariance matrix
#' @param NNarray n X (m + 1) matrix representing the nearest neighbor indices
#' among previous observations. This is typically the return of
#' GpGp::find_ordered_nn
#' @return upper Cholesky of the inverse of `covMat`
#' @examples
#' library(GpGp)
#' n1 <- 10
#' n2 <- 10
#' n <- n1 * n2
#' locs <- as.matrix(expand.grid((1:n1) / n1, (1:n2) / n2))
#' covparms <- c(2, 0.3, 0)
#' cov_mat <- GpGp::matern15_isotropic(covparms, locs)
#' m <- 30
#' NNarray <- GpGp::find_ordered_nn(locs, m = m)
#' # Vecchia approx --------------------------------
#' U_Vecc <- get_sp_inv_chol(cov_mat, NNarray)
#' U <- solve(chol(cov_mat))
#' cat("Frobenius norm of the difference is", sqrt(sum((U - U_Vecc)^2)))
#'
#' @export
get_sp_inv_chol <- function(covMat, NNarray) {
  n <- nrow(covMat)
  inv_chol <- matrix(0, n, n)
  for (i in 1:n) {
    idx <- sort(NNarray[i, ])
    idx <- idx[!is.na(idx)]
    nnz <- length(idx)
    tmp_mat <- solve(covMat[idx, idx])
    inv_chol[idx, i] <- tmp_mat[, nnz] / sqrt(tmp_mat[nnz, nnz])
  }
  return(inv_chol)
}
