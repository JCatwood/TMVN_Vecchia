get_sp_inv_chol <- function(covMat, NNarray){
  n <- nrow(covMat)
  inv_chol <- matrix(0, n, n)
  for(i in 1 : n){
    idx <- sort(NNarray[i, ])
    idx <- idx[!is.na(idx)]
    nnz <- length(idx)
    tmp_mat <- solve(covMat[idx, idx])
    inv_chol[idx, i] <- tmp_mat[, nnz] / sqrt(tmp_mat[nnz, nnz])
  }
  return(inv_chol)
}
