#' Compute the conditional mean multiplier and the conditional variance
#'   under the Vecchia approximation.
#' 
vecc_cond_mean_var <- function(covMat, NNarray){
  n <- nrow(covMat)
  m <- ncol(NNarray)
  cond_mean_coeff <- matrix(0, n, m)
  cond_var <- rep(NA, n)
  cond_var[1] <- covMat[1, 1]
  A <- matrix(0, n, n)
  for(i in 2 : n){
    ind_cond <- NNarray[2 : min(i, (m + 1))]
    cov_mat_sub_inv <- solve(covMat[ind_cond, ind_cond])
    cov_vec_sub <- covMat[i, ind_cond]
    cond_var[i] <- covMat[i, i] - 
      t(cov_vec_sub) %*% cov_mat_sub_inv %*% cov_vec_sub
    cond_mean_coeff[i, min(i - 1, m)] <- t(cov_vec_sub) %*% cov_mat_sub_inv
    A[i, ind_cond] <- cond_mean_coeff[i, min(i - 1, m)]
  }
  return(list(cond_mean_coeff = cond_mean_coeff, cond_var = cond_var, 
              nn = NNarray, A = A))
}