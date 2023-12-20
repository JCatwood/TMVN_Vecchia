library(Matrix)
library(GpGp)

# Compute the conditional mean multiplier and the conditional variance
#   under the Vecchia approximation.
#
vecc_cond_mean_var <- function(covMat, NNarray) {
  n <- nrow(covMat)
  m <- ncol(NNarray) - 1
  cond_mean_coeff <- matrix(0, n, m)
  cond_var <- rep(NA, n)
  cond_var[1] <- covMat[1, 1]
  A <- matrix(0, n, n)
  for (i in 2:n) {
    ind_cond <- NNarray[i, 2:min(i, (m + 1))]
    cov_mat_sub_inv <- solve(covMat[ind_cond, ind_cond])
    cov_vec_sub <- covMat[i, ind_cond]
    cond_var[i] <- covMat[i, i] - as.numeric(
      t(cov_vec_sub) %*% cov_mat_sub_inv %*% cov_vec_sub
    )
    cond_mean_coeff[i, 1:min(i - 1, m)] <- t(cov_vec_sub) %*% cov_mat_sub_inv
    A[i, ind_cond] <- cond_mean_coeff[i, 1:min(i - 1, m)]
  }
  return(list(
    cond_mean_coeff = cond_mean_coeff, cond_var = cond_var,
    nn = NNarray, A = A
  ))
}


# Compute the conditional mean multiplier and the conditional variance
#   under the Vecchia approximation. Different from `vecc_cond_mean_var`, the
#   returned A is a sparse matrix.
#
vecc_cond_mean_var_sp <- function(NNarray, covMat = NULL, locs = NULL,
                                  covName = NULL, covParms = NULL) {
  n <- nrow(NNarray)
  m <- ncol(NNarray) - 1
  if (any(NNarray[, 1] != 1:n)) {
    stop("Unexpected NNarray: first col is not 1 : n\n")
  }
  cond_mean_coeff <- matrix(0, n, m)
  cond_var <- rep(NA, n)
  if (is.null(locs)) {
    use_locs <- F
    cov_func <- function(ind) {
      covMat[ind, ind, drop = F]
    }
  } else {
    use_locs <- T
    cov_func_GpGp <- utils:: getFromNamespace(covName, "GpGp")
    cov_func <- function(ind) {
      cov_func_GpGp(covParms, locs[ind, , drop = F])
    }
  }
  cond_var[1] <- cov_func(1)[1, 1]
  for (i in 2:n) {
    ind_cond <- NNarray[i, 1:min(i, (m + 1))]
    cov_mat_sub <- cov_func(ind_cond)
    cov_mat_sub_inv <- solve(cov_mat_sub[-1, -1])
    cov_vec_sub <- cov_mat_sub[1, -1]
    cond_var[i] <- cov_mat_sub[1, 1] - as.numeric(
      t(cov_vec_sub) %*% cov_mat_sub_inv %*% cov_vec_sub
    )
    cond_mean_coeff[i, 1:min(i - 1, m)] <- t(cov_vec_sub) %*% cov_mat_sub_inv
  }
  nnz_A <- (m + 1) * m / 2 + (n - m) * (m + 1) # num of non-zero
  A_row_inds <- rep(0, nnz_A)
  A_col_inds <- rep(0, nnz_A)
  A_vals <- rep(0, nnz_A)
  ind <- 1
  # iteration through rows
  for (i in 1:n) {
    nnz_A_i <- min(m + 1, i)
    A_row_inds[ind:(ind + nnz_A_i - 1)] <- i
    # first col ind is i
    A_col_inds[ind:(ind + nnz_A_i - 1)] <- NNarray[i, 1:nnz_A_i]
    # first A[i, i] should be zero
    A_vals[ind] <- 0
    if (i > 1) {
      A_vals[(ind + 1):(ind + nnz_A_i - 1)] <-
        cond_mean_coeff[i, 1:min(i - 1, m)]
    }
    ind <- ind + nnz_A_i
  }
  # create sparse A
  A <- Matrix::sparseMatrix(i = A_row_inds, j = A_col_inds, x = A_vals, dims = c(n, n))
  return(list(
    cond_mean_coeff = cond_mean_coeff, cond_var = cond_var,
    nn = NNarray, A = A
  ))
}


# TEST-------------------------------
#
# ## example spatial covariance matrices --------------------------------
# library(GpGp)
# library(VeccTMVN)
# set.seed(123)
# n1 <- 10
# n2 <- 10
# n <- n1*n2
# locs <- as.matrix(expand.grid((1 : n1) / n1, (1 : n2) / n2))
# covparms <- c(2, 0.3, 0)
# cov_mat <- matern15_isotropic(covparms, locs)
#
# ## ordering and NN --------------------------------
# m <- 30
# ord <- order_maxmin(locs)
# locs_ord <- locs[ord, , drop = FALSE]
# cov_mat_ord <- matern15_isotropic(covparms, locs_ord)
# NNarray <- find_ordered_nn(locs_ord, m = m)
#
# ## Vecchia approx --------------------------------
# vecc_cond_mean_var_obj <- vecc_cond_mean_var(cov_mat_ord, NNarray)
# vecc_cond_mean_var_obj_sp <- vecc_cond_mean_var_sp(NNarray, cov_mat_ord)
# vecc_cond_mean_var_obj_sp_locs <- vecc_cond_mean_var_sp(
#   NNarray, locs = locs_ord, covName = "matern15_isotropic", covParms = covparms)
# cat("Using covmat, F-norm of difference between dense and sparse A matrices is",
#     sqrt(sum((as.matrix(vecc_cond_mean_var_obj_sp$A) -
#                 vecc_cond_mean_var_obj$A)^2)), "\n")
# cat("Using locs, F-norm of difference between dense and sparse A matrices is",
#     sqrt(sum((as.matrix(vecc_cond_mean_var_obj_sp$A) -
#                 vecc_cond_mean_var_obj_sp_locs$A)^2)), "\n")
# cat("Cond var error", sqrt(sum((vecc_cond_mean_var_obj_sp_locs$cond_var -
#                                   vecc_cond_mean_var_obj$cond_var)^2)), "\n")
