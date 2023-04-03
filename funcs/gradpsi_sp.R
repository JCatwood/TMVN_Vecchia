library(Matrix)
library(VeccTMVN)


H11_mul <- function(veccCondMeanVarObj, dPsi, D, x) {
  as.vector(t(dPsi / D / D * as.vector(veccCondMeanVarObj$A %*% x)) %*%
    (veccCondMeanVarObj$A))
}


H12_mul <- function(veccCondMeanVarObj, dPsi, D, x) {
  as.vector(t((x + dPsi * x) / D) %*% veccCondMeanVarObj$A) - x / D
}


H21_mul <- function(veccCondMeanVarObj, dPsi, D, x) {
  as.vector(veccCondMeanVarObj$A %*% x) / D * (1 + dPsi) - x / D
}


H22_mul <- function(veccCondMeanVarObj, dPsi, D, x) {
  ((1 + dPsi) * x)
}


HInv11_mul <- function(veccCondMeanVarObj, dPsi, D, V, x) {
  -as.vector(solve(V, solve(t(V), x)))
}


HInv12_mul <- function(veccCondMeanVarObj, dPsi, D, V, x) {
  x <- x / (1 + dPsi)
  x <- H12_mul(veccCondMeanVarObj, dPsi, D, x)
  as.vector(solve(V, solve(t(V), x)))
}


HInv21_mul <- function(veccCondMeanVarObj, dPsi, D, V, x) {
  x <- solve(V, solve(t(V), x))
  x <- H21_mul(veccCondMeanVarObj, dPsi, D, x)
  as.vector(x / (1 + dPsi))
}


HInv22_mul <- function(veccCondMeanVarObj, dPsi, D, V, x) {
  x <- x / (1 + dPsi)
  y <- x
  x <- H12_mul(veccCondMeanVarObj, dPsi, D, x)
  x <- -solve(V, solve(t(V), x))
  x <- H21_mul(veccCondMeanVarObj, dPsi, D, x)
  x <- x / (1 + dPsi)
  as.vector(x + y)
}


# Gradient, Newton step, and Hessian multiply gradient of the psi function in
#   Idea V, computed using sparse A. For this function, `xAndBeta` should be a
#     vector of length 2 * n and all the returned vectors (matrices) are of
#     dimension 2 * n. In other words, beta_n and x_n are taken into
#     consideration. Meanwhile, based on Idea 5, it is obvious that beta_n = 0
#     and we don't need to know the value of x_n anyways
grad_idea5_sp <- function(xAndBeta, veccCondMeanVarObj, a, b) {
  n <- length(a)
  x <- xAndBeta[1:n]
  beta <- xAndBeta[(n + 1):(2 * n)]
  D <- sqrt(veccCondMeanVarObj$cond_var)
  mu_c <- as.vector(veccCondMeanVarObj$A %*% x)
  a_tilde_shift <- (a - mu_c) / D - beta
  b_tilde_shift <- (b - mu_c) / D - beta
  log_diff_cdf <- TruncatedNormal::lnNpr(a_tilde_shift, b_tilde_shift)
  pl <- exp(-0.5 * a_tilde_shift^2 - log_diff_cdf) / sqrt(2 * pi)
  pu <- exp(-0.5 * b_tilde_shift^2 - log_diff_cdf) / sqrt(2 * pi)
  Psi <- pl - pu
  # compute grad ------------------------------------------------
  dpsi_dx <- as.vector(t(veccCondMeanVarObj$A) %*%
    (beta / D + Psi / D)) - beta / D
  dpsi_dbeta <- beta - (x - mu_c) / D + Psi
  return(c(dpsi_dx, dpsi_dbeta))
}


# Gradient, Newton step, and Hessian multiply gradient of the psi function in
#   Idea V, computed using sparse A. For this function, `xAndBeta` should be a
#     vector of length 2 * n and all the returned vectors (matrices) are of
#     dimension 2 * n. In other words, beta_n and x_n are taken into
#     consideration. Meanwhile, based on Idea 5, it is obvious that beta_n = 0
#     and we don't need to know the value of x_n anyways
grad_jacprod_jacsolv_idea5 <- function(xAndBeta, veccCondMeanVarObj, a, b,
                                       retJac = F) {
  n <- length(a)
  x <- xAndBeta[1:n]
  beta <- xAndBeta[(n + 1):(2 * n)]
  D <- sqrt(veccCondMeanVarObj$cond_var)
  mu_c <- as.vector(veccCondMeanVarObj$A %*% x)
  a_tilde_shift <- (a - mu_c) / D - beta
  b_tilde_shift <- (b - mu_c) / D - beta
  log_diff_cdf <- TruncatedNormal::lnNpr(a_tilde_shift, b_tilde_shift)
  pl <- exp(-0.5 * a_tilde_shift^2 - log_diff_cdf) / sqrt(2 * pi)
  pu <- exp(-0.5 * b_tilde_shift^2 - log_diff_cdf) / sqrt(2 * pi)
  Psi <- pl - pu
  # compute grad ------------------------------------------------
  dpsi_dx <- as.vector(t(veccCondMeanVarObj$A) %*%
    (beta / D + Psi / D)) - beta / D
  dpsi_dbeta <- beta - (x - mu_c) / D + Psi
  # compute dPsi ------------------------------------------------
  a_tilde_shift[is.infinite(a_tilde_shift)] <- 0
  b_tilde_shift[is.infinite(b_tilde_shift)] <- 0
  dPsi <- (-Psi^2) + a_tilde_shift * pl - b_tilde_shift * pu
  # compute Hessian (usually for testing) -------------------------------
  if (retJac) {
    dpsi_dx_dx <- t(veccCondMeanVarObj$A) %*%
      (dPsi / D / D * veccCondMeanVarObj$A)
    dpsi_dx_dbeta <- t(dPsi / D * veccCondMeanVarObj$A) -
      t((diag(rep(1, n)) - veccCondMeanVarObj$A) / D)
    dpsi_dbeta_dbeta <- diag(1 + dPsi)
    H <- rbind(
      cbind(dpsi_dx_dx, dpsi_dx_dbeta),
      cbind(t(dpsi_dx_dbeta), dpsi_dbeta_dbeta)
    )
  }
  # compute Hessian \cdot grad -------------------------------------------
  H11_dpsi_dx <- H11_mul(veccCondMeanVarObj, dPsi, D, dpsi_dx)
  H12_dpsi_dbeta <- H12_mul(veccCondMeanVarObj, dPsi, D, dpsi_dbeta)
  H21_dpsi_dx <- H21_mul(veccCondMeanVarObj, dPsi, D, dpsi_dx)
  H22_dpsi_dbeta <- H22_mul(veccCondMeanVarObj, dPsi, D, dpsi_dbeta)
  # compute V -------------------------------------------
  # L = D^{-1} A, lower tri
  L <- veccCondMeanVarObj$A / D
  # First entry of each col in L should be diag entry
  if (any(L@i[(L@p[-(n + 1)] + 1)] != c(0:(n - 1)))) {
    stop("Error L is not lower-triangular\n")
  }
  # L = D^{-1} A - D^{-1}
  L@x[(L@p[-(n + 1)] + 1)] <- -1 / D
  # "precision" matrix P \approx L^{\top} L. The triangular part of P is
  #   assumed to have the same sparsity as L, only upper-tri of P is stored
  P_col_ind <- L@i + 1
  P_row_ind <- rep(1:n, diff(L@p))
  P_vals <- sp_mat_mul_query(P_row_ind, P_col_ind, L@i, L@p, L@x)
  # P = L^{\top} L + D^{-1} (I + dPsi)^{-1} D^{-1} - D^{-2}
  P_diag_ind <- P_col_ind == P_row_ind
  P_vals[P_diag_ind] <- P_vals[P_diag_ind] + (1 / (1 + dPsi) - 1) / (D^2)
  P <- sparseMatrix(
    i = P_row_ind, j = P_col_ind, x = P_vals, dims = c(n, n),
    symmetric = T
  )
  # call ic0 from GPVecchia, this changes P matrix!
  V <- GPvecchia::ichol(P)
  # V should be upper-tri
  if (!(Matrix::isTriangular(V) &&
    attr(Matrix::isTriangular(V), "kind") == "U")) {
    stop("Returned V is not an upper-tri matrix\n")
  }
  # compute Jac^{-1} \cdot grad -------------------------------------------
  HInv11_dpsi_dx <- HInv11_mul(
    veccCondMeanVarObj, dPsi, D, V,
    dpsi_dx
  )
  HInv12_dpsi_dbeta <- HInv12_mul(
    veccCondMeanVarObj, dPsi, D, V,
    dpsi_dbeta
  )
  HInv21_dpsi_dx <- HInv21_mul(
    veccCondMeanVarObj, dPsi, D, V,
    dpsi_dx
  )
  HInv22_dpsi_dbeta <- HInv22_mul(
    veccCondMeanVarObj, dPsi, D, V,
    dpsi_dbeta
  )

  rslt <- list(
    grad = c(dpsi_dx, dpsi_dbeta),
    jac_grad = c(
      H11_dpsi_dx + H12_dpsi_dbeta,
      H21_dpsi_dx + H22_dpsi_dbeta
    ),
    jac_inv_grad = c(
      HInv11_dpsi_dx + HInv12_dpsi_dbeta,
      HInv21_dpsi_dx + HInv22_dpsi_dbeta
    )
  )
  if (retJac) {
    rslt[["jac"]] <- H
  }
  return(rslt)
}


# # TEST-------------------------------
#
#
# ## example MVN probabilities --------------------------------
# library(GpGp)
# source("inv_chol.R")
# source("vecc_cond_mean_var.R")
# source("gradpsi.R")
# set.seed(123)
# n1 <- 10
# n2 <- 10
# n <- n1 * n2
# locs <- as.matrix(expand.grid((1:n1) / n1, (1:n2) / n2))
# covparms <- c(2, 0.3, 0)
# cov_mat <- matern15_isotropic(covparms, locs)
# a_list <- list(rep(-Inf, n), rep(-1, n), -runif(n) * 2)
# b_list <- list(rep(-2, n), rep(1, n), runif(n) * 2)
#
# ## ordering and NN --------------------------------
# m <- 30
# ord <- order_maxmin(locs)
# locs_ord <- locs[ord, , drop = FALSE]
# cov_mat_ord <- matern15_isotropic(covparms, locs_ord)
# a_list_ord <- lapply(a_list, function(x) {
#   x[ord]
# })
# b_list_ord <- lapply(b_list, function(x) {
#   x[ord]
# })
# NNarray <- find_ordered_nn(locs_ord, m = m)
#
# ## Vecchia approx --------------------------------
# U <- get_sp_inv_chol(cov_mat_ord, NNarray)
# cov_mat_Vecc <- solve(U %*% t(U))
# vecc_cond_mean_var_obj <- vecc_cond_mean_var_sp(cov_mat_ord, NNarray)
#
# ## Compare dpsi ------------------------------
# for (i in 1:length(a_list_ord)) {
#   a_ord <- a_list_ord[[i]]
#   b_ord <- b_list_ord[[i]]
#   x0_padded <- runif(2 * n)
#   x0_padded[2 * n] <- 0
#   x0 <- x0_padded[-c(n, 2 * n)]
#   grad_ds <- grad_idea5(x0, vecc_cond_mean_var_obj, a_ord, b_ord)
#   jac_ds <- jac_idea5(x0, vecc_cond_mean_var_obj, a_ord, b_ord)
#   solve_obj <- grad_jacprod_jacsolv_idea5(x0_padded, vecc_cond_mean_var_obj,
#     a_ord, b_ord,
#     retJac = T
#   )
#   cat("grad error: ", sum(abs(grad_ds - solve_obj$grad[-c(n, 2 * n)])), "\n")
#   cat(
#     "jac error: ", sum(abs(jac_ds - solve_obj$jac[-c(n, 2 * n), -c(n, 2 * n)])),
#     "\n"
#   )
#   H_inv <- solve(solve_obj$jac)
#   jac_grad_test <- as.vector(solve_obj$jac %*% solve_obj$grad)
#   jac_inv_grad_test <- as.vector(H_inv %*% solve_obj$grad)
#   cat("jac_grad error: ", sum(abs(jac_grad_test - solve_obj$jac_grad)), "\n")
#   cat("jac_inv_grad error: ", sum(abs(jac_inv_grad_test - solve_obj$jac_inv_grad)), "\n")
# }
