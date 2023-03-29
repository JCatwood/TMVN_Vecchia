library(Matrix)
library(VeccTMVN)


H11_mul <- function(veccCondMeanVarObj, dPsi, D, x){
  as.vector(t(dPsi / D / D * as.vector(veccCondMeanVarObj$A %*% x)) %*%
              (veccCondMeanVarObj$A))
}


H12_mul <- function(veccCondMeanVarObj, dPsi, D, x){
  as.vector(t((x + dPsi * x) / D) %*% veccCondMeanVarObj$A) - x / D
}


H21_mul <- function(veccCondMeanVarObj, dPsi, D, x){
  as.vector(veccCondMeanVarObj$A %*% x) / D * (1 + dPsi) - x / D
}


H22_mul <- function(veccCondMeanVarObj, dPsi, D, x){
  ((1 + dPsi) * x)
}


HInv11_mul <- function(veccCondMeanVarObj, dPsi, D, V, x){
  - solve(t(V), solve(V, x))
}


HInv12_mul <- function(veccCondMeanVarObj, dPsi, D, V, x){
  x <- x / (1 + dPsi)
  x <- H12_mul(veccCondMeanVarObj, dPsi, D, x)
  solve(t(V), solve(V, x))
}


HInv21_mul <- function(veccCondMeanVarObj, dPsi, D, V, x){
  x <- solve(t(V), solve(V, x))
  x <- H21_mul(veccCondMeanVarObj, dPsi, D, x)
  x / (1 + dPsi)
}


HInv22_mul <- function(veccCondMeanVarObj, dPsi, D, V, x){
  x <- x / (1 + dPsi)
  y <- x
  x <- H12_mul(veccCondMeanVarObj, dPsi, D, x)
  x <- - solve(t(V), solve(V, x))
  x <- H21_mul(veccCondMeanVarObj, dPsi, D, x)
  x <- x / (1 + dPsi)
  x + y
}


# Gradient, Newton step, and Hessian multiply gradient of the psi function in
#   Idea V, computed using sparse A
grad_jacprod_jacsolv_idea5 <- function(xAndBeta, veccCondMeanVarObj, a, b) {
  n <- length(a)
  x <- rep(0, n)
  beta <- rep(0, n)
  x[1 : (n - 1)] <- xAndBeta[1 : (n - 1)]
  beta[1 : (n - 1)] <- xAndBeta[n : (2 * n - 2)]
  D <- sqrt(veccCondMeanVarObj$cond_var)
  mu_c <- as.vector(veccCondMeanVarObj$A %*% x)
  a_tilde_shift <- (a - mu_c) / D - beta
  b_tilde_shift <- (b - mu_c) / D - beta
  log_diff_cdf <- TruncatedNormal::lnNpr(a_tilde_shift, b_tilde_shift)
  pl <- exp(-0.5 * a_tilde_shift ^ 2 - log_diff_cdf) / sqrt(2 * pi)
  pu <- exp(-0.5 * b_tilde_shift ^ 2 - log_diff_cdf) / sqrt(2 * pi)
  Psi <- pl - pu

  dpsi_dx_padded <- as.vector(t(veccCondMeanVarObj$A) %*%
                                (beta / D + Psi / D)) - beta / D  # length n
  dpsi_dbeta_padded <- beta - (x - mu_c) / D + Psi  # length n
  dpsi_dx_padded[n] <- 0
  dpsi_dbeta_padded[n] <- 0

  a_tilde_shift[is.infinite(a_tilde_shift)] <- 0
  b_tilde_shift[is.infinite(b_tilde_shift)] <- 0
  dPsi <- (-Psi ^ 2) + a_tilde_shift * pl - b_tilde_shift * pu

  dpsi_dx_dx <- t(veccCondMeanVarObj$A[, -n]) %*%
    (dPsi / D / D * veccCondMeanVarObj$A[, -n])
  dpsi_dx_dbeta <- t(dPsi / D * veccCondMeanVarObj$A)[-n, -n] -
    t((diag(rep(1, n - 1)) - veccCondMeanVarObj$A[-n, -n]) / D[-n])
  dpsi_dbeta_dbeta <- diag(1 + dPsi[-n])
  H <- rbind(cbind(dpsi_dx_dx, dpsi_dx_dbeta),
        cbind(t(dpsi_dx_dbeta), dpsi_dbeta_dbeta))

  # compute Jac \cdot grad -------------------------------------------
  H11_dpsi_dx <- H11_mul(veccCondMeanVarObj, dPsi, D, dpsi_dx_padded)[-n]
  H12_dpsi_dbeta <- H12_mul(veccCondMeanVarObj, dPsi, D, dpsi_dbeta_padded)[-n]
  H21_dpsi_dx <- H21_mul(veccCondMeanVarObj, dPsi, D, dpsi_dx_padded)[-n]
  H22_dpsi_dbeta <- H22_mul(veccCondMeanVarObj, dPsi, D, dpsi_dbeta_padded)[-n]

  # compute V -------------------------------------------
  # U = t(D^{-1} A), upper tri
  U <- t(veccCondMeanVarObj$A / D)
  # First entry of each col in U should be diag entry ############### TBD
  if(any(U@i[(U@p[-(n + 1)] + 1)] != c(1 : n)))
    stop("Error U is not lower-triangular\n")
  # U = D^{-1} A - D^{-1}
  U@x[(U@p[-(n + 1)] + 1)] <- - 1 / D
  # "precision" matrix P \approx U U^{\top} as the triangular part of P is
  #   assumed to have the same sparsity as U, only upper-tri of P is stored
  P_row_ind <- U@i
  P_col_ind <- rep(1 : n, diff(U@p))
  P_vals <- sp_mat_mul_query(P_row_ind, P_col_ind, U@i, U@p, U@x)
  # P = U U^{\top} + D^{-1} (I + dPsi)^{-1} D^{-1} - D^{-2}
  P_diag_ind <- P_col_ind == P_row_ind
  P_vals[P_diag_ind] <- P_vals[P_diag_ind] + (1 / (1 + dPsi) - 1) / (D^2)
  P <- sparseMatrix(i = P_row_ind, j = P_col_ind, x = P_vals, dims = c(n, n),
                    triangular = T)
  # call ic0 from GPVecchia
  V <- GPvecchia::ichol(P)  # V should be upper-tri
  if (!(Matrix::isTriangular(V) &&
        attr(Matrix::isTriangular(V), "kind") == "U"))
    stop("Returned V is not an upper-tri matrix\n")

  # compute Jac^{-1} \cdot grad -------------------------------------------
  HInv11_dpsi_dx <- HInv11_mul(veccCondMeanVarObj, dPsi, D, V,
                                   dpsi_dx_padded)[-n]
  HInv12_dpsi_dbeta <- HInv12_mul(veccCondMeanVarObj, dPsi, D, V,
                                      dpsi_dbeta_padded)[-n]
  HInv21_dpsi_dx <- HInv21_mul(veccCondMeanVarObj, dPsi, D, V,
                                   dpsi_dx_padded)[-n]
  HInv22_dpsi_dbeta <- HInv22_mul(veccCondMeanVarObj, dPsi, D, V,
                                         dpsi_dbeta_padded)[-n]
}


# TEST-------------------------------


## example MVN probabilities --------------------------------
library(GpGp)
source("inv_chol.R")
source("vecc_cond_mean_var.R")
source("gradpsi.R")
set.seed(123)
n1 <- 10
n2 <- 10
n <- n1*n2
locs <- as.matrix(expand.grid((1 : n1) / n1, (1 : n2) / n2))
covparms <- c(2, 0.3, 0)
cov_mat <- matern15_isotropic(covparms, locs)
a_list <- list(rep(-Inf, n), rep(-1, n), -runif(n) * 2)
b_list <- list(rep(-2, n), rep(1, n), runif(n) * 2)

## ordering and NN --------------------------------
m <- 30
ord <- order_maxmin(locs)
locs_ord <- locs[ord, , drop = FALSE]
cov_mat_ord <- matern15_isotropic(covparms, locs_ord)
a_list_ord <- lapply(a_list, function(x){x[ord]})
b_list_ord <- lapply(b_list, function(x){x[ord]})
NNarray <- find_ordered_nn(locs_ord, m = m)

## Vecchia approx --------------------------------
U <- get_sp_inv_chol(cov_mat_ord, NNarray)
cov_mat_Vecc <- solve(U %*% t(U))
vecc_cond_mean_var_obj <- vecc_cond_mean_var_sp(cov_mat_ord, NNarray)

## Compare dpsi ------------------------------
for(i in 1 : length(a_list_ord)){
  a_ord <- a_list_ord[[i]]
  b_ord <- b_list_ord[[i]]
  x0 <- runif(2 * n - 2)
  grad_ds <- grad_idea5(x0, vecc_cond_mean_var_obj, a_ord, b_ord)
  jac_ds <- jac_idea5(x0, vecc_cond_mean_var_obj, a_ord, b_ord)
  solve_obj <- grad_jacprod_jacsolv_idea5(x0, vecc_cond_mean_var_obj,
                                          a_ord, b_ord)
}
