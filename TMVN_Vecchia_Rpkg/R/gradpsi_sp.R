library(Matrix)


H11_mul <- function(veccCondMeanVarObj, dPsi, D, x) {
  as.vector(Matrix::t(dPsi / D / D * as.vector(veccCondMeanVarObj$A %*% x)) %*%
    (veccCondMeanVarObj$A))
}


H12_mul <- function(veccCondMeanVarObj, dPsi, D, x) {
  as.vector(Matrix::t((x + dPsi * x) / D) %*% veccCondMeanVarObj$A) - x / D
}


H21_mul <- function(veccCondMeanVarObj, dPsi, D, x) {
  as.vector(veccCondMeanVarObj$A %*% x) / D * (1 + dPsi) - x / D
}


H22_mul <- function(veccCondMeanVarObj, dPsi, D, x) {
  ((1 + dPsi) * x)
}


HInv11_mul <- function(veccCondMeanVarObj, dPsi, D, V, S, x) {
  -as.vector(solve(V, solve(Matrix::t(V), x))) / S
}


HInv12_mul <- function(veccCondMeanVarObj, dPsi, D, V, S, x) {
  x <- x / (1 + dPsi)
  x <- H12_mul(veccCondMeanVarObj, dPsi, D, x)
  as.vector(solve(V, solve(Matrix::t(V), x))) / S
}


HInv21_mul <- function(veccCondMeanVarObj, dPsi, D, V, S, x) {
  x <- solve(V, solve(Matrix::t(V), x)) / S
  x <- H21_mul(veccCondMeanVarObj, dPsi, D, x)
  as.vector(x / (1 + dPsi))
}


HInv22_mul <- function(veccCondMeanVarObj, dPsi, D, V, S, x) {
  x <- x / (1 + dPsi)
  y <- x
  x <- H12_mul(veccCondMeanVarObj, dPsi, D, x)
  x <- -solve(V, solve(Matrix::t(V), x)) / S
  x <- H21_mul(veccCondMeanVarObj, dPsi, D, x)
  x <- x / (1 + dPsi)
  as.vector(x + y)
}


# Gradient, Newton step, and Hessian multiply gradient of the psi function in
# Idea V, computed using sparse A. For this function, `xAndBeta` should be a
# vector of length 2 * n and all the returned vectors (matrices) are of
# dimension 2 * n. In other words, beta_n and x_n are taken into
# consideration. Meanwhile, based on Idea 5, it is obvious that beta_n = 0
# and we don't need to know the value of x_n anyways
grad_jacprod_jacsolv_idea5 <- function(xAndBeta, veccCondMeanVarObj, a, b,
                                       retJac = F, retProd = T, retSolv = T,
                                       VAdj = T, verbose = T, mu = rep(0, length(a))) {
  n <- length(a)
  x <- xAndBeta[1:n]
  beta <- xAndBeta[(n + 1):(2 * n)]
  D <- sqrt(veccCondMeanVarObj$cond_var)
  mu_c <- as.vector(veccCondMeanVarObj$A %*% (x - mu)) + mu
  a_tilde_shift <- (a - mu_c) / D - beta
  b_tilde_shift <- (b - mu_c) / D - beta
  log_diff_cdf <- TruncatedNormal::lnNpr(a_tilde_shift, b_tilde_shift)
  pl <- exp(-0.5 * a_tilde_shift^2 - log_diff_cdf) / sqrt(2 * pi)
  pu <- exp(-0.5 * b_tilde_shift^2 - log_diff_cdf) / sqrt(2 * pi)
  Psi <- pl - pu
  # compute grad ------------------------------------------------
  dpsi_dx <- as.vector(Matrix::t(veccCondMeanVarObj$A) %*%
    (beta / D + Psi / D)) - beta / D
  dpsi_dbeta <- beta - (x - mu_c) / D + Psi
  # compute dPsi ------------------------------------------------
  if (retProd | retSolv | retJac) {
    a_tilde_shift[is.infinite(a_tilde_shift)] <- 0
    b_tilde_shift[is.infinite(b_tilde_shift)] <- 0
    dPsi <- (-Psi^2) + a_tilde_shift * pl - b_tilde_shift * pu
  }
  # compute Hessian (usually for testing) -------------------------------
  if (retJac) {
    dpsi_dx_dx <- Matrix::t(veccCondMeanVarObj$A) %*%
      (dPsi / D / D * veccCondMeanVarObj$A)
    dpsi_dx_dbeta <- Matrix::t(dPsi / D * veccCondMeanVarObj$A) -
      Matrix::t((diag(rep(1, n)) - veccCondMeanVarObj$A) / D)
    dpsi_dbeta_dbeta <- diag(1 + dPsi)
    H <- rbind(
      cbind(dpsi_dx_dx, dpsi_dx_dbeta),
      cbind(Matrix::t(dpsi_dx_dbeta), dpsi_dbeta_dbeta)
    )
  }
  # compute Hessian \cdot grad -------------------------------------------
  if (retProd | (retSolv & VAdj)) {
    H11_dpsi_dx <- H11_mul(veccCondMeanVarObj, dPsi, D, dpsi_dx)
    H12_dpsi_dbeta <- H12_mul(veccCondMeanVarObj, dPsi, D, dpsi_dbeta)
    H21_dpsi_dx <- H21_mul(veccCondMeanVarObj, dPsi, D, dpsi_dx)
    H22_dpsi_dbeta <- H22_mul(veccCondMeanVarObj, dPsi, D, dpsi_dbeta)
  }
  # compute Hessian^{-1} \cdot grad -------------------------------------------
  if (retSolv) {
    # compute V -------------------------------------------
    ichol_succ <- T
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
    P <- Matrix::sparseMatrix(
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
    if (any(is.na(V)) |
      min(diag(V)) / max(diag(V)) < sqrt(.Machine$double.eps)) {
      if (verbose) {
        cat("ichol failed, using diagonal adjustment on L instead\n")
      }
      ichol_succ <- F
    }
    if (!ichol_succ) {
      # increase the absolute values of the diag coeffs of L and use it as V s.t.
      #   diag(V^{\top} V) = L^{\top} L + D^{-1} (I + dPsi)^{-1} D^{-1} - D^{-2}
      L_row_ind <- L@i + 1
      L_col_ind <- rep(1:n, diff(L@p))
      L_diag_ind <- L_row_ind == L_col_ind
      V_vals <- L@x
      V_vals[L_diag_ind] <- sign(V_vals[L_diag_ind]) *
        sqrt(V_vals[L_diag_ind]^2 + (1 / (1 + dPsi) - 1) / (D^2))
      V <- Matrix::sparseMatrix(
        i = L_row_ind, j = L_col_ind, x = V_vals, dims = c(n, n),
        triangular = T
      )
      if (!(Matrix::isTriangular(V) &&
        attr(Matrix::isTriangular(V), "kind") == "L")) {
        stop("Returned V is not an lower-tri matrix\n")
      }
    }
    # compute S matrix for adjusting V^{\top} V ----------------------------------
    if (VAdj) {
      v1 <- as.vector(Matrix::t(as.vector(V %*% dpsi_dx)) %*% V) # V^{\top} V dpsi_dx
      v2 <- -H11_dpsi_dx +
        H12_mul(veccCondMeanVarObj, dPsi, D, H21_dpsi_dx / (1 + dPsi))
      S <- v2 / v1
      S_mask <- abs(S) < 1e-8 | abs(S) > 1e8 | is.na(S)
      if (any(S_mask)) {
        if (verbose) {
          cat("Abnormal values in S found, setting them to 1\n")
        }
        S[S_mask] <- 1
      }
      if (retJac) {
        H_hat <- rbind(
          cbind(
            -Matrix::t(V) %*% V %*% diag(S) +
              dpsi_dx_dbeta %*% solve(dpsi_dbeta_dbeta) %*% Matrix::t(dpsi_dx_dbeta),
            dpsi_dx_dbeta
          ),
          cbind(Matrix::t(dpsi_dx_dbeta), dpsi_dbeta_dbeta)
        )
      }
    } else {
      S <- rep(1, n)
    }
    # compute Jac^{-1} \cdot grad -------------------------------------------
    HInv11_dpsi_dx <- HInv11_mul(
      veccCondMeanVarObj, dPsi, D, V, S,
      dpsi_dx
    )
    HInv12_dpsi_dbeta <- HInv12_mul(
      veccCondMeanVarObj, dPsi, D, V, S,
      dpsi_dbeta
    )
    HInv21_dpsi_dx <- HInv21_mul(
      veccCondMeanVarObj, dPsi, D, V, S,
      dpsi_dx
    )
    HInv22_dpsi_dbeta <- HInv22_mul(
      veccCondMeanVarObj, dPsi, D, V, S,
      dpsi_dbeta
    )
  }
  # build return list ----------------------------------------------
  rslt <- list(grad = c(dpsi_dx, dpsi_dbeta))
  if (retProd) {
    rslt[["jac_grad"]] <- c(
      H11_dpsi_dx + H12_dpsi_dbeta,
      H21_dpsi_dx + H22_dpsi_dbeta
    )
  }
  if (retSolv) {
    rslt[["jac_inv_grad"]] <- c(
      HInv11_dpsi_dx + HInv12_dpsi_dbeta,
      HInv21_dpsi_dx + HInv22_dpsi_dbeta
    )
  }
  if (retJac) {
    rslt[["jac"]] <- H
    if (VAdj & retSolv) {
      rslt[["jac_hat"]] <- H_hat
    }
  }
  return(rslt)
}

# Wrapper around the C function psi so that the first argument is `x`
# The inputs `x`, `beta`, `a`, `b`, and `mu` should be of the same length
# This function is used for finding optimal (argmax) `x` given beta
# This function is intended to be used inside `mvnrnd_wrap`
psi_wrapper <- function(x, beta, veccCondMeanVarObj, NN, a, b,
                        mu) {
  psi(
    a, b, NN, mu, veccCondMeanVarObj$cond_mean_coeff,
    sqrt(veccCondMeanVarObj$cond_var), beta, x
  )
}

# Computed dpsi/dx
# The inputs `x`, `beta`, `a`, `b`, and `mu` should be of the same length
# This function is used for finding optimal (argmax) `x` given beta
# This function is intended to be used inside `mvnrnd_wrap`
dpsi_dx <- function(x, beta, veccCondMeanVarObj, NN, a, b, mu) {
  n <- length(a)
  D <- sqrt(veccCondMeanVarObj$cond_var)
  mu_c <- as.vector(veccCondMeanVarObj$A %*% (x - mu)) + mu
  a_tilde_shift <- (a - mu_c) / D - beta
  b_tilde_shift <- (b - mu_c) / D - beta
  log_diff_cdf <- TruncatedNormal::lnNpr(a_tilde_shift, b_tilde_shift)
  pl <- exp(-0.5 * a_tilde_shift^2 - log_diff_cdf) / sqrt(2 * pi)
  pu <- exp(-0.5 * b_tilde_shift^2 - log_diff_cdf) / sqrt(2 * pi)
  Psi <- pl - pu
  # compute grad ------------------------------------------------
  dpsi_dx <- as.vector(Matrix::t(veccCondMeanVarObj$A) %*%
    (beta / D + Psi / D)) - beta / D
  return(dpsi_dx)
}

# # TEST-------------------------------
#
#
# ## example MVN probabilities --------------------------------
# library(GpGp)
# library(Matrix)
# library(VeccTMVN)
# set.seed(123)
# n1 <- 10
# n2 <- 10
# n <- n1 * n2
# locs <- as.matrix(expand.grid((1:n1) / n1, (1:n2) / n2))
# covparms <- c(2, 0.3, 0)
# cov_mat <- matern15_isotropic(covparms, locs)
# a_list <- list(rep(-Inf, n), rep(-1, n), -runif(n) * 2)
# b_list <- list(rep(-2, n), rep(1, n), runif(n) * 2)
# mu <- runif(n)
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
# mu_ord <- mu[ord]
# NNarray <- find_ordered_nn(locs_ord, m = m)
#
# ## Vecchia approx --------------------------------
# U <- get_sp_inv_chol(cov_mat_ord, NNarray)
# cov_mat_Vecc <- solve(U %*% Matrix::t(U))
# vecc_cond_mean_var_obj <- vecc_cond_mean_var_sp(NNarray, cov_mat_ord)
#
# ## Compare dpsi ------------------------------
# for (i in 1:length(a_list_ord)) {
#   a_ord <- a_list_ord[[i]]
#   b_ord <- b_list_ord[[i]]
#   x0_padded <- runif(2 * n)
#   x0_padded[2 * n] <- 0
#   x0 <- x0_padded[-c(n, 2 * n)]
#   a_ord_demean <- a_ord - mu_ord
#   b_ord_demean <- b_ord - mu_ord
#   x0_padded_demean <- x0_padded - c(mu_ord, rep(0, n))
#   x0_demean <- x0_padded_demean[-c(n, 2 * n)]
#   grad_ds <- grad_idea5(x0_demean, vecc_cond_mean_var_obj, a_ord_demean, b_ord_demean)
#   jac_ds <- jac_idea5(x0_demean, vecc_cond_mean_var_obj, a_ord_demean, b_ord_demean)
#   solve_obj <- grad_jacprod_jacsolv_idea5(x0_padded, vecc_cond_mean_var_obj,
#                                           a_ord, b_ord,
#                                           retJac = T, mu = mu_ord
#   )
#   cat("grad error: ", sum(abs(grad_ds - solve_obj$grad[-c(n, 2 * n)])), "\n")
#   cat(
#     "jac error: ", sum(abs(jac_ds - solve_obj$jac[-c(n, 2 * n), -c(n, 2 * n)])),
#     "\n"
#   )
#   H_inv <- solve(solve_obj$jac)
#   jac_grad_test <- as.vector(solve_obj$jac %*% solve_obj$grad)
#   jac_grad_test_2 <- as.vector(Matrix::t(solve_obj$jac_hat) %*% solve_obj$grad)
#   jac_inv_grad_test <- as.vector(H_inv %*% solve_obj$grad)
#   jac_inv_grad_test_2 <- as.vector(solve(solve_obj$jac_hat) %*% solve_obj$grad)
#   cat("jac_grad error: ", sum(abs(jac_grad_test - solve_obj$jac_grad)), "\n")
#   cat("jac_grad error 2: ", sum(abs(jac_grad_test_2 - solve_obj$jac_grad)), "\n")
#   cat("jac_inv_grad error: ", sum(abs(jac_inv_grad_test - solve_obj$jac_inv_grad)), "\n")
#   cat("jac_inv_grad error 2: ", sum(abs(jac_inv_grad_test_2 - solve_obj$jac_inv_grad)), "\n")
# }
