library(GpGp)
library(VeccTMVN)
library(nleqslv)

## local pmvn func for testing --------------------------------
pmvn <- function(a, b, vecc_cond_mean_var_obj) {
  solv_idea_5 <- nleqslv(rep(0, 2 * n - 2),
    fn = grad_idea5,
    jac = jac_idea5,
    veccCondMeanVarObj = vecc_cond_mean_var_obj,
    a = a, b = b,
    global = "pwldog",
    method = "Newton",
    control = list(maxit = 500L)
  )
  cat("nleqslv finish code is", solv_idea_5$termcd, "\n")
  beta <- rep(0, n)
  beta[1:n - 1] <- solv_idea_5$x[n:(2 * n - 2)]

  exp_psi <- sample_psi_idea5_cpp(vecc_cond_mean_var_obj, a, b,
    beta = beta, N_level1 = 12,
    N_level2 = 1e4
  )
  est_tilt_quasi <- mean(exp_psi)
  err_tilt_quasi <- sd(exp_psi) / sqrt(N_level1)
  cat("est_tilt_quasi", est_tilt_quasi, "err_tilt_quasi", err_tilt_quasi, "\n")
}

## example MVN probabilities --------------------------------
n1 <- 30
n2 <- 30
n <- n1 * n2
locs <- as.matrix(expand.grid((1:n1) / n1, (1:n2) / n2))
covparms <- c(2, 0.1, 0)
cov_mat <- matern15_isotropic(covparms, locs)
a_list <- list(rep(-Inf, n), rep(-Inf, n), rep(-Inf, n))
b_list <- list(rep(-2, n), rep(1, n), -runif(n) * 2)

## Compute MVN probs --------------------------------
N_level1 <- 12 # Level 1 MC size
N_level2 <- 1e4 # Level 2 MC size
m <- 30 # num of nearest neighbors
for (i in 1:length(a_list)) {
  ### ordering based on integration limits --------------------------------
  pnorm_at_a <- pnorm(a_list[[i]], sd = sqrt(covparms[1]))
  pnorm_at_b <- pnorm(b_list[[i]], sd = sqrt(covparms[1]))
  ord <- order(pnorm_at_b - pnorm_at_a, decreasing = F)
  locs_ord <- locs[ord, , drop = F]
  cov_mat_ord <- matern15_isotropic(covparms, locs_ord)
  a_ord <- a_list[[i]][ord]
  b_ord <- b_list[[i]][ord]
  ### NN and Vecchia approx --------------------------------
  NNarray <- find_ordered_nn(locs_ord, m = m)
  vecc_cond_mean_var_obj <- vecc_cond_mean_var(cov_mat_ord, NNarray)
  ### scale s.t. diag(L) is 1 -----------------------------
  diag_L_Vecc <- sqrt(vecc_cond_mean_var_obj$cond_var)
  a_ord_scale <- a_ord / diag_L_Vecc
  b_ord_scale <- b_ord / diag_L_Vecc
  cov_mat_ord_scale <- t(t(cov_mat_ord / diag_L_Vecc) / diag_L_Vecc)
  vecc_cond_mean_var_obj_scale <- vecc_cond_mean_var(cov_mat_ord_scale, NNarray)
  ### compare ---------------------
  pmvn(a_ord, b_ord, vecc_cond_mean_var_obj)
  pmvn(a_ord_scale, b_ord_scale, vecc_cond_mean_var_obj_scale)
}
