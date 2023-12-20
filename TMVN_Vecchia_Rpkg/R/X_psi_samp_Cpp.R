# Sample from the proposal density in Idea V and compute psi for each sample.
#   Notice that `E[exp(psi)]` is the MVN probability. Zero mean is assumed.
# The `truncnorm` package uses accept-reject sampling and seems to be able to
#   sample from tail truncation although I haven't verified its accuracy in
#   tail sampling.
# Input:
#   veccCondMeanVarObj - contains information of the conditional mean
#     coefficient, the conditional variance, and the NN array of the Vecchia
#     approximation
#   a - lower bound vector for TMVN
#   b - upper bound vector for TMVN
#   beta - parameter of the proposal density
#   N_level1 - First level Monte Carlo sample size
#   N_level2 - Second level Monte Carlo sample size
#   mu - Mean of the MVN prob
# Return a list of two:
#   [[1]] a vector of length N_level1, representing the exp(psi) values
#   [[2]] a vector of length N_level1, representing the exponents
#
sample_psi_idea5_cpp <- function(veccCondMeanVarObj, a, b,
                                 beta = rep(0, length(a)), N_level1 = 10,
                                 N_level2 = 1000, mu = rep(0, length(a))) {
  cond_sd <- sqrt(veccCondMeanVarObj$cond_var)
  exp_psi <- mvndns(a, b, veccCondMeanVarObj$nn, mu,
    veccCondMeanVarObj$cond_mean_coeff,
    cond_sd, beta,
    NLevel1 = N_level1, NLevel2 = N_level2
  )
  return(exp_psi)
}


# # TEST -------------------------------------------------------
# library(GpGp)
# library(TruncatedNormal)
# library(mvtnorm)
# library(nleqslv)
# library(VeccTMVN)
# library(tlrmvnmvt)
#
# ## example MVN probabilities --------------------------------
# n1 <- 10
# n2 <- 10
# n <- n1*n2
# locs <- as.matrix(expand.grid((1 : n1) / n1, (1 : n2) / n2))
# covparms <- c(2, 0.3, 0)
# cov_mat <- matern15_isotropic(covparms, locs)
# a_list <- list(rep(-Inf, n), rep(-1, n), -runif(n) * 2)
# b_list <- list(rep(-2, n), rep(1, n), runif(n) * 2)
#
# ## Compute MVN probs --------------------------------
# N_level1 <- 12  # Level 1 MC size
# N_level2 <- 1e4 # Level 2 MC size
# m <- 30  # num of nearest neighbors
# for(i in 1 : length(a_list)){
#   ### ordering based on integration limits --------------------------------
#   pnorm_at_a <- pnorm(a_list[[i]], sd = sqrt(covparms[1]))
#   pnorm_at_b <- pnorm(b_list[[i]], sd = sqrt(covparms[1]))
#   ord <- order(pnorm_at_b - pnorm_at_a, decreasing = F)
#   locs_ord <- locs[ord, , drop = F]
#   cov_mat_ord <- matern15_isotropic(covparms, locs_ord)
#   a_ord <- a_list[[i]][ord]
#   b_ord <- b_list[[i]][ord]
#   ### NN and Vecchia approx --------------------------------
#   NNarray <- find_ordered_nn(locs_ord, m = m)
#   U <- get_sp_inv_chol(cov_mat_ord, NNarray)
#   cov_mat_Vecc <- solve(U %*% t(U))
#   vecc_cond_mean_var_obj <- vecc_cond_mean_var(cov_mat_ord, NNarray)
#   ### Find proposal parameters -------------------------
#   solv_idea_5 <- nleqslv(rep(0, 2 * n - 2),
#                          fn = grad_idea5,
#                          jac = jac_idea5,
#                          veccCondMeanVarObj = vecc_cond_mean_var_obj,
#                          a = a_ord, b = b_ord,
#                          global = "pwldog",
#                          method = "Newton",
#                          control = list(maxit = 500L))
#   cat("nleqslv finish code is", solv_idea_5$termcd, "\n")
#   beta <- rep(0, n)
#   beta[1 : n - 1] <- solv_idea_5$x[n : (2 * n - 2)]
#   ### Compute MVN prob with idea V -----------------------
#   exp_psi <- sample_psi_idea5_cpp(vecc_cond_mean_var_obj, a_ord, b_ord,
#                           beta = beta, N_level1 = N_level1,
#                           N_level2 = N_level2)
#   exp_psi <- exp_psi[[1]] * exp(exp_psi[[2]])
#   est_tilt_quasi <- mean(exp_psi)
#   err_tilt_quasi <- sd(exp_psi) / sqrt(N_level1)
#
#   ### Compute MVN prob with other methods -----------------------
#   est_TN <- TruncatedNormal::pmvnorm(
#     rep(0, n), cov_mat_Vecc, lb = a_ord, ub = b_ord)
#   est_TLR <- tlrmvnmvt::pmvn(a_ord, b_ord, sigma = cov_mat_Vecc)
#   est_Genz <- mvtnorm::pmvnorm(a_ord, b_ord, sigma = cov_mat_Vecc)
#   cat("est_tilt_quasi", est_tilt_quasi, "err_tilt_quasi", err_tilt_quasi, "\n",
#       "est_TN", est_TN, "err_TN", attributes(est_TN)$relerr * est_TN, "\n",
#       "est_TLR", est_TLR, "err_TLR", attributes(est_TLR)$error, "\n",
#       "est_Genz", est_Genz, "err_Genz", attributes(est_Genz)$error, "\n"
#   )
# }
