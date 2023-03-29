library(truncnorm)
library(TruncatedNormal)
library(randtoolbox)


#' Sample from the proposal density in Idea V and compute psi for each sample.
#'   Notice that `E[exp(psi)]` is the MVN probability. Zero mean is assumed.
#' The `truncnorm` package uses accept-reject sampling and seems to be able to
#'   sample from tail truncation although I haven't verified its accuracy in
#'   tail sampling.
#' Input:
#'   N - number of samples to draw
#'   veccCondMeanVarObj - contains information of the conditional mean
#'     coefficient, the conditional variance, and the NN array of the Vecchia
#'     approximation
#'   a - lower bound vector for TMVN
#'   b - upper bound vector for TMVN
#'   beta - parameter of the proposal density
#' Return the a vector of length N, representing the psi values
#'
sample_psi_idea5 <- function(N, veccCondMeanVarObj, a, b,
                             beta = rep(0, length(x)), usePseudo = F){
  n <- length(a)
  X <- matrix(NA, n, N)
  lnNpr_sum <- rep(0, N)
  inner_prod <- rep(0, N)
  if(usePseudo)
    x <- t(as.matrix(randtoolbox::sobol(
      N, dim = n, init =TRUE, scrambling = 1, seed=ceiling(1e6*runif(1)))))
  for(i in 1 : n){
    ind <- veccCondMeanVarObj$nn[i, -1]
    sd <- sqrt(veccCondMeanVarObj$cond_var[i]) # scalar
    mu <- apply(veccCondMeanVarObj$cond_mean_coeff[i, ] * X[ind, ], 2, sum,
                na.rm = T) # vector of length N
    a_tilde_i_mbeta <- (a[i] - mu) / sd - beta[i]
    b_tilde_i_mbeta <- (b[i] - mu) / sd - beta[i]
    if(usePseudo){
      X[i, ] <- norminvp(x[i, ], a_tilde_i_mbeta, b_tilde_i_mbeta) * sd + mu +
        beta[i] * sd
    }else{
      X[i, ] <- rtruncnorm(n = 1, a = a[i], b = b[i], mean = mu + beta[i] * sd,
                           sd = sd)
    }

    lnNpr_sum <- lnNpr_sum +
      TruncatedNormal::lnNpr(a_tilde_i_mbeta, b_tilde_i_mbeta)
    inner_prod <- inner_prod + (X[i, ] - mu) * beta[i] / sd

  }
  psi <- 0.5 * sum(beta^2) - inner_prod + lnNpr_sum
  return(psi)
}


# # TEST -------------------------------------------------------
# library(GpGp)
# library(TruncatedNormal)
# library(mvtnorm)
# library(nleqslv)
# source("inv_chol.R")
# source("vecc_cond_mean_var.R")
# source("gradpsi.R")
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
# N <- 1e4  # MC sample size
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
#   psi <- sample_psi_idea5(N, vecc_cond_mean_var_obj, a_ord, b_ord,
#                           beta = beta, usePseudo = T)
#   est_tilt_quasi <- mean(exp(psi))
#   err_tilt_quasi <- sd(apply(matrix(exp(psi), nrow = 10), 1, mean)) / sqrt(10)
#   psi <- sample_psi_idea5(N, vecc_cond_mean_var_obj, a_ord, b_ord,
#                           beta = beta, usePseudo = F)
#   est_tilt <- mean(exp(psi))
#   err_tilt <- sd(apply(matrix(exp(psi), nrow = 10), 1, mean)) / sqrt(10)
#   psi <- sample_psi_idea5(N, vecc_cond_mean_var_obj, a_ord, b_ord,
#                           beta = rep(0, n), usePseudo = T)
#   est_quasi <- mean(exp(psi))
#   err_quasi <- sd(apply(matrix(exp(psi), nrow = 10), 1, mean)) / sqrt(10)
#   psi <- sample_psi_idea5(N, vecc_cond_mean_var_obj, a_ord, b_ord,
#                           beta = rep(0, n), usePseudo = F)
#   est <- mean(exp(psi))
#   err <- sd(apply(matrix(exp(psi), nrow = 10), 1, mean)) / sqrt(10)
#
#   ### Compute MVN prob with other methods -----------------------
#   est_TN <- TruncatedNormal::pmvnorm(
#     rep(0, n), cov_mat_Vecc, lb = a_ord, ub = b_ord)
#   est_TLR <- tlrmvnmvt::pmvn(a_ord, b_ord, sigma = cov_mat_Vecc)
#   est_Genz <- mvtnorm::pmvnorm(a_ord, b_ord, sigma = cov_mat_Vecc)
#   cat("est_tilt_quasi", est_tilt_quasi, "err_tilt_quasi", err_tilt_quasi, "\n",
#       "est_tilt", est_tilt, "err_tilt", err_tilt, "\n",
#       "est_quasi", est_quasi, "err_quasi", err_quasi, "\n",
#       "est", est, "err", err, "\n",
#       "est_TN", est_TN, "err_TN", attributes(est_TN)$relerr * est_TN, "\n",
#       "est_TLR", est_TLR, "err_TLR", attributes(est_TLR)$error, "\n",
#       "est_Genz", est_Genz, "err_Genz", attributes(est_Genz)$error, "\n"
#       )
# }
