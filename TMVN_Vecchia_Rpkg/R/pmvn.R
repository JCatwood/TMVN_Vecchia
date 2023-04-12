library(truncnorm)


#' Compute multivariate normal (MVN) probabilities that have spatial covariance
#' matrices using Vecchia approximation
#'
#' @param lower lower bound vector for TMVN
#' @param upper upper bound vector for TMVN
#' @param mean MVN mean
#' @param locs location (feature) matrix n X d
#' @param covName covariance function name from the `GpGp` package
#' @param covParms parameters for `covName`
#' @param m Vecchia conditioning set size
#' @param sigma dense covariance matrix, not needed when `locs` is not null
#' @param NLevel1 first level Monte Carlo sample size
#' @param NLevel2 second level Monte Carlo sample size
#' @param verbose verbose or not
#' @return estimated MVN probability and estimation error
#'
pmvn <- function(lower, upper, mean, locs, covName = "matern15_isotropic",
                 covParms = c(1.0, 0.1, 0.0), m = 30, sigma = NULL,
                 NLevel1 = 12, NLevel2 = 1e4, verbose = F) {
  # standardize the input MVN prob -----------------------------
  lower <- lower - mean
  upper <- upper - mean
  if (is.null(sigma)) {
    n <- nrow(locs)
    use_sigma <- F
    margin_sd <- sqrt(covParms[1])
    upper <- upper / margin_sd
    lower <- lower / margin_sd
    covParms[1] <- 1
  } else {
    n <- nrow(sigma)
    use_sigma <- T
    margin_sd <- sqrt(diag(sigma))
    upper <- upper / margin_sd
    lower <- lower / margin_sd
    sigma <- t(t(sigma / margin_sd) / margin_sd)
  }
  tmvn_prob_1D <- pnorm(upper) - pnorm(lower)
  if (any(tmvn_prob_1D < 0)) {
    stop("Invalid MVN probability. Truncated marginal
         probabilities have negative value(s)\n")
  }
  lower_upper <- matrix(0, n, 2)
  lower_upper[, 1] <- lower
  lower_upper[, 2] <- upper
  # reorder based on tmvn_prob_1D --------------------------------
  ord <- order(tmvn_prob_1D, decreasing = F)
  lower_ord <- lower_upper[ord, 1]
  upper_ord <- lower_upper[ord, 2]
  if (use_sigma) {
    sigma_ord <- sigma[ord, ord, drop = F]
  } else {
    locs_ord <- locs[ord, , drop = F]
  }
  # find nearest neighbors for Vecchia --------------------------------
  if (use_sigma) {
    NN <- find_nn_corr(sigma_ord, m)
  } else {
    NN <- GpGp::find_ordered_nn(locs_ord, m)
  }
  # find Vecchia approx object -----------------------------------
  if (use_sigma) {
    vecc_obj <- vecc_cond_mean_var_sp(NN, covMat = sigma_ord)
  } else {
    vecc_obj <- vecc_cond_mean_var_sp(NN,
      locs = locs_ord, covName = covName,
      covParms = covParms
    )
  }
  # find tilting parameter beta -----------------------------------
  trunc_expect <- etruncnorm(lower_ord, upper_ord)
  x0 <- c(trunc_expect, rep(0, n))
  solv_idea_5_sp <- optim(
    x0,
    fn = function(x, ...) {
      ret <- grad_jacprod_jacsolv_idea5(x, ...,
        retJac = F,
        retProd = F, retSolv = F
      )
      0.5 * sum((ret$grad)^2)
    },
    gr = function(x, ...) {
      ret <- grad_jacprod_jacsolv_idea5(x, ...,
        retJac = F,
        retProd = T, retSolv = F
      )
      ret$jac_grad
    },
    method = "L-BFGS-B",
    veccCondMeanVarObj = vecc_obj,
    a = lower_ord, b = upper_ord, verbose = verbose,
    lower = c(lower_ord, rep(-Inf, n)), upper = c(upper_ord, rep(Inf, n))
  )
  if (verbose) {
    cat(
      "Gradient norm at the optimal beta is", sqrt(2 * solv_idea_5_sp$value),
      "\n"
    )
  }
  if (any(solv_idea_5_sp$par[1:n] < lower_ord) ||
    any(solv_idea_5_sp$par[1:n] > upper_ord)) {
    warning("Optimal x is outside the integration region during minmax tilting\n")
  }
  beta <- solv_idea_5_sp$par[(n + 1):(2 * n)]
  # compute MVN probs and est error ---------------------------------
  exp_psi <- sample_psi_idea5_cpp(vecc_obj, lower_ord, upper_ord,
    beta = beta, N_level1 = N_level1,
    N_level2 = N_level2
  )
  est_prob <- mean(exp_psi)
  est_prob_err <- sd(exp_psi) / sqrt(N_level1)
  attr(est_prob, "error") <- est_prob_err
  return(est_prob)
}


# TEST -------------------------------------------------------

# library(GpGp)
# library(TruncatedNormal)
# library(mvtnorm)
# library(VeccTMVN)
#
#
# ## example MVN probabilities --------------------------------
# set.seed(123)
# n1 <- 10
# n2 <- 10
# n <- n1 * n2
# locs <- as.matrix(expand.grid((1:n1) / n1, (1:n2) / n2))
# covparms <- c(2, 0.3, 0)
# cov_mat <- matern15_isotropic(covparms, locs)
# a_list <- list(rep(-Inf, n), rep(-1, n), -runif(n) * 2 - 4)
# b_list <- list(rep(-2, n), rep(1, n), -runif(n) * 2)
#
# ## Compute MVN probs --------------------------------
# N_level1 <- 12 # Level 1 MC size
# N_level2 <- 1e4 # Level 2 MC size
# m <- 30 # num of nearest neighbors
# for (i in 1:length(a_list)) {
#   a <- a_list[[i]]
#   b <- b_list[[i]]
#   ### Compute MVN prob with idea V -----------------------
#   time_Vecc <- system.time(est_Vecc <- VeccTMVN::pmvn(a, b, 0, locs,
#                                                       covName = "matern15_isotropic",
#                                                       covParms = covparms, m = m, verbose = F
#   ))[[3]]
#   ### Compute MVN prob with other methods -----------------------
#   time_TN <- system.time(est_TN <- TruncatedNormal::pmvnorm(
#     rep(0, n), cov_mat,
#     lb = a, ub = b
#   ))[[3]]
#   time_TLR <- system.time(
#     est_TLR <- tlrmvnmvt::pmvn(a, b, sigma = cov_mat)
#   )[[3]]
#   time_Genz <- system.time(
#     est_Genz <- mvtnorm::pmvnorm(a, b, sigma = cov_mat)
#   )[[3]]
#   cat(
#     "est_Vecc", est_Vecc, "err_Vecc", attributes(est_Vecc)$error, "time_Vecc",
#     time_Vecc, "\n"
#   )
#   cat(
#     "est_TN", est_TN, "err_TN", attributes(est_TN)$relerr * est_TN, "time_TN",
#     time_TN, "\n"
#   )
#   cat(
#     "est_TLR", est_TLR, "err_TLR", attributes(est_TLR)$error, "time_TLR",
#     time_TLR, "\n"
#   )
#   cat(
#     "est_Genz", est_Genz, "err_Genz", attributes(est_Genz)$error, "time_Genz",
#     time_Genz, "\n"
#   )
# }
