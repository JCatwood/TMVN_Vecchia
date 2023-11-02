library(GpGp)
library(truncnorm)

#' Simulate censored multivariate normal (MVN) as censored locations using the
#'   Vecchia approximation
#'
#' @param lower lower bound vector for TMVN
#' @param upper upper bound vector for TMVN
#' @param mean MVN mean
#' @param locs location (feature) matrix n X d
#' @param covName covariance function name from the `GpGp` package
#' @param covParms parameters for `covName`
#' @param m Vecchia conditioning set size
#' @param sigma dense covariance matrix, not needed when `locs` is not null
#' @param N number of samples required
#' @param verbose verbose level
#' @return n X N matrix of generated samples
#'
#' @export
mvrandn <- function(lower, upper, mean, locs = NULL,
                    covName = "matern15_isotropic", covParms = c(1.0, 0.1, 0.0),
                    m = 30, sigma = NULL, N = 1e3, verbose = F) {
  ## standardize the input MVN prob -----------------------------
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
  if (any(lower < -10)) {
    lower[lower < -10] <- -10
  }
  if (any(upper < -10)) {
    upper[upper < -10] <- -10
  }
  if (any(lower > 10)) {
    lower[lower > 10] <- 10
  }
  if (any(upper > 10)) {
    upper[upper > 10] <- 10
  }
  if (any(upper < lower)) {
    stop("Invalid MVN probability. Truncated marginal
         probabilities have negative value(s)\n")
  }
  lower_upper <- matrix(0, n, 2)
  lower_upper[, 1] <- lower
  lower_upper[, 2] <- upper
  lower <- lower_upper[, 1]
  upper <- lower_upper[, 2]
  ## reorder --------------------------------
  if (use_sigma) {
    ord <- Vecc_reorder(lower, upper, m, covMat = sigma)$order
  } else {
    ord <- Vecc_reorder(
      lower, upper, m, locs, covName, covParms
    )$order
  }
  lower <- lower[ord]
  upper <- upper[ord]
  if (use_sigma) {
    sigma <- sigma[ord, ord, drop = F]
  } else {
    locs <- locs[ord, , drop = F]
  }
  ## find nearest neighbors for Vecchia --------------------------------
  if (use_sigma) {
    NN <- find_nn_corr(sigma, m)
  } else {
    NN <- GpGp::find_ordered_nn(locs, m)
  }
  ## find Vecchia approx object -----------------------------------
  if (use_sigma) {
    vecc_obj <- vecc_cond_mean_var_sp(NN, covMat = sigma)
    idx <- which(vecc_obj$cond_var < 0.01)
    if (length(idx) > 0) {
      diag(sigma)[idx] <- diag(sigma)[idx] + 0.01
      vecc_obj <- vecc_cond_mean_var_sp(NN, covMat = sigma)
    }
  } else {
    vecc_obj <- vecc_cond_mean_var_sp(NN,
      locs = locs, covName = covName,
      covParms = covParms
    )
    if (any(vecc_obj$cond_var < 0.01)) {
      covParms[length(covParms)] <- 0.01 # nugget
      vecc_obj <- vecc_cond_mean_var_sp(NN,
        locs = locs, covName = covName,
        covParms = covParms
      )
    }
  }
  ## compute MVN probs and est error ---------------------------------
  samp_Vecc_ord <- mvnrnd_wrap(
    lower, upper,
    mu = 0,
    NN = NN, veccObj = vecc_obj, N = N, verbose = verbose
  )
  ord_rev <- integer(n)
  ord_rev[ord] <- 1:n
  samp_Vecc <- margin_sd * samp_Vecc_ord[ord_rev, ] + mean
  return(samp_Vecc)
}


# TEST -------------------------------------------------------
# library(GpGp)
# library(mvtnorm)
# library(TruncatedNormal)
# library(VeccTMVN)
# set.seed(123)
# n1 <- 10
# n2 <- 10
# n <- n1 * n2
# locs <- as.matrix(expand.grid((1:n1) / n1, (1:n2) / n2))
# covparms <- c(2, 0.1, 0)
# mu <- rep(1, n)
# N <- 1000
# cov_mat <- matern15_isotropic(covparms, locs)
# a <- rep(-1, n)
# b <- rep(-0, n)
# samp_TN <- TruncatedNormal::mvrandn(
#   a, b, cov_mat,
#   n = N, mu = mu
# )
# samp_Vecc <- mvrandn(
#   a, b, mu, locs, "matern15_isotropic", covparms,
#   m = 30, N = N, verbose = T
# )
# ##  histogram for verification -------------------
# par(mfrow = c(1, 2))
# hist(samp_Vecc, main = "Vecc Samples")
# hist(samp_TN, main = "TN Samples")
# image(matrix(samp_TN[, 1], n1, n2))
# image(matrix(samp_Vecc[, 1], n1, n2))
