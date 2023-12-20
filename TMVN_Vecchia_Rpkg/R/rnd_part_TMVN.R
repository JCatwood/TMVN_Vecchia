library(GpGp)
library(truncnorm)

#' Simulate censored multivariate normal (MVN) as censored locations using the
#'   Vecchia approximation
#'
#' @param locs location (feature) matrix n X d
#' @param indCensor indices of locations that have only censored observations
#' @param y observed (not censored) values, of length n
#' @param bCensor upper bound, above which observations are not censored,
#' can be different for different locations, of length 1 or n
#' @param covName covariance function name from the `GpGp` package
#' @param covParms parameters for `covName`
#' @param m Vecchia conditioning set size
#' @param N number of samples required
#' @param verbose verbose level
#' @param reorder whether to Vecchia univariate variable reordering
#' @return n X N matrix of generated samples
#'
#' @export
ptmvrandn <- function(locs, indCensor, y, bCensor,
                      covName = NULL, covParms = NULL, m = 30,
                      N = 1e3, verbose = T, reorder = T) {
  ## standardize variance ------------------------------
  parm_sd <- sqrt(covParms[1])
  y <- y / parm_sd
  bCensor <- bCensor / parm_sd
  covParms[1] <- 1
  ## extract and separate observed and censored data ---------------------------
  n <- nrow(locs)
  n_obs <- n - length(indCensor)
  n_censor <- length(indCensor)
  if (n_obs < 1 | n_censor < 1) {
    stop("loglk_censor_MVN should be called with
         non-empty censor/observed data")
  }
  locs_obs <- locs[-indCensor, , drop = F]
  locs_censor <- locs[indCensor, , drop = F]
  y_obs <- y[-indCensor]
  if (length(bCensor) > 1) {
    b_censor <- bCensor[indCensor]
  } else {
    b_censor <- rep(bCensor, n_censor)
  }
  ## reorder --------------------------------
  locs <- rbind(locs_obs, locs_censor)
  if (reorder) {
    ord <- Vecc_reorder(
      c(y_obs, rep(-Inf, n_censor)), c(y_obs, b_censor),
      m,
      locs = locs, covName = covName,
      covParms = covParms
    )$order
  } else {
    ord <- 1:n
  }
  if (any(ord[(n_obs + 1):n] <= n_obs)) {
    warning("Vecc_reorder failed\n")
  } else {
    odr_censor <- ord[(n_obs + 1):n] - n_obs
    locs_censor <- locs_censor[odr_censor, , drop = F]
    b_censor <- b_censor[odr_censor]
    locs <- rbind(locs_obs, locs_censor)
  }
  ## NN and Vecchia approx obj --------------------------------
  NN <- GpGp::find_ordered_nn(locs, m)
  vecc_obj <- vecc_cond_mean_var_sp(NN,
    locs = locs, covName = covName,
    covParms = covParms
  )
  ## vecc approx object for censored data --------------------------------
  y_obs_padded <- c(y_obs, rep(0, n_censor))
  cond_mean <- as.vector(vecc_obj$A %*% y_obs_padded)[(n_obs + 1):n]
  a <- rep(-10, n_censor)
  b <- b_censor
  if (any(a > b)) {
    stop("b_censor is too small (smaller than 10 std) \n")
  }
  cond_var <- vecc_obj$cond_var[(n_obs + 1):n]
  NN <- NN[(n_obs + 1):n, , drop = F] - n_obs
  NN_mask <- NN <= 0
  NN_mask[is.na(NN_mask)] <- F
  cond_mean_coeff <- vecc_obj$cond_mean_coeff[(n_obs + 1):n, , drop = F]
  cond_mean_coeff[NN_mask[, -1, drop = F]] <- 0
  NN[NN_mask] <- 1
  A <- vecc_obj$A[(n_obs + 1):n, (n_obs + 1):n, drop = F]
  vecc_obj_censor <- list(
    cond_mean_coeff = cond_mean_coeff, cond_var = cond_var,
    nn = NN, A = A
  )
  ## compute MVN probs and est error ---------------------------------
  samp_Vecc_ord <- mvnrnd_wrap(
    a, b,
    mu = cond_mean,
    NN = NN, veccObj = vecc_obj_censor, N = N, verbose = 1
  ) * parm_sd
  ord_rev <- integer(n_censor)
  ord_rev[odr_censor] <- 1:n_censor
  samp_Vecc <- matrix(0, n, N)
  samp_Vecc[-indCensor, ] <- y_obs
  samp_Vecc[indCensor, ] <- samp_Vecc_ord[ord_rev, ]
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
# N <- 1000
# cov_mat <- matern15_isotropic(covparms, locs)
# y <- as.vector(t(chol(cov_mat)) %*% rnorm(n))
# b_censor <- 1
# ind_censor <- which(y < b_censor)
# ind_obs <- which(!(y < b_censor))
# tmp_mat <- cov_mat[ind_censor, ind_obs] %*% solve(cov_mat[ind_obs, ind_obs])
# cond_mean <- as.vector(tmp_mat %*% y[ind_obs])
# cond_cov_mat <- cov_mat[ind_censor, ind_censor] -
#   tmp_mat %*% cov_mat[ind_obs, ind_censor]
# samp_TN <- TruncatedNormal::mvrandn(
#   rep(-Inf, length(cond_mean)), rep(b_censor, length(cond_mean)),
#   cond_cov_mat,
#   n = N, mu = cond_mean
# )
# samp_Vecc <- ptmvrandn(
#   locs, ind_censor, y, b_censor, "matern15_isotropic", covparms,
#   m = 30, N = N
# )
# ##  histogram for verification -------------------
# par(mfrow = c(1, 2))
# hist(samp_Vecc[ind_censor, ], main = "Vecc Samples")
# hist(samp_TN, main = "TN Samples")
# samp_TN_padded <- matrix(0, n, N)
# samp_TN_padded[ind_obs, ] <- y[ind_obs]
# samp_TN_padded[ind_censor, ] <- samp_TN
# image(matrix(samp_TN_padded[, 1], n1, n2))
# image(matrix(samp_Vecc[, 1], n1, n2))
