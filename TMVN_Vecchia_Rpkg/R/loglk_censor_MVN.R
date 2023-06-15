library(GpGp)
library(truncnorm)

#' Compute censored multivariate normal (MVN) log-probabilities that have
#' spatial covariance matrices using Vecchia approximation
#'
#' @param locs location (feature) matrix n X d
#' @param indCensor indices of locations that have only censored observations
#' @param y observed (not censored) values, of length n
#' @param bCensor upper bound, above which observations are not censored,
#' can be different for different locations, of length 1 or n
#' @param covName covariance function name from the `GpGp` package
#' @param covParms parameters for `covName`
#' @param m Vecchia conditioning set size
#' @param NLevel1 first level Monte Carlo sample size
#' @param NLevel2 second level Monte Carlo sample size
#' @param verbose verbose level
#' @return estimated MVN probability and estimation error
#'

loglk_censor_MVN <- function(locs, indCensor, y, bCensor,
                             covName = NULL, covParms = NULL, m = 30,
                             NLevel1 = 10, NLevel2 = 1e3,
                             verbose = T) {
  # extract and separate observed and censored data --------------------------------
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
  # maxmin order --------------------------------
  # odr_obs <- GpGp::order_maxmin(locs_obs)
  # odr_censor <- GpGp::order_maxmin(locs_censor)
  # locs_obs <- locs_obs[odr_obs, , drop = F]
  # locs_censor <- locs_censor[odr_censor, , drop = F]
  # y_obs <- y_obs[odr_obs]
  # b_censor <- b_censor[odr_censor]
  locs <- rbind(locs_obs, locs_censor)
  # NN and Vecchia approx obj --------------------------------
  NN <- GpGp::find_ordered_nn(locs, m)
  vecc_obj <- vecc_cond_mean_var_sp(NN,
    locs = locs, covName = covName,
    covParms = covParms
  )
  # log pdf for observed data --------------------------------
  loglik_pdf <- GpGp::vecchia_meanzero_loglik(
    covParms, covName, y_obs,
    locs_obs, NN[1:n_obs, , drop = F]
  )$loglik
  # vecc approx object for censored data --------------------------------
  y_obs_padded <- c(y_obs, rep(0, n_censor))
  cond_mean <- as.vector(vecc_obj$A %*% y_obs_padded)[(n_obs + 1):n]
  a <- rep(-Inf, n_censor)
  b <- b_censor
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
  # find tilting parameter beta -----------------------------------
  trunc_expect <- etruncnorm(a, b / sqrt(covParms[1]))
  x0 <- c(trunc_expect, rep(0, n_censor))
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
    veccCondMeanVarObj = vecc_obj_censor,
    a = a, b = b, mu = cond_mean, verbose = verbose,
    lower = c(a, rep(-Inf, n)), upper = c(b, rep(Inf, n)),
    control = list(maxit = 500)
  )
  if (verbose) {
    cat(
      "Gradient norm at the optimal beta is", sqrt(2 * solv_idea_5_sp$value),
      "\n"
    )
  }
  beta <- solv_idea_5_sp$par[(n_censor + 1):(2 * n_censor)]
  # compute MVN probs and est error ---------------------------------
  exp_psi <- sample_psi_idea5_cpp(vecc_obj_censor, a, b,
    beta = beta, NLevel1, NLevel2, mu = cond_mean
  )
  est_prob <- mean(exp_psi)
  # est_prob_err <- sd(exp_psi) / sqrt(NLevel1)
  # attr(est_prob, "error") <- est_prob_err
  return(log(est_prob) + loglik_pdf)
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
# cov_mat <- matern15_isotropic(covparms, locs)
# y <- as.vector(t(chol(cov_mat)) %*% rnorm(n))
# odr <- order(y, decreasing = T)
# y <- y[odr]
# locs <- locs[odr, , drop = F]
# cov_mat <- cov_mat[odr, odr]
# b_censor <- 1
# ind_censor <- which(y < b_censor)
# ind_obs <- which(!(y < b_censor))
# logpdf <- mvtnorm::dmvnorm(y[ind_obs],
#   sigma = cov_mat[ind_obs, ind_obs],
#   log = T
# )
# tmp_mat <- cov_mat[ind_censor, ind_obs] %*% solve(cov_mat[ind_obs, ind_obs])
# cond_mean <- as.vector(tmp_mat %*% y[ind_obs])
# cond_cov_mat <- cov_mat[ind_censor, ind_censor] -
#   tmp_mat %*% cov_mat[ind_obs, ind_censor]
# logcdf <- log(TruncatedNormalBeta::pmvnorm(cond_mean, cond_cov_mat, ub = b_censor))
# loglk_Vecc <- loglk_censor_MVN(locs, ind_censor, y, b_censor, "matern15_isotropic",
#   covparms,
#   m = 50
# )
# 
# ranges <- seq(from = 0.07, to = 0.1, by = 0.001)
# loglk_Vecc_vec <- rep(0, length(ranges))
# idx <- 1
# for (myrange in ranges) {
#   covparms[2] <- myrange
#   set.seed(123)
#   loglk_Vecc_vec[idx] <- loglk_censor_MVN(locs, ind_censor, y, b_censor, "matern15_isotropic",
#     covparms,
#     m = 50
#   )
#   idx <- idx + 1
# }
# 
# plot(ranges, loglk_Vecc_vec)
