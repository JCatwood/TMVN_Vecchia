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
#' @export
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
  # reorder --------------------------------
  locs <- rbind(locs_obs, locs_censor)
  # odr <- FIC_reorder_univar(c(y_obs, rep(-Inf, n_censor)), c(y_obs, b_censor),
  #                           m, locs = locs, covName = covName,
  #                           covParms = covParms)
  # if(any(odr[(n_obs + 1) : n] <= n_obs)){
  #   warning("FIC_reorder_univar failed\n")
  # }else{
  #   odr_censor <- odr[(n_obs + 1) : n] - n_obs
  #   locs_censor <- locs_censor[odr_censor, , drop = F]
  #   b_censor <- b_censor[odr_censor]
  #   locs <- rbind(locs_obs, locs_censor)
  # }
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
  trunc_expect <- truncnorm::etruncnorm(a, b, mean = cond_mean, sd = sqrt(covParms[1]))
  x0 <- c(trunc_expect, rep(0, n_censor))
  solv_idea_5_sp <- stats::optim(
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
  exponent <- min(exp_psi[[2]])
  log_est_prob <- exponent +
    log(mean(exp_psi[[1]] * exp(exp_psi[[2]] - exponent)))
  return(log_est_prob + loglik_pdf)
}
