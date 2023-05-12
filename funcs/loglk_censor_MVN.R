library(GpGp)

#' Compute censored multivariate normal (MVN) log-probabilities that have
#' spatial covariance matrices using Vecchia approximation
#'
#' @param locs location (feature) matrix n X d
#' @param indCensor indices of locations that have only censored observations
#' @param yObs observed (not censored) values, of length n
#' @param bCensor upper bound, above which observations are not censored,
#' can be different for different locations, of length 1 or n
#' @param covName covariance function name from the `GpGp` package
#' @param covParms parameters for `covName`
#' @param m Vecchia conditioning set size
#' @param ... can be the following parameters that will be passed to `pmvn`
#' @param NLevel1 first level Monte Carlo sample size
#' @param NLevel2 second level Monte Carlo sample size
#' @param verbose verbose or not
#' @return estimated MVN probability and estimation error
#'

loglk_censor_MVN <- function(locs, indCensor, yObs, bCensor,
                             covName = NULL, covParms = NULL, m = 30, ...) {
  n <- nrow(locs)
  locs_obs <- locs[-indCensor]
  obcs_censor <- locs[indCensor]
  odr_censor <- GpGp::order_maxmin(obcs_censor)
}
