require(mvPot)
require(CDFNormalAproxPackNoC)
#' Function implementing VCDF provided by \url{https://www.tandfonline.com/doi/full/10.1080/00949655.2021.2016759}
#' @param upper upper integration limits, assuming lower integration limits are negative infinity
#' @param mean mean parameter of the MVN probability
#' @param sigma covariance matrix parameter of the MVN probability
#' @param p controling the Monte Carlo sample size
#' @return the estimate of the logarithm of the MVN probability
CDFNormalAproxPackNoC.pmvn <- function(upper, mean, sigma, p = 499) {
  corr <- cov2cor(sigma)
  ### Using 50 neighbors sequentially and 10 joint estimation
  joint <- 10
  use <- 50
  latticelist50 <- lapply(X = 2:(use * joint + joint), FUN = function(x) {
    genVecQMC(p, (x))
  })
  pmvNorm3(
    M = mean, Sigma = sigma, Corr = corr, upper = upper, use = use,
    lattice = latticelist50, joint = joint
  )
}
