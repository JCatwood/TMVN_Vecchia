require(mvPot)
require(CDFNormalAproxPackNoC)
CDFNormalAproxPackNoC.pmvn <- function(upper, mean, sigma){
  corr <- cov2cor(sigma)
  ###Using 50 neighbors sequentially and 10 joint estimation
  joint <- 10
  use <- 50
  p <- 499
  latticelist50 <- lapply(X = 2:(use*joint+joint),FUN = function(x){
    genVecQMC(p, (x))})
  pmvNorm3(M = mean,Sigma = sigma,Corr=corr,upper = upper,use = use,
           lattice = latticelist50,joint = joint)
}
