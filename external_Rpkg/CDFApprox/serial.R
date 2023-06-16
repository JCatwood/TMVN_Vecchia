require(mvtnorm)
require(mvPot)
require(doParallel)
require(foreach)
require(CDFNormalAproxPackNoC)


###Defining areas locations
############################################################
A <- expand.grid(1:100,1:100)
A <- A[sample(nrow(A)),]
###Calculating covariance and correlation matrix using#####
###exponential function#
############################################################
D <- as.matrix(dist(A,upper = T,diag = T))
sigma <- 1
rho <- 10
Sigma <- sigma^2*exp(-D/rho)
Corr <- cov2cor(Sigma)
###Creating mean vector
############################################################
M <- rep(0,nrow(A))
###Creating upper limit
############################################################
X <- M
###Creating lattice used by the QuasiMonteCarlo approach
############################################################
###Full distribution
p <- 3607
latticeRule <- genVecQMC(3607, (nrow(A) ))
###Using 50 neighbors sequentially and 10 joint estimation
joint <- 10
use <- 50
p <- 499
latticelist50 <- lapply(X = 2:(use*joint+joint),FUN = function(x){
  genVecQMC(p, (x))})
############################################################
###Running functions
############################################################
###Creating place to save time
RunningTime <- matrix(data = NA_real_,nrow = 5,ncol = 2)
RunningTime <- as.data.frame(RunningTime)
names(RunningTime) <- c("Full","50")
Results <- matrix(data = NA_real_,nrow = 5,ncol = 2)
Results <- as.data.frame(RunningTime)
names(Results) <- c("Full","50")

for(i in 1:5){
  
  ############################################################
  ###Full distribution
  RunningTime$Full[i]<-system.time(Results$Full[i] <- log(mvtNormQuasiMonteCarlo(latticeRule$primeP, X, Sigma, latticeRule$genVec)[[1]]))[[3]]
  ###n by n distribution
  RunningTime$R[i]<-system.time(Results$R[i] <- pmvNorm3(M = M,Sigma = Sigma,Corr=Corr,upper = X,use = 50,
                                                         lattice = latticelist50,joint = joint))[[3]]
}
