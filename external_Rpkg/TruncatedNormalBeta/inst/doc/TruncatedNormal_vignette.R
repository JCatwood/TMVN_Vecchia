## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align="center", 
  fig.width = 5, 
  fig.height = 5
)

## ----setup, echo = FALSE------------------------------------------------------
library(TruncatedNormal)
par( bty = "l", pch = 20, xaxs = "i", yaxs = "i")
set.seed(0)

## ----simutrunc----------------------------------------------------------------
library(TruncatedNormal)
set.seed(1234)
sigma <- matrix(c(1,0.9,0.9,1), ncol = 2)
mu <- c(-3, 0)
u <- c(-6, Inf)
A <- matrix(c(1,-1,0,1), ncol = 2, byrow = TRUE)
# Sample truncated Gaussian variables and back-transforms
Y <- rtmvnorm(n = 1e2, mu = c(A %*% mu), sigma = A %*% sigma %*% t(A), ub = u)
X <- t(solve(A) %*% t(Y))
plot(X, panel.first = abline(a = 6, b = 1, col = 2),
     xlab = expression(x[1]), ylab = expression(x[2]), 
     xlim = c(-8,0), ylim = c(-5,5))
# Compare with unconstrained samples
points(rtmvnorm(n=1e2, mu = mu, sigma = sigma), col = 4) 

## ----rareproba----------------------------------------------------------------
d <- 1000
sigma <- 0.5 * (diag(d) + matrix(1, d, d))
est <- pmvnorm(sigma = sigma, lb = rep(0, d), type = "qmc", B = 1e4)
print(est)
#Compare est with exact value by computing relative error
abs(est - 1/(d+1))*(d+1)



## ----highdim,  fig.align="center", fig.width = 5, fig.height = 5--------------
d <- 60
sigma <- 0.1 * diag(d) + 0.9 * matrix(1, d, d)
l <- (1:d)/(4*d); u <- l + 2
X <- rtmvt(n = 1e4, sigma = sigma, lb = l, ub = u, df = 3)
boxplot(t(X) ~ as.factor(1:d), xlab = "dimension index", 
        ylab = expression(X["i"]))


## ----normqprec----------------------------------------------------------------
l <- 9; u <- 9.5
hist(rtnorm(n = 1e4, lb = l, ub = u), 
     xlim = c(9,9.5), xaxs = "i", main = "", xlab = "x")
# Now compare speed of the two methods
timing <- matrix(0, ncol = 2, nrow = 20)
for(i in 1:20){
  timing[i,] <- c(
    system.time(rtnorm(n = 1e5, lb = l, ub = u, method = "fast"))[3],
    system.time(rtnorm(n = 1e5, lb = l, ub = u, method = "invtransfo"))[3]
  )
}
colMeans(timing)

## ----probitneph---------------------------------------------------------------
# Exact Bayesian Posterior Simulation Example.

data("lupus"); # load lupus data
Y <- lupus[,1]; # response data
X <- as.matrix(lupus[,-1])  # construct design matrix
n <- nrow(X)
d <- ncol(X)
X <- diag(2*Y-1) %*% X; # incorporate response into design matrix
nusq <- 10000; # prior scale parameter
C <- solve(diag(d)/nusq + crossprod(X))
sigma <- diag(n) + nusq*tcrossprod(X) # this is covariance of Z given beta
est <- pmvnorm(sigma = sigma, lb = 0) 
# estimate acceptance probability of crude Monte Carlo
print(attributes(est)$upbnd/est[1])
# reciprocal of acceptance probability
Z <- rtmvnorm(sigma = sigma, n = 1e3, lb = rep(0, n))
 # sample exactly from auxiliary distribution
beta <- rtmvnorm(n = nrow(Z), sigma = C) + Z %*% X %*% C
 # simulate beta given Z and plot boxplots of marginals
boxplot(beta, ylab = expression(beta))
# plot the boxplots of the marginal distribution of the betas
print(colMeans(beta)) # output the posterior means
 

## ----tobit--------------------------------------------------------------------
data(mroz, package = "TruncatedNormal")
#Censored observations denote Yc, Yu for uncensored
Y <- mroz[,"whrs"]
X <- cbind(1, as.matrix(mroz[,-1]), I(mroz[,"exp"]^2))
n <- nrow(X); d <- ncol(X)
uncens <- Y > 0
Yu <- Y[uncens]; Yc <- Y[!uncens]
Xu <- X[uncens,]; Xc <- X[!uncens,]
invXtXu <- solve(crossprod(Xu))
sigma <- diag(nrow(Xc)) + Xc %*% invXtXu %*% t(Xc)
s <- sqrt(c(t(Yu) %*% (diag(nrow(Xu))- Xu %*% invXtXu %*% t(Xu)) %*% Yu))
# least squares residual variance estimate
nu <- nrow(Xu) - (d - 1) # degrees of freedom
beta_hat <- invXtXu %*% crossprod(Xu, Yu)
Yc_hat <- c(Xc %*% beta_hat) # fitted values
l <- sqrt(nu) * Yc_hat/s # upper threshold for censoring is zero
# simulate (Z,R) from a truncated Student
B <- 1e3
TR <- tregress(n = B, lb = l, ub = rep(Inf, length(l)), sigma = sigma, df = nu)
R <- TR$R
Z <- t(TR$Z)
# Reverse the mapping (beta,sigma) -> (Z,R) 
sig <- s/R # posterior of sigma
C <- solve(crossprod(Xu) + crossprod(Xc))
beta <- matrix(0, nrow = B, ncol = d)
for(i in 1:B){
    W <- Yc_hat - sig[i]*Z[,i] # auxiliary variables
    beta[i,] <- c(C %*% (crossprod(Xu, Yu) + crossprod(Xc, W))) + 
                    sig[i]*rtmvnorm(sigma = C, n = 1)
}
colnames(beta) <- colnames(X)
# Boxplots of the marginal posterior distribution
boxplot(beta[,-1], las = 2, ylab = expression(beta))
# Plot marginal means and standard deviations
summary(beta)

