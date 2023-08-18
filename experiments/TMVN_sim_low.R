rm(list = ls())
library(GpGp)
library(mvtnorm)
library(TruncatedNormal)
library(VeccTMVN)
set.seed(123)
# input generation
n1 <- 30
n2 <- 30
n <- n1 * n2
locs <- as.matrix(expand.grid((1:n1) / n1, (1:n2) / n2))
covparms <- c(1.0, 0.1, 0.01)
mu <- rep(0, n)
N <- 1000
cov_mat <- matern15_isotropic(covparms, locs)
a <- rep(-Inf, n)
b <- rep(-0, n)
time_Vecc_joint <- system.time(samp_Vecc_joint <- mvrandn(
  a, b, mu, locs, "matern15_isotropic", covparms,
  m = 30, N = N, verbose = T
))[[3]]
time_TN <- system.time(samp_TN <- TruncatedNormal::mvrandn(
  a, b, cov_mat,
  n = N, mu = mu
))[[3]]
if (!file.exists("results")) {
  dir.create("results")
}
save(n1, n2, time_TN, time_Vecc_joint, samp_Vecc_joint, samp_TN,
  file = "results/TMVN_sim_low.RData"
)
##  histogram for verification -------------------
par(mfrow = c(1, 2))
hist(samp_Vecc_joint, breaks = 30, main = "Vecc Samples")
hist(samp_TN, breaks = 30, main = "TN Samples")
image(matrix(samp_TN[, 1], n1, n2))
image(matrix(samp_Vecc_joint[, 1], n1, n2))
