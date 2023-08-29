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
# time_Vecc_joint <- system.time(samp_Vecc_joint <- mvrandn(
#   a, b, mu, locs, "matern15_isotropic", covparms,
#   m = 30, N = N, verbose = T
# ))[[3]]
# time_TN <- system.time(samp_TN <- TruncatedNormal::mvrandn(
#   a, b, cov_mat,
#   n = N, mu = mu
# ))[[3]]
# if (!file.exists("results")) {
#   dir.create("results")
# }
# save(n1, n2, time_TN, time_Vecc_joint, samp_Vecc_joint, samp_TN,
#   file = "results/TMVN_sim_low.RData"
# )
##  histogram for verification -------------------
load("results/TMVN_sim_low.RData")
if (!file.exists("plots")) {
  dir.create("plots")
}
pdf(file = "plots/TMVN_sim_low.pdf", width = 15, height = 5)
par(mfrow = c(1, 3))
hist(samp_Vecc_joint,
  breaks = 30, main = "VMET Samples", cex.lab = 1.3,
  cex.axis = 1.3
)
hist(samp_TN,
  breaks = 30, main = "MET Samples", cex.lab = 1.3,
  cex.axis = 1.3
)
lim_qq <- range(samp_TN)
plot(sort(samp_Vecc_joint)[seq(from = 1, to = N * n, by = 10)],
  sort(samp_TN)[seq(from = 1, to = N * n, by = 10)],
  xlab = "VMET samples",
  ylab = "MET samples",
  xlim = lim_qq, ylim = lim_qq,
  cex.lab = 1.3,
  cex.axis = 1.3
)
abline(a = 0, b = 1, col = "red", lty = "dashed")
dev.off()
image(matrix(samp_TN[, 1], n1, n2))
image(matrix(samp_Vecc_joint[, 1], n1, n2))
