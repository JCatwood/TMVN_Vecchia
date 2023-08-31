rm(list = ls())
library(GpGp)
library(mvtnorm)
library(TruncatedNormal)
library(VeccTMVN)
# generate TMVN realization ---------------------
set.seed(123)
n1 <- 30
n2 <- 30
n <- n1 * n2
m <- 40
N <- 1e3
locs <- as.matrix(expand.grid((1:n1) / n1, (1:n2) / n2))
covparms <- c(1, 0.1, 0.01)
cov_mat <- matern15_isotropic(covparms, locs)
y <- as.vector(t(chol(cov_mat)) %*% rnorm(n))
b_censor <- 0
ind_censor <- which(y < b_censor)
ind_obs <- which(!(y < b_censor))
n_censor <- length(ind_censor)
n_obs <- n - n_censor
cov_name <- "matern15_isotropic"
# joint inference with VeccTMVN --------------------------
time_joint_Vecc <- system.time(
  samp_Vecc <- ptmvrandn(locs, ind_censor, y, b_censor, cov_name,
    covparms,
    m = m, N = N
  )
)[[3]]
mean_est <- rowMeans(samp_Vecc)
plot(y, mean_est)
# generate the north-west corner separately --------------------------
y_obs <- y[ind_obs]
locs_obs <- locs[ind_obs, , drop = F]
locs_censor <- locs[ind_censor, , drop = F]
locs_censor_mask <- (locs_censor[, 1]) < 0.6 & (locs_censor[, 2] > 0.4)
locs_censor <- locs_censor[locs_censor_mask, , drop = F]
locs_northwest <- rbind(locs_obs, locs_censor)
time_northwest_Vecc <- system.time(
  samp_Vecc_northwest <- ptmvrandn(
    locs_northwest,
    c((1 + n_obs):(n_obs + nrow(locs_censor))),
    c(y_obs, rep(NA, nrow(locs_censor))),
    b_censor, cov_name, covparms,
    m = m, N = N
  )
)[[3]]
# plot north-west corner samples --------------------------
mask_interest_northwest <- (locs_northwest[, 1]) < 0.5 &
  (locs_northwest[, 2] > 0.5)
samp_intest_northwest <-
  samp_Vecc_northwest[mask_interest_northwest, , drop = F]
mast_interest_all <- (locs[, 1] < 0.5) & (locs[, 2] > 0.5)
samp_interest_all <- samp_Vecc[mast_interest_all, , drop = F]
lim_qq <- c(min(samp_interest_all, samp_intest_northwest), 0)
if (!file.exists("plots")) {
  dir.create("plots")
}
pdf(file = "plots/PTMVN_sim_low.pdf", width = 15, height = 5)
par(mfrow = c(1, 3))
hist(samp_interest_all[samp_interest_all < 0],
  breaks = 30, main = "Global Sampling", cex.lab = 1.3,
  cex.axis = 1.3, xlab = NULL
)
hist(samp_intest_northwest[samp_intest_northwest < 0],
  breaks = 30, main = "Regional Sampling", cex.lab = 1.3,
  cex.axis = 1.3, xlab = NULL
)
plot(sort(samp_interest_all), sort(samp_intest_northwest),
  xlab = "Global Simulation",
  ylab = "Regional Simulation",
  xlim = lim_qq, ylim = lim_qq,
  cex.lab = 1.3,
  cex.axis = 1.3
)
abline(a = 0, b = 1, col = "red", lty = "dashed")
dev.off()
