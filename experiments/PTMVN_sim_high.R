rm(list = ls())
library(GpGp)
library(mvtnorm)
library(TruncatedNormal)
library(VeccTMVN)
library(CensSpatial)
# generate TMVN realization ---------------------
set.seed(123)
# size of a 2D grid, with n1 horizontal and n2 vertical bands
n1 <- 80
n2 <- 80
n <- n1 * n2
n_test <- 500 # size of testing dataset
m <- 40 # the tuning parameter controling conditioning set sizes
N <- 1e3 # number of samples to generate
locs <- as.matrix(expand.grid((1:n1) / n1, (1:n2) / n2))
# generate testing locations from [0, 0.5] X [0.5, 1.0]
locs_test <- matrix(runif(n_test * 2), n_test, 2)
locs_test[, 1] <- locs_test[, 1] * 0.5
locs_test[, 2] <- locs_test[, 2] * 0.5 + 0.5
covparms <- c(1, 0.1, 0.03)
cov_mat <- matern15_isotropic(covparms, rbind(locs, locs_test))
y_all <- as.vector(t(chol(cov_mat)) %*% rnorm(n + n_test))
y_test <- y_all[(n + 1):(n + n_test)]
y <- y_all[1:n] # training responses
b_censor <- 0 # level of detection, below which the responses are censored
ind_censor <- which(y < b_censor)
ind_obs <- which(!(y < b_censor))
n_censor <- length(ind_censor)
n_obs <- n - n_censor
cov_name <- "matern15_isotropic"
# generate the north-west corner separately --------------------------
# extract all observed responses
y_obs <- y[ind_obs]
y_censor <- y[ind_censor]
locs_obs <- locs[ind_obs, , drop = F]
locs_censor <- locs[ind_censor, , drop = F]
# extract the censored responses only from [0, 0.6] X [0.4, 1.0], which envelops
# the testing responses in [0, 0.5] X [0.5, 1.0]
locs_censor_mask <- (locs_censor[, 1]) < 0.6 & (locs_censor[, 2] > 0.4)
locs_censor <- locs_censor[locs_censor_mask, , drop = F]
y_censor <- y_censor[locs_censor_mask]
# combine all observed responses and the censored responses from [0, 0.6] X [0.4, 1.0] 
locs_northwest <- rbind(locs_obs, locs_censor)
y_northwest <- c(y_obs, y_censor)
# substitute the censored responses by the level of detection to be used by our 
# competitor method
y_northwest_aug <- y_northwest
y_northwest_aug[(n_obs + 1):length(y_northwest)] <- b_censor
# draw samples of partially censored GP ---------------------
time_northwest_Vecc <- system.time(
  samp_Vecc_northwest <- ptmvrandn(
    locs_northwest,
    c((1 + n_obs):(n_obs + nrow(locs_censor))),
    c(y_obs, rep(NA, nrow(locs_censor))),
    b_censor, cov_name, covparms,
    m = m, N = N
  )
)[[3]]
if (!file.exists("results")) {
  dir.create("results")
}
save(time_northwest_Vecc, locs_northwest, samp_Vecc_northwest,
  locs_test, y_test, locs_obs, y_obs, covparms, n_test, n_obs,
  file = "results/PTMVN_sim_high.RData"
)
# save the results. Once you have the results, you can comment out the 
# draw samples of partially censored GP section to be faster

# generate samples with SAEMSCL ----------------------------
# draw samples using the CensSpatial package, which I only installed on my Linux VM
# the following code may not run on MacOS as CensSpatial is hard to install on Mac
y_aug <- c(y_obs, rep(b_censor, nrow(locs_censor)))
cc <- c(rep(F, length(y_obs)), rep(T, nrow(locs_censor)))
time_SAEMSCL <- system.time(est_SAEMSCL <- SAEMSCL(cc, y_aug,
  cens.type = "left", trend = "other", x = matrix(1, length(y_aug), 1),
  coords = locs_northwest,
  kappa = 1.5, M = 50,
  perc = 0.25, MaxIter = 10, pc = 0.2, cov.model = "matern",
  fix.nugget = F, nugget = covparms[3] + 0.5,
  inits.sigmae = covparms[1], inits.phi = covparms[2], search = TRUE,
  lower = covparms[2], upper = covparms[2] + 1e-4
))[3]  # not scalable to this problem size
save(time_SAEMSCL, est_SAEMSCL,
  file = "results/PTMVN_sim_high_SAEMSCL.RData"
)
# save the results. Once you have the results, you can comment out the 
# generate samples with SAEMSCL section to be faster

# plot north-west corner samples --------------------------
load("results/PTMVN_sim_high.RData")
mask_interest_northwest <- (locs_northwest[, 1]) < 0.5 &
  (locs_northwest[, 2] > 0.5)
samp_intest_northwest <-
  samp_Vecc_northwest[mask_interest_northwest, , drop = F]
y_interest_northwest <- y_northwest[mask_interest_northwest]
hist(samp_intest_northwest)
plot(rowMeans(samp_intest_northwest), y_interest_northwest)
# predict at locs_test ----------------------
# prediction with ordinary kriging based on samples generated by VMET,
# CensSpatial, and level-of-detection supplemented GP. LOD GP doesn't have 
# randomness, hence generating only one sample
cov_mat <- matern15_isotropic(covparms, rbind(locs_obs, locs_test))
pred_GP <- as.vector(cov_mat[(n_obs + 1):(n_obs + n_test), 1:n_obs] %*%
  solve(cov_mat[1:n_obs, 1:n_obs], y_obs))
n_northwest <- nrow(locs_northwest)
cov_mat <- matern15_isotropic(covparms, rbind(locs_northwest, locs_test))
pred_GP_aug <- as.vector(cov_mat[(n_northwest + 1):(n_northwest + n_test), 1:n_northwest] %*%
  solve(cov_mat[1:n_northwest, 1:n_northwest], y_northwest_aug))
cov_mat <- matern15_isotropic(covparms, rbind(locs_northwest, locs_test))

pred_PTMVN <- rowMeans(
  cov_mat[(n_northwest + 1):(n_northwest + n_test), 1:n_northwest] %*%
    solve(cov_mat[1:n_northwest, 1:n_northwest], samp_Vecc_northwest)
)
sqrt(mean((y_test - pred_GP)^2))
sqrt(mean((y_test - pred_PTMVN)^2))
sqrt(mean((y_test - pred_GP_aug)^2))
