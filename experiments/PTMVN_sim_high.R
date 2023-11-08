rm(list = ls())
library(GpGp)
library(mvtnorm)
library(TruncatedNormal)
library(VeccTMVN)
library(CensSpatial)
# generate TMVN realization ---------------------
set.seed(123)
n1 <- 80
n2 <- 80
n <- n1 * n2
n_test <- 500
m <- 40
N <- 1e3
locs <- as.matrix(expand.grid((1:n1) / n1, (1:n2) / n2))
locs_test <- matrix(runif(n_test * 2), n_test, 2)
locs_test[, 1] <- locs_test[, 1] * 0.5
locs_test[, 2] <- locs_test[, 2] * 0.5 + 0.5
covparms <- c(1, 0.1, 0.03)
cov_mat <- matern15_isotropic(covparms, rbind(locs, locs_test))
y_all <- as.vector(t(chol(cov_mat)) %*% rnorm(n + n_test))
y_test <- y_all[(n + 1):(n + n_test)]
y <- y_all[1:n]
b_censor <- 0
ind_censor <- which(y < b_censor)
ind_obs <- which(!(y < b_censor))
n_censor <- length(ind_censor)
n_obs <- n - n_censor
cov_name <- "matern15_isotropic"
# generate the north-west corner separately --------------------------
y_obs <- y[ind_obs]
y_censor <- y[ind_censor]
locs_obs <- locs[ind_obs, , drop = F]
locs_censor <- locs[ind_censor, , drop = F]
locs_censor_mask <- (locs_censor[, 1]) < 0.6 & (locs_censor[, 2] > 0.4)
locs_censor <- locs_censor[locs_censor_mask, , drop = F]
y_censor <- y_censor[locs_censor_mask]
locs_northwest <- rbind(locs_obs, locs_censor)
y_northwest <- c(y_obs, y_censor)
y_northwest_aug <- y_northwest
y_northwest_aug[(n_obs + 1):length(y_northwest)] <- b_censor
# time_northwest_Vecc <- system.time(
#   samp_Vecc_northwest <- ptmvrandn(
#     locs_northwest,
#     c((1 + n_obs):(n_obs + nrow(locs_censor))),
#     c(y_obs, rep(NA, nrow(locs_censor))),
#     b_censor, cov_name, covparms,
#     m = m, N = N
#   )
# )[[3]]
# # save results ---------------------------------
# if (!file.exists("results")) {
#   dir.create("results")
# }
# save(time_northwest_Vecc, locs_northwest, samp_Vecc_northwest,
#   locs_test, y_test, locs_obs, y_obs, covparms, n_test, n_obs,
#   file = "results/PTMVN_sim_high.RData"
# )
# generate samples with SAEMSCL ----------------------------
# y_aug <- c(y_obs, rep(b_censor, nrow(locs_censor)))
# cc <- c(rep(F, length(y_obs)), rep(T, nrow(locs_censor)))
# time_SAEMSCL <- system.time(est_SAEMSCL <- SAEMSCL(cc, y_aug,
#   cens.type = "left", trend = "other", x = matrix(1, length(y_aug), 1),
#   coords = locs_northwest,
#   kappa = 1.5, M = 50,
#   perc = 0.25, MaxIter = 10, pc = 0.2, cov.model = "matern",
#   fix.nugget = F, nugget = covparms[3] + 0.5,
#   inits.sigmae = covparms[1], inits.phi = covparms[2], search = TRUE,
#   lower = covparms[2], upper = covparms[2] + 1e-4
# ))[3]  # not scalable to this problem size
# save(time_SAEMSCL, est_SAEMSCL,
#   file = "results/PTMVN_sim_high_SAEMSCL.RData"
# )
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
