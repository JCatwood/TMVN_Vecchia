rm(list = ls())
library(GpGp)
library(mvtnorm)
library(TruncatedNormal)
library(VeccTMVN)
library(CensSpatial)
# generate TMVN realization ---------------------
set.seed(321)
# size of a 2D grid, with n1 horizontal and n2 vertical bands
n1 <- 30
n2 <- 30
n <- n1 * n2
m <- 30 # the tuning parameter controling conditioning set sizes
N <- 1e3 # number of samples to generate
locs <- as.matrix(expand.grid((1:n1) / n1, (1:n2) / n2))
covparms <- c(1, 0.1, 0.01)
cov_mat <- matern15_isotropic(covparms, locs)
y <- as.vector(t(chol(cov_mat)) %*% rnorm(n))
b_censor <- 0 # level of detection, below which the responses are censored
ind_censor <- which(y < b_censor)
ind_obs <- which(!(y < b_censor))
n_censor <- length(ind_censor)
n_obs <- n - n_censor
cov_name <- "matern15_isotropic"
# joint inference with VeccTMVN --------------------------
# generate samples underlying the censored responses with VMET and time it
time_joint_Vecc <- system.time(
  samp_Vecc <- ptmvrandn(locs, ind_censor, y, b_censor, cov_name,
    covparms,
    m = m, N = N
  )
)[[3]]
mean_est <- rowMeans(samp_Vecc) # mean prediction
plot(y, mean_est)
# generate the north-west corner separately --------------------------
# next generate samples underlying the censored responses in the north-west corner
# conditional on all observed responses
y_obs <- y[ind_obs]
locs_obs <- locs[ind_obs, , drop = F]
locs_censor <- locs[ind_censor, , drop = F]
locs_censor_mask <- (locs_censor[, 1]) < 0.6 & (locs_censor[, 2] > 0.4)
# extract censored responses in the north-west corner
locs_censor <- locs_censor[locs_censor_mask, , drop = F]
locs_northwest <- rbind(locs_obs, locs_censor)
y_censor <- y[ind_censor]
y_censor <- y_censor[locs_censor_mask]
# generate samples underlying the censored responses with VMET and time it
time_northwest_Vecc <- system.time(
  samp_Vecc_northwest <- ptmvrandn(
    locs_northwest,
    c((1 + n_obs):(n_obs + nrow(locs_censor))),
    c(y_obs, rep(NA, nrow(locs_censor))),
    b_censor, cov_name, covparms,
    m = m, N = N
  )
)[[3]]
# compare the RMSE computed using all observed responses or only the observed
# responses in the north-west corner
RMSE_region <- sqrt(mean((y_censor -
  rowMeans(samp_Vecc_northwest)[
    (1 + n_obs):(n_obs + nrow(locs_censor))
  ])^2))
RMSE_global <- sqrt(mean((y_censor -
  rowMeans(samp_Vecc)[ind_censor][locs_censor_mask])^2))

# predict with level-of-detection supplemented GP ---------------------
RMSE_GP_aug <- sqrt(mean((y_censor - b_censor)^2))

# predict with SAEMSCL -------------------------------
# draw samples using the CensSpatial package, which I only installed on my Linux VM
# the following code may not run on MacOS as CensSpatial is hard to install on Mac
y_aug <- c(y_obs, rep(b_censor, nrow(locs_censor)))
cc <- c(rep(F, length(y_obs)), rep(T, nrow(locs_censor)))
time_SAEMSCL <- system.time(est_SAEMSCL <- SAEMSCL(cc, y_aug,
  cens.type = "left", trend = "other", x = matrix(1, length(y_aug), 1),
  coords = locs_northwest,
  kappa = 1.5, M = 500,
  perc = 0.25, MaxIter = 10, pc = 0.2, cov.model = "matern",
  fix.nugget = TRUE, nugget = covparms[3],
  inits.sigmae = covparms[1], inits.phi = covparms[2], search = TRUE,
  lower = covparms[2], upper = covparms[2] + 1e-4
))[3]
time_SAEMSCL_pred <- system.time(
  pred_SAEMSCL <- predSCL(
    matrix(1, nrow(locs_censor), 1),
    locs_censor, est_SAEMSCL
  )
)[3]
RMSE_GP_SAEMSCL <- sqrt(mean((y_censor - pred_SAEMSCL$prediction)^2))
save(time_SAEMSCL, time_SAEMSCL_pred, est_SAEMSCL, pred_SAEMSCL, RMSE_GP_SAEMSCL,
  file = "results/PTMVN_sim_low.RData"
)
load("results/PTMVN_sim_low.RData")
# save the results. Once you have the results, you can comment out the 
# predict with SAEMSCL section to be faster

# plot north-west corner samples --------------------------
library(ggplot2)
mask_interest_northwest <- (locs_northwest[, 1]) < 0.5 &
  (locs_northwest[, 2] > 0.5)
samp_intest_northwest <-
  samp_Vecc_northwest[mask_interest_northwest, , drop = F]
mask_interest_all <- (locs[, 1] < 0.5) & (locs[, 2] > 0.5)
samp_interest_all <- samp_Vecc[mask_interest_all, , drop = F]
lim_qq <- c(min(samp_interest_all, samp_intest_northwest), 0)
if (!file.exists("plots")) {
  dir.create("plots")
}
ggplot(data = data.frame(val = samp_interest_all[samp_interest_all < 0])) +
  geom_histogram(
    mapping = aes(x = val),
    breaks = seq(from = -2.6, to = 0, by = 0.1)
  ) +
  scale_x_continuous(name = "Global Sampling") +
  scale_y_continuous(name = "Frequency") +
  theme(
    text = element_text(size = 14), legend.position = "none"
  )
ggsave(paste0("plots/PTMVN_sim_low_global.pdf"),
  width = 5,
  height = 5
)

ggplot(data = data.frame(val = samp_intest_northwest[samp_intest_northwest < 0])) +
  geom_histogram(
    mapping = aes(x = val),
    breaks = seq(from = -2.6, to = 0, by = 0.1)
  ) +
  scale_x_continuous(name = "Regional Sampling") +
  scale_y_continuous(name = "Frequency") +
  theme(
    text = element_text(size = 14), legend.position = "none"
  )
ggsave(paste0("plots/PTMVN_sim_low_region.pdf"),
  width = 5,
  height = 5
)

ggplot(data = data.frame(
  x = sort(samp_intest_northwest[samp_intest_northwest < 0]),
  y = sort(samp_interest_all[samp_interest_all < 0])
)) +
  geom_point(mapping = aes(x = x, y = y)) +
  scale_x_continuous(
    name = "Regional samples",
    limits = c(range(samp_interest_all[samp_interest_all < 0]))
  ) +
  scale_y_continuous(
    name = "Global samples",
    limits = c(range(samp_interest_all[samp_interest_all < 0]))
  ) +
  theme(
    text = element_text(size = 14), legend.position = "none"
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")
ggsave(paste0("plots/PTMVN_sim_low_qqplot.pdf"),
  width = 5,
  height = 5
)
# compute RMSE -----------------------------------------
cat("RMSE for global sim", RMSE_global, "\n")
cat("RMSE for regional sim", RMSE_region, "\n")
cat("RMSE for LOD-GP", RMSE_GP_aug, "\n")
cat("RMSE for SAEMSCL", RMSE_GP_SAEMSCL, "\n")
