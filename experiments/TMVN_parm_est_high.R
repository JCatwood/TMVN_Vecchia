rm(list = ls())

library(GpGp)
library(mvtnorm)
library(TruncatedNormal)
library(VeccTMVN)
set.seed(123)
n1 <- 80
n2 <- 80
n <- n1 * n2
locs <- as.matrix(expand.grid((1:n1) / n1, (1:n2) / n2))
covparms <- c(1, 0.1, 0.03)
cov_mat <- matern15_isotropic(covparms, locs)
y <- as.vector(t(chol(cov_mat)) %*% rnorm(n))
b_censor <- 0
ind_censor <- which(y < b_censor)
ind_obs <- which(!(y < b_censor))
y_aug <- y
y_aug[ind_censor] <- 0
n_censor <- length(ind_censor)
n_obs <- n - n_censor
cov_name <- "matern15_isotropic"

# reorder if needed ---------------------------
# y_obs <- y[ind_obs]
# y_censor <- y[ind_censor] # also need to adjust b_censor if needed
# locs_obs <- locs[ind_obs, , drop = F]
# locs_censor <- locs[ind_censor, , drop = F]
# y <- c(y_obs, y_censor)
# locs <- rbind(locs_obs, locs_censor)
# ind_censor <- (1:length(ind_censor)) + length(ind_obs)
# ind_obs <- 1:length(ind_obs)
# odr <- Vecc_reorder(c(y_obs, rep(-Inf, n_censor)),
#                     c(y_obs, rep(b_censor, n_censor)),
#   30,
#   locs = locs, covName = cov_name,
#   covParms = covparms
# )$order
# if (any(odr[(n_obs + 1):n] <= n_obs)) {
#   warning("Vecchia_univar failed\n")
# } else {
#   locs <- locs[odr, , drop = F]
#   y <- y[odr]
# }

# compute the likelihood surface ------------------------------------
# ranges <- seq(from = 0.05, to = 0.25, by = 0.005)
# loglk_Vecc_vec <- rep(0, length(ranges))
# loglk_GPGP_vec <- rep(0, length(ranges))
# idx <- 1
# NN <- GpGp::find_ordered_nn(locs, m = 50)
# for (myrange in ranges) {
#   covparms[2] <- myrange
#   set.seed(123)
#   loglk_Vecc_vec[idx] <- loglk_censor_MVN(
#     locs, ind_censor, y, b_censor, cov_name,
#     covparms,
#     m = 50
#   )
#   loglk_GPGP_vec[idx] <- GpGp::vecchia_meanzero_loglik(
#     covparms, cov_name, y_aug,
#     locs,
#     NNarray = NN
#   )$loglik
#   idx <- idx + 1
# }
# save(ranges, loglk_Vecc_vec, loglk_GPGP_vec,
#   file = "results/censored_normal_sim.RData"
# )

# plot the likelihood surface ------------------------------------
load(paste0("results/censored_normal_sim.RData"))
library(ggplot2)
library(tidyr)
library(scales)
mycol <- hue_pal()(2)
mydf <- data.frame(
  range = ranges,
  loglk_censor_MVN = loglk_Vecc_vec - min(loglk_Vecc_vec),
  loglk_GP = loglk_GPGP_vec - min(loglk_GPGP_vec)
)
mydf <- pivot_longer(mydf, c(2:3), names_to = "model", values_to = "loglk")
max_range_CMVN <- ranges[which.max(loglk_Vecc_vec)]
max_range_GP <- ranges[which.max(loglk_GPGP_vec)]
ggplot(mydf, aes(x = range, y = loglk, group = model, color = model)) +
  geom_line() +
  geom_vline(xintercept = 0.1, linetype = "dotdash", color = "black") +
  geom_vline(xintercept = max_range_CMVN, linetype = "dashed", color = mycol[1]) +
  geom_vline(xintercept = max_range_GP, linetype = "dashed", color = mycol[2]) +
  scale_x_continuous(name = "Range Parameter") +
  scale_y_continuous(name = "Log-likelihood") +
  scale_color_discrete(labels = c("censored MVN", "GP")) +
  theme(
    legend.position = c(0.8, 0.9), legend.text = element_text(size = 12),
    legend.title = element_blank(), text = element_text(size = 14)
  )
if (!file.exists("plots")) {
  dir.create("plots")
}
ggsave(paste0("plots/censored_normal_sim.pdf"), width = 5, height = 5)
