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
ggplot(data = data.frame(val = as.vector(samp_Vecc_joint))) +
  geom_histogram(
    mapping = aes(x = val),
    breaks = seq(from = -6.2, to = 0, by = 0.2)
  ) +
  scale_x_continuous(name = "VMET Samples") +
  scale_y_continuous(name = "Frequency") +
  theme(
    text = element_text(size = 14), legend.position = "none"
  )
ggsave(paste0("plots/TMVN_sim_low_vecc_joint.pdf"),
  width = 5,
  height = 5
)

ggplot(data = data.frame(val = as.vector(samp_TN))) +
  geom_histogram(
    mapping = aes(x = val),
    breaks = seq(from = -6.2, to = 0, by = 0.2)
  ) +
  scale_x_continuous(name = "MET Samples") +
  scale_y_continuous(name = "Frequency") +
  theme(
    text = element_text(size = 14), legend.position = "none"
  )
ggsave(paste0("plots/TMVN_sim_low_TN.pdf"),
  width = 5,
  height = 5
)

ggplot(data = data.frame(
  x = sort(samp_Vecc_joint)[seq(from = 1, to = N * n, by = 10)],
  y = sort(samp_TN)[seq(from = 1, to = N * n, by = 10)]
)) +
  geom_point(mapping = aes(x = x, y = y)) +
  scale_x_continuous(name = "VMET samples", limits = c(range(samp_TN))) +
  scale_y_continuous(name = "MET samples", limits = c(range(samp_TN))) +
  theme(
    text = element_text(size = 14), legend.position = "none"
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")
ggsave(paste0("plots/TMVN_sim_low_qqplot.pdf"),
  width = 5,
  height = 5
)
