rm(list = ls())

# Simulation -------------------------------
library(VeccTMVN)
library(tlrmvnmvt)
library(TruncatedNormal)
library(mvtnorm)
library(GpGp)
source("../funcs/CDFNormalAproxPackNoC_pmvn.R")
source("../funcs/utils.R")
## MVN prob gen funcs ----------------------------------------
prob1_gen <- function(n, d, retDenseCov = F) {
  locs <- grid_gen(n, d)$grid
  a <- rep(-Inf, n)
  b <- rep(0, n)
  cov_parms <- c(1.0, 0.1, 0.01)
  cov_name <- "matern15_isotropic"
  cov_mat <- get(cov_name)(cov_parms, locs)
  # odr <- TruncatedNormal::cholperm(cov_mat, a, b)$perm
  # a <- a[odr]
  # b <- b[odr]
  # locs <- locs[odr, , drop = F]
  # cov_mat <- get(cov_name)(cov_parms, locs)
  return(list(
    a = a, b = b, locs = locs, cov_parms = cov_parms,
    cov_name = cov_name, cov_mat = cov_mat
  ))
}
prob2_gen <- function(n, d, retDenseCov = F) {
  locs <- latin_gen(n, d)
  a <- rep(-Inf, n)
  b <- -runif(n, 0, 2)
  cov_parms <- c(1.0, 0.1, 0.01)
  cov_name <- "matern15_isotropic"
  cov_mat <- get(cov_name)(cov_parms, locs)
  return(list(
    a = a, b = b, locs = locs, cov_parms = cov_parms,
    cov_name = cov_name, cov_mat = cov_mat
  ))
}
prob3_gen <- function(n, d, retDenseCov = F) {
  locs <- grid_gen(n, d)$grid
  a <- rep(-1, n)
  b <- rep(1, n)
  cov_parms <- c(1.0, 0.1, 0.01)
  cov_name <- "matern15_isotropic"
  cov_mat <- get(cov_name)(cov_parms, locs)
  return(list(
    a = a, b = b, locs = locs, cov_parms = cov_parms,
    cov_name = cov_name, cov_mat = cov_mat
  ))
}
## Prob setups -----------------------
set.seed(123)
n <- 900
d <- 2
m_vec <- seq(from = 10, to = 90, by = 20)
N_SOV_TLR <- seq(from = 1e4, to = 9e4, by = 2e4)
N_TN_VCDF <- seq(from = 6e3, to = 1e4, by = 1e3)
prob_ind <- 3
order_mtd <- 3
prob_obj <- get(paste0("prob", prob_ind, "_gen"))(n, d, retDenseCov = T)
a <- prob_obj$a
b <- prob_obj$b
locs <- prob_obj$locs
cov_mat <- prob_obj$cov_mat
cov_name <- prob_obj$cov_name
cov_parms <- prob_obj$cov_parms
z_order <- tlrmvnmvt::zorder(locs)
## Iteratively compute the same MVN prob -----------------------
niter <- 30
# time_df <- data.frame(matrix(
#   NA, niter,
#   2 * length(N_SOV_TLR) + 2 * length(N_TN_VCDF) + length(m_vec)
# ))
# prob_df <- data.frame(matrix(
#   NA, niter,
#   2 * length(N_SOV_TLR) + 2 * length(N_TN_VCDF) + length(m_vec)
# ))
# for (i in 1:niter) {
#   est_Vecc <- rep(NA, length(m_vec))
#   time_Vecc <- rep(NA, length(m_vec))
#   for (j in 1:length(m_vec)) {
#     cat("VeccTMVN", j, "\n")
#     ### Compute MVN prob with idea V -----------------------
#     m <- m_vec[j]
#     time_Vecc[j] <- system.time(est_Vecc[j] <- VeccTMVN::pmvn(a, b, 0,
#       locs = locs, covName = cov_name,
#       reorder = order_mtd, covParms = cov_parms,
#       m = m, verbose = T,
#       NLevel1 = 10, NLevel2 = 1e3, m_ord = m
#     ))[[3]]
#   }
#   ### Compute MVN prob with other methods -----------------------
#   est_TLR <- rep(NA, length(N_SOV_TLR))
#   time_TLR <- rep(NA, length(N_SOV_TLR))
#   for (j in 1:length(N_SOV_TLR)) {
#     cat("TLR", j, "\n")
#     N <- N_SOV_TLR[j]
#     err_obj <- try(
#       time_TLR[j] <- system.time(
#         est_TLR[j] <- tlrmvnmvt::pmvn(a[z_order], b[z_order],
#           sigma = cov_mat[z_order, z_order],
#           algorithm = tlrmvnmvt::TLRQMC(N = round(N / 20), m = sqrt(n), epsl = 1e-6)
#         )
#       )[[3]]
#     )
#     if (class(err_obj) == "try-error") {
#       time_TLR[j] <- NA
#       est_TLR[j] <- NA
#     }
#   }
#   est_SOV <- rep(NA, length(N_SOV_TLR))
#   time_SOV <- rep(NA, length(N_SOV_TLR))
#   for (j in 1:length(N_SOV_TLR)) {
#     cat("SOV", j, "\n")
#     N <- N_SOV_TLR[j]
#     time_SOV[j] <- system.time(
#       est_SOV[j] <- tlrmvnmvt::pmvn(a, b,
#         sigma = cov_mat,
#         algorithm = tlrmvnmvt::GenzBretz(N = round(N / 20))
#       )
#     )[[3]]
#   }
#   est_TN <- rep(NA, length(N_TN_VCDF))
#   time_TN <- rep(NA, length(N_TN_VCDF))
#   for (j in 1:length(N_TN_VCDF)) {
#     cat("TN", j, "\n")
#     N <- N_TN_VCDF[j]
#     time_TN[j] <- system.time(est_TN[j] <- TruncatedNormal::pmvnorm(
#       rep(0, n), cov_mat,
#       lb = a, ub = b, B = N
#     ))[[3]]
#   }
#   est_Nascimento <- rep(NA, length(N_TN_VCDF))
#   time_Nascimento <- rep(NA, length(N_TN_VCDF))
#   for (j in 1:length(N_TN_VCDF)) {
#     cat("VCDF", j, "\n")
#     N <- N_TN_VCDF[j]
#     err_obj <- try(
#       time_Nascimento[j] <- system.time(
#         est_Nascimento[j] <- exp(CDFNormalAproxPackNoC.pmvn(b,
#           mean = rep(0, n),
#           sigma = cov_mat, p = round(N / 20)
#         ))
#       )[[3]]
#     )
#     if (class(err_obj) == "try-error") {
#       time_Nascimento[j] <- NA
#       est_Nascimento[j] <- NA
#     }
#   }
#   ### save results ------------------------
#   time_df[i, ] <- c(time_Vecc, time_TN, time_TLR, time_SOV, time_Nascimento)
#   prob_df[i, ] <- c(est_Vecc, est_TN, est_TLR, est_SOV, est_Nascimento)
# }

# TB removed
load(paste0(
  "results/Vecc_bias_lowdim_exp",
  prob_ind, "_", order_mtd, ".RData"
))
# TB removed

est_bmark <- rep(NA, 100)
for (j in 1:100) {
  cat("Benchmark TN iter", j, "\n")
  est_bmark[j] <- TruncatedNormal::pmvnorm(
    rep(0, n), cov_mat,
    lb = a, ub = b, B = max(N_TN_VCDF) * 10
  )
}
if (!file.exists("results")) {
  dir.create("results")
}
save(m_vec, N_SOV_TLR, N_TN_VCDF, time_df, prob_df, est_bmark, file = paste0(
  "results/Vecc_bias_lowdim_exp",
  prob_ind, "_", order_mtd, ".RData"
))

# Plotting -----------------------------------
load(paste0(
  "results/Vecc_bias_lowdim_exp",
  prob_ind, "_", order_mtd, ".RData"
))
library(ggplot2)
library(tidyr)
library(scales)
time_vs_err_plt <- function(probDf, timeDf, benchmark = NULL) {
  col_names <- c(
    paste("VMET", m_vec, sep = "_"),
    paste("MET", N_TN_VCDF, sep = "_"),
    paste("TLR", N_SOV_TLR, sep = "_"),
    paste("SOV", N_SOV_TLR, sep = "_"),
    paste("VCDF", N_TN_VCDF, sep = "_")
  )
  colnames(probDf) <- col_names
  colnames(timeDf) <- col_names
  colname_benchmark <- paste("MET", max(N_TN_VCDF), sep = "_")
  if (is.null(benchmark)) {
    benchmark <- mean(probDf[[colname_benchmark]])
  } else {
    benchmark <- mean(benchmark)
  }
  rmse <- sqrt(colMeans((probDf - benchmark)^2))
  secs <- colMeans(timeDf)
  mydf <- data.frame(
    time = secs, rmse = rmse,
    method = gsub("[_0-9]*", "", col_names)
  )
  # log2 trans to x and y ticks are applied by default
  breaks <- exp(seq(
    from = log(min(mydf$rmse) / 2),
    to = log(max(mydf$rmse)), length.out = 5
  ))
  ggplot(data = mydf, mapping = aes(x = time, y = rmse)) +
    geom_point(mapping = aes(colour = method, shape = method), size = 3) +
    geom_line(mapping = aes(group = method, colour = method)) +
    scale_y_continuous(
      name = "RMSE", trans = "log2", breaks = breaks,
      labels = signif(breaks, digits = 2)
    ) +
    scale_x_continuous(
      name = "time (seconds)", trans = "log2"
    ) +
    ggtitle(paste("Scenario", prob_ind)) +
    theme(
      text = element_text(size = 14),
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
}
if (!file.exists("plots")) {
  dir.create("plots")
}
time_vs_err_plt(prob_df, time_df, est_bmark)
ggsave(
  paste0("plots/err_vs_time_lowdim_", prob_ind, "_", order_mtd, ".pdf"),
  width = 5,
  height = 5
)
