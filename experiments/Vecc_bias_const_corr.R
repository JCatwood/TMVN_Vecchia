rm(list = ls())

# Simulation -------------------------------
library(VeccTMVN)
library(tlrmvnmvt)
library(TruncatedNormal)
library(mvtnorm)
library(GpGp)
library(scales)
source("../funcs/CDFNormalAproxPackNoC_pmvn.R")
## MVN prob gen funcs ----------------------------------------
#' Generate MVN probabilities with a constant covariace matrix 
#' @param n number of locations, also the dimension of the MVN probability
#' @param rho the constant correlation
#' @return a list containing:
#' \\item{a}{lower integration limits}
#' \\item{b}{upper integration limits}
#' \\item{covmat}{the generate covariance matrix for the MVN probability}
prob_gen <- function(n, rho = 0.5) {
  a <- rep(-Inf, n)
  b <- rep(0, n)
  covmat <- matrix(rho, n, n)
  diag(covmat) <- 1.0
  return(list(
    a = a, b = b, covmat = covmat
  ))
}
## Prob setups -----------------------
set.seed(123)
n <- 900
m_vec <- seq(from = 10, to = 100, by = 10) # see the change of bias as m increases
# see the change of bias as Monte Carlo sample size increases
N_SOV_TLR <- seq(from = 1e4, to = 9e4, by = 2e4)
N_TN_VCDF <- seq(from = 6e3, to = 1e4, by = 1e3)
prob_obj <- prob_gen(n)
a <- prob_obj$a
b <- prob_obj$b
covmat <- prob_obj$covmat
## Iteratively compute the same MVN prob -----------------------
niter <- 30
time_df <- data.frame(matrix(
  NA, niter,
  2 * length(N_SOV_TLR) + 2 * length(N_TN_VCDF) + length(m_vec) * 2
))
prob_df <- data.frame(matrix(
  NA, niter,
  2 * length(N_SOV_TLR) + 2 * length(N_TN_VCDF) + length(m_vec) * 2
))
# MVN probability estimation ------------------------------
# Estimate MVN probabilities for `niter` times using VMET, VMET with multi-level
# Monte Carlo, SOV, TLR SOV, VCDF, and MET
# Please refer to the paper https://arxiv.org/abs/2311.09426 for the details of
# each method
for (i in 1:niter) {
  cat("i = ", i, "\n")
  ## VMET ---------------------------
  est_Vecc <- rep(NA, length(m_vec))
  time_Vecc <- rep(NA, length(m_vec))
  for (j in 1:length(m_vec)) {
    m <- m_vec[j]
    time_Vecc[j] <- system.time(est_Vecc[j] <- VeccTMVN::pmvn(
      a, b, 0,
      sigma = covmat, reorder = 0, m = m,
      NLevel1 = 10, NLevel2 = 1e3
    ))[[3]]
  }
  cat("VeccTMVN done\n")
  ## VMET MLMC ---------------------------
  est_Vecc_MLMC <- rep(NA, length(m_vec))
  time_Vecc_MLMC <- rep(NA, length(m_vec))
  for (j in 1:length(m_vec)) {
    m <- m_vec[j]
    time_Vecc_MLMC[j] <- system.time(est_Vecc_MLMC[j] <- VeccTMVN::pmvn_MLMC(
      a, b, 0,
      sigma = covmat, reorder = 0, m1 = m, m2 = 2 * m,
      NLevel1 = 10, NLevel2 = 1e3
    ))[[3]]
  }
  cat("VeccTMVN MLMC done\n")
  ## TLR method -----------------------
  est_TLR <- rep(NA, length(N_SOV_TLR))
  time_TLR <- rep(NA, length(N_SOV_TLR))
  for (j in 1:length(N_SOV_TLR)) {
    cat("TLR", j, "\n")
    N <- N_SOV_TLR[j]
    err_obj <- try(
      time_TLR[j] <- system.time(
        est_TLR[j] <- tlrmvnmvt::pmvn(a, b,
          sigma = covmat,
          algorithm = tlrmvnmvt::TLRQMC(N = round(N / 20), m = sqrt(n), epsl = 1e-6)
        )
      )[[3]]
    )
    if (class(err_obj) == "try-error") {
      time_TLR[j] <- NA
      est_TLR[j] <- NA
    }
  }
  ## SOV method -----------------------
  est_SOV <- rep(NA, length(N_SOV_TLR))
  time_SOV <- rep(NA, length(N_SOV_TLR))
  for (j in 1:length(N_SOV_TLR)) {
    cat("SOV", j, "\n")
    N <- N_SOV_TLR[j]
    time_SOV[j] <- system.time(
      est_SOV[j] <- tlrmvnmvt::pmvn(a, b,
        sigma = covmat,
        algorithm = tlrmvnmvt::GenzBretz(N = round(N / 20))
      )
    )[[3]]
  }
  ## MET method -----------------------
  est_TN <- rep(NA, length(N_TN_VCDF))
  time_TN <- rep(NA, length(N_TN_VCDF))
  for (j in 1:length(N_TN_VCDF)) {
    cat("TN", j, "\n")
    N <- N_TN_VCDF[j]
    time_TN[j] <- system.time(est_TN[j] <- TruncatedNormal::pmvnorm(
      rep(0, n), covmat,
      lb = a, ub = b, B = N
    ))[[3]]
  }
  ## VCDF method -----------------------
  est_Nascimento <- rep(NA, length(N_TN_VCDF))
  time_Nascimento <- rep(NA, length(N_TN_VCDF))
  for (j in 1:length(N_TN_VCDF)) {
    cat("VCDF", j, "\n")
    N <- N_TN_VCDF[j]
    err_obj <- try(
      time_Nascimento[j] <- system.time(
        est_Nascimento[j] <- exp(CDFNormalAproxPackNoC.pmvn(b,
          mean = rep(0, n),
          sigma = covmat, p = round(N / 20)
        ))
      )[[3]]
    )
    if (class(err_obj) == "try-error") {
      time_Nascimento[j] <- NA
      est_Nascimento[j] <- NA
    }
  }
  ### save results ------------------------
  time_df[i, ] <- c(time_Vecc, time_Vecc_MLMC, time_TN, time_TLR, time_SOV, time_Nascimento)
  prob_df[i, ] <- c(est_Vecc, est_Vecc_MLMC, est_TN, est_TLR, est_SOV, est_Nascimento)
}
# assign colors to methods to align different figures
col_names <- c(
  paste("VMET", m_vec, sep = "_"),
  paste("VMET-ML", m_vec, sep = "_"),
  paste("MET", N_TN_VCDF, sep = "_"),
  paste("TLR", N_SOV_TLR, sep = "_"),
  paste("SOV", N_SOV_TLR, sep = "_"),
  paste("VCDF", N_TN_VCDF, sep = "_")
)
colors <- hue_pal()(6)
colors <- c(
  rep(colors[1], length(m_vec)),
  rep(colors[2], length(m_vec)),
  rep(colors[3], length(N_TN_VCDF)),
  rep(colors[4], length(N_SOV_TLR)),
  rep(colors[5], length(N_SOV_TLR)),
  rep(colors[6], length(N_TN_VCDF))
)
colnames(time_df) <- col_names
colnames(prob_df) <- col_names
if (!file.exists("results")) {
  dir.create("results")
}
save(m_vec, time_df, prob_df, col_names, colors, file = paste0(
  "results/Vecc_bias_const_corr_exp.RData"
))
# save the results. Once you have the results, you can comment out the 
# MVN probability estimation section to obtain the plots faster

# Plotting -----------------------------------
load(paste0(
  "results/Vecc_bias_const_corr_exp.RData"
))
library(ggplot2)
library(tidyr)
library(scales)
six_colors <- hue_pal()(6)
six_shapes <- c(1:6)
names(six_colors) <- c("VMET", "VMET-ML", "MET", "TLR", "SOV", "VCDF")
names(six_shapes) <- c("VMET", "VMET-ML", "MET", "TLR", "SOV", "VCDF")
time_vs_err_plt <- function(probDf, timeDf, colNames) {
  # colname_benchmark <- "MET"
  # benchmark <- mean(probDf[[colname_benchmark]])
  benchmark <- 1 / (n + 1)
  rmse <- sqrt(colMeans((probDf - benchmark)^2))
  secs <- colMeans(timeDf)
  mydf <- data.frame(
    time = secs, rmse = rmse,
    method = gsub("[_0-9]*", "", colNames)
  )
  mydf$method <- factor(mydf$method, levels = unique(mydf$method))
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
      name = "time (seconds)"
    ) +
    scale_color_manual(values = six_colors) +
    scale_shape_manual(values = six_shapes) +
    theme(
      text = element_text(size = 16),
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.position = c(0.75, 0.7)
    )
}
if (!file.exists("plots")) {
  dir.create("plots")
}
time_vs_err_plt(
  prob_df[, 1:(length(m_vec) * 2)],
  time_df[, 1:(length(m_vec) * 2)],
  col_names[1:(length(m_vec) * 2)]
)
ggsave(paste0("plots/err_vs_time_const_corr_VMET_only.pdf"),
  width = 6,
  height = 5
)
time_vs_err_plt(prob_df, time_df, col_names)
ggsave(paste0("plots/err_vs_time_const_corr.pdf"),
  width = 6,
  height = 5
)
