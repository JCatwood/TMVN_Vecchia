rm(list = ls())

# Simulation -------------------------------
library(VeccTMVN)
library(tlrmvnmvt)
library(TruncatedNormal)
library(mvtnorm)
library(GpGp)
source("../funcs/CDFNormalAproxPackNoC_pmvn.R")
## MVN prob gen funcs ----------------------------------------
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
m_vec <- seq(from = 10, to = 100, by = 10)
prob_obj <- prob_gen(n)
a <- prob_obj$a
b <- prob_obj$b
covmat <- prob_obj$covmat
## Iteratively compute the same MVN prob -----------------------
niter <- 30
time_df <- data.frame(matrix(NA, niter, 4 + length(m_vec)))
prob_df <- data.frame(matrix(NA, niter, 4 + length(m_vec)))
# for (i in 1:niter) {
#   cat("i = ", i, "\n")
#   est_Vecc <- rep(NA, length(m_vec))
#   time_Vecc <- rep(NA, length(m_vec))
#   for (j in 1:length(m_vec)) {
#     ### Compute MVN prob with idea V -----------------------
#     m <- m_vec[j]
#     time_Vecc[j] <- system.time(est_Vecc[j] <- VeccTMVN::pmvn(
#       a, b, 0, sigma = covmat, reorder = 0, m = m, 
#       NLevel1 = 10, NLevel2 = 1e3
#     ))[[3]]
#   }
#   cat("VeccTMVN done\n")
#   ### Compute MVN prob with other methods -----------------------
#   err_obj <- try(
#     time_TLR <- system.time(
#       est_TLR <- tlrmvnmvt::pmvn(a, b,
#         sigma = covmat,
#         algorithm = tlrmvnmvt::TLRQMC(N = 500, m = sqrt(n), epsl = 1e-6)
#       )
#     )[[3]]
#   )
#   if (class(err_obj) == "try-error") {
#     time_TLR <- NA
#     est_TLR <- NA
#   }
#   cat("TLRMVN done\n")
#   time_SOV <- system.time(
#     est_SOV <- tlrmvnmvt::pmvn(a, b,
#       sigma = covmat,
#       algorithm = tlrmvnmvt::GenzBretz(N = 500)
#     )
#   )[[3]]
#   cat("SOV done\n")
#   time_TN <- system.time(est_TN <- TruncatedNormal::pmvnorm(
#     rep(0, n), covmat,
#     lb = a, ub = b, B = 1e4
#   ))[[3]]
#   cat("TN done\n")
#   time_Nascimento <- system.time(
#     est_Nascimento <- exp(CDFNormalAproxPackNoC.pmvn(b,
#       mean = rep(0, n),
#       sigma = covmat
#     ))
#   )[[3]]
#   cat("VCDF done\n")
#   ### save results ------------------------
#   time_df[i, ] <- c(time_Vecc, time_TN, time_TLR, time_SOV, time_Nascimento)
#   prob_df[i, ] <- c(est_Vecc, est_TN, est_TLR, est_SOV, est_Nascimento)
# }
# if (!file.exists("results")) {
#   dir.create("results")
# }
# save(m_vec, time_df, prob_df, file = paste0(
#   "results/Vecc_bias_const_corr_exp.RData"
# ))

# Plotting -----------------------------------
load(paste0(
  "results/Vecc_bias_const_corr_exp.RData"
))
library(ggplot2)
library(tidyr)
library(scales)
box_plt_low_dim <- function(mydf, yName = NULL,
                            yTrans = "identity") {
  mtd_names <- c(paste0("m", m_vec), "MET", "TLR", "SOV", "VCDF")
  i <- which(mtd_names == "MET")
  tmp <- mydf[, i]
  mydf[, i] <- mydf[, length(mtd_names)]
  mydf[, length(mtd_names)] <- tmp
  tmp <- mtd_names[i]
  mtd_names[i] <- mtd_names[length(mtd_names)]
  mtd_names[length(mtd_names)] <- tmp
  colnames(mydf) <- mtd_names
  mydf_pivot <- pivot_longer(mydf, cols = 1:ncol(mydf), names_to = "method")
  hex <- hue_pal()(5)
  if (yTrans == "identity") {
    breaks <- seq(from = min(mydf_pivot$value), to = max(mydf_pivot$value), length.out = 5)
  } else {
    breaks <- exp(seq(
      from = log(min(mydf_pivot$value)),
      to = log(max(mydf_pivot$value)), length.out = 5
    ))
  }

  ggplot(mydf_pivot, aes(x = method, y = value)) +
    geom_boxplot(fill = c(rep(hex[1], length(m_vec)), hex[2:5])) +
    scale_x_discrete(limits = mtd_names) +
    scale_y_continuous(
      name = yName, trans = yTrans, breaks = breaks,
      labels = signif(breaks, digits = 2)
    ) +
    ggtitle(paste("Constant correlation ", covmat[1, 2])) +
    theme(
      text = element_text(size = 14), legend.position = "none",
      axis.title.x = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
}
box_plt_low_dim_Botev_only <- function(mydf, yLim = NULL, yName = NULL,
                                       yTrans = "identity") {
  mtd_names <- c(paste0("m = ", m_vec), "MET")
  colnames(mydf) <- mtd_names
  mydf_pivot <- pivot_longer(mydf, cols = 1:ncol(mydf), names_to = "method")
  ggplot(mydf_pivot, aes(x = method, y = value)) +
    scale_x_discrete(limits = mtd_names, name = NULL) +
    scale_y_continuous(limits = yLim, name = yName, trans = yTrans) +
    geom_boxplot()
}
if (!file.exists("plots")) {
  dir.create("plots")
}
# box_plt_low_dim(prob_df, yName = "log-probability estimates", yTrans = "log2")
box_plt_low_dim(prob_df, yName = "probability estimates")
ggsave(paste0("plots/Vecc_bias_const_corr_exp.pdf"),
  width = 5,
  height = 5
)
# box_plt_low_dim_Botev_only(prob_df[, 1:(length(m_vec) + 1)],
#   yName = "MVN prob"
# )
# ggsave(paste0("plots/prob_const_corr_exp.pdf"), width = 5, height = 5)
box_plt_low_dim(time_df, yName = "Time (seconds)", yTrans = "log2")
ggsave(paste0("plots/time_const_corr_exp.pdf"), width = 5, height = 5)
