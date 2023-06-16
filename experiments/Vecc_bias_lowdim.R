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
## Prob setups -----------------------
set.seed(123)
n <- 900
d <- 2
m_vec <- seq(from = 10, to = 30, by = 10)
prob_ind <- 1
prob_obj <- get(paste0("prob", prob_ind, "_gen"))(n, d, retDenseCov = T)
a <- prob_obj$a
b <- prob_obj$b
locs <- prob_obj$locs
cov_mat <- prob_obj$cov_mat
cov_name <- prob_obj$cov_name
cov_parms <- prob_obj$cov_parms
## Iteratively compute the same MVN prob -----------------------
niter <- 30
time_df <- data.frame(matrix(NA, niter, 3 + length(m_vec)))
prob_df <- data.frame(matrix(NA, niter, 3 + length(m_vec)))
for (i in 1:niter) {
  est_Vecc <- rep(NA, length(m_vec))
  time_Vecc <- rep(NA, length(m_vec))
  for (j in 1:length(m_vec)) {
    ### Compute MVN prob with idea V -----------------------
    m <- m_vec[j]
    time_Vecc[j] <- system.time(est_Vecc[j] <- VeccTMVN::pmvn(a, b, 0,
      locs = locs, covName = cov_name,
      reorder = 1, covParms = cov_parms,
      m = m, verbose = T,
      NLevel1 = 10, NLevel2 = 1e3,
    ))[[3]]
  }
  ### Compute MVN prob with other methods -----------------------
  time_TN <- system.time(est_TN <- TruncatedNormal::pmvnorm(
    rep(0, n), cov_mat,
    lb = a, ub = b, B = 1e4
  ))[[3]]
  time_TLR <- system.time(
    est_TLR <- tlrmvnmvt::pmvn(a, b,
      sigma = cov_mat,
      algorithm = tlrmvnmvt::GenzBretz(N = 500)
    )
  )[[3]]
  time_Nascimento <- system.time(
    est_Nascimento <- exp(CDFNormalAproxPackNoC.pmvn(b,
      mean = rep(0, n),
      sigma = cov_mat
    ))
  )[[3]]
  ### save results ------------------------
  time_df[i, ] <- c(time_Vecc, time_TN, time_TLR, time_Nascimento)
  prob_df[i, ] <- c(est_Vecc, est_TN, est_TLR, est_Nascimento)
}
if (!file.exists("results")) {
  dir.create("results")
}
save(m_vec, time_df, prob_df, file = paste0(
  "results/Vecc_bias_lowdim_exp",
  prob_ind, ".RData"
))

# Plotting -----------------------------------
load(paste0(
  "results/Vecc_bias_lowdim_exp",
  prob_ind, ".RData"
))
library(ggplot2)
library(tidyr)
box_plt_low_dim <- function(mydf, yLim = NULL, yName = NULL,
                            yTrans = "identity") {
  mtd_names <- c(paste0("m = ", m_vec), "ET", "TLR", "VCDF")
  colnames(mydf) <- mtd_names
  mydf_pivot <- pivot_longer(mydf, cols = 1:ncol(mydf), names_to = "method")
  ggplot(mydf_pivot, aes(x = method, y = value)) +
    scale_x_discrete(limits = mtd_names, name = NULL) +
    scale_y_continuous(limits = yLim, name = yName, trans = yTrans) +
    geom_boxplot()
}
box_plt_low_dim_Botev_only <- function(mydf, yLim = NULL, yName = NULL,
                                       yTrans = "identity") {
  mtd_names <- c(paste0("m = ", m_vec), "ET")
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
box_plt_low_dim(prob_df, yName = "log MVN prob", yTrans = "log2")
ggsave(paste0("plots/logprob_lowdim_exp", prob_ind, ".pdf"),
  width = 5,
  height = 5
)
box_plt_low_dim_Botev_only(prob_df[, 1:(length(m_vec) + 1)],
  yName = "MVN prob"
)
ggsave(paste0("plots/prob_lowdim_exp", prob_ind, ".pdf"), width = 5, height = 5)
box_plt_low_dim(time_df, yName = "time (seconds)", yTrans = "log2")
ggsave(paste0("plots/time_lowdim_exp", prob_ind, ".pdf"), width = 5, height = 5)
