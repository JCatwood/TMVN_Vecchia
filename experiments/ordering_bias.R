rm(list = ls())

# Simulation -------------------------------
library(VeccTMVN)
library(TruncatedNormalBeta)
library(TruncatedNormal)
library(GpGp)
source("../funcs/utils.R")
## MVN prob gen funcs ----------------------------------------
prob1_gen <- function(n, d, retDenseCov = F, reorder = F) {
  locs <- grid_gen(n, d)$grid
  a <- rep(-Inf, n)
  b <- rep(0, n)
  cov_parms <- c(1.0, 0.1, 0.01)
  cov_name <- "matern15_isotropic"
  cov_mat <- get(cov_name)(cov_parms, locs)
  if (reorder) {
    odr <- TruncatedNormal::cholperm(cov_mat, a, b)$perm
    a <- a[odr]
    b <- b[odr]
    locs <- locs[odr, , drop = F]
    cov_mat <- get(cov_name)(cov_parms, locs)
  }
  return(list(
    a = a, b = b, locs = locs, cov_parms = cov_parms,
    cov_name = cov_name, cov_mat = cov_mat
  ))
}
prob2_gen <- function(n, d, retDenseCov = F, reorder = F) {
  locs <- latin_gen(n, d)
  a <- rep(-Inf, n)
  b <- -runif(n, 0, 2)
  cov_parms <- c(1.0, 0.1, 0.01)
  cov_name <- "matern15_isotropic"
  cov_mat <- get(cov_name)(cov_parms, locs)
  if (reorder) {
    odr <- TruncatedNormal::cholperm(cov_mat, a, b)$perm
    a <- a[odr]
    b <- b[odr]
    locs <- locs[odr, , drop = F]
    cov_mat <- get(cov_name)(cov_parms, locs)
  }
  return(list(
    a = a, b = b, locs = locs, cov_parms = cov_parms,
    cov_name = cov_name, cov_mat = cov_mat
  ))
}
prob3_gen <- function(n, d, retDenseCov = F, reorder = F) {
  locs <- grid_gen(n, d)$grid
  a <- rep(-1, n)
  b <- rep(1, n)
  cov_parms <- c(1.0, 0.1, 0.01)
  cov_name <- "matern15_isotropic"
  cov_mat <- get(cov_name)(cov_parms, locs)
  if (reorder) {
    odr <- TruncatedNormal::cholperm(cov_mat, a, b)$perm
    a <- a[odr]
    b <- b[odr]
    locs <- locs[odr, , drop = F]
    cov_mat <- get(cov_name)(cov_parms, locs)
  }
  return(list(
    a = a, b = b, locs = locs, cov_parms = cov_parms,
    cov_name = cov_name, cov_mat = cov_mat
  ))
}
## Prob setups -----------------------
set.seed(123)
n <- 900
d <- 2
m <- 30
prob_ind <- 1
prob_obj <- get(paste0("prob", prob_ind, "_gen"))(n, d, retDenseCov = T,
  reorder = F)
prob_obj_reorder <- get(paste0("prob", prob_ind, "_gen"))(n, d, retDenseCov = T,
  reorder = T)
## Iteratively compute the same MVN prob -----------------------
niter <- 30
time_df <- data.frame(matrix(NA, niter, 4))
prob_df <- data.frame(matrix(NA, niter, 4))
for (i in 1:niter) {
  ### Compute MVN prob with idea V -----------------------
  time_Vecc <- system.time(est_Vecc <- VeccTMVN::pmvn(prob_obj$a, prob_obj$b, 0,
    locs = prob_obj$locs, covName = prob_obj$cov_name,
    reorder = 0, covParms = prob_obj$cov_parms,
    m = m, verbose = T,
    NLevel1 = 10, NLevel2 = 1e3
  ))[[3]]
  time_Vecc_reorder <- system.time(
    est_Vecc_reorder <- VeccTMVN::pmvn(
      prob_obj_reorder$a, prob_obj_reorder$b, 0,
      locs = prob_obj_reorder$locs, covName = prob_obj_reorder$cov_name,
      reorder = 0, covParms = prob_obj_reorder$cov_parms,
      m = m, verbose = T,
      NLevel1 = 10, NLevel2 = 1e3
    )
  )[[3]]
  ### Compute MVN prob with other methods -----------------------
  time_TN <- system.time(est_TN <- TruncatedNormalBeta::pmvnorm(
    rep(0, n), prob_obj$cov_mat,
    lb = prob_obj$a, ub = prob_obj$b, B = 1e4
  ))[[3]]
  time_TN_reorder <- system.time(est_TN_reorder <- TruncatedNormalBeta::pmvnorm(
    rep(0, n), prob_obj_reorder$cov_mat,
    lb = prob_obj_reorder$a, ub = prob_obj_reorder$b, B = 1e4
  ))[[3]]
  ### save results ------------------------
  time_df[i, ] <- c(time_Vecc, time_Vecc_reorder, time_TN, time_TN_reorder)
  prob_df[i, ] <- c(est_Vecc, est_Vecc_reorder, est_TN, est_TN_reorder)
}
if (!file.exists("results")) {
  dir.create("results")
}
save(m_vec, time_df, prob_df, file = paste0(
  "results/ordering_bias",
  prob_ind, ".RData"
))

# Plotting -----------------------------------
load(paste0(
  "results/ordering_bias",
  prob_ind, ".RData"
))
library(ggplot2)
library(tidyr)
box_plt_low_dim <- function(mydf, yLim = NULL, yName = NULL,
                            yTrans = "identity") {
  mtd_names <- c("Vecc", "Vecc_reorder", "ET", "ET_reorder")
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
box_plt_low_dim(prob_df, yName = "MVN prob")
ggsave(paste0("plots/ordering_bias", prob_ind, ".pdf"),
  width = 5,
  height = 5
)
