rm(list = ls())

# Simulation -------------------------------
library(VeccTMVN)
library(TruncatedNormalBeta)
library(GpGp)
library(R.utils)
source("../funcs/utils.R")
## MVN prob gen funcs ----------------------------------------
prob1_gen <- function(n, d, m) {
  locs <- grid_gen(n, d)$grid
  a <- rep(-Inf, n)
  b <- rep(0, n)
  cov_parms <- c(1.0, 0.1, 0.01)
  cov_name <- "matern15_isotropic"
  cov_mat <- get(cov_name)(cov_parms, locs)
  odr_FIC <- FIC_reorder_univar(a, b, m, covMat = cov_mat)
  odr_FIC100 <- FIC_reorder_univar(a, b, 100, covMat = cov_mat)
  odr_Vecc <- Vecc_reorder(a, b, m, covMat = cov_mat)$order
  odr_univar <- TruncatedNormalBeta::cholperm(cov_mat, a, b)$perm
  return(list(
    a = a, b = b, locs = locs, cov_parms = cov_parms,
    cov_name = cov_name, cov_mat = cov_mat,
    odr_FIC = odr_FIC, odr_FIC100 = odr_FIC100, odr_Vecc = odr_Vecc,
    odr_univar = odr_univar
  ))
}
prob2_gen <- function(n, d, m) {
  locs <- latin_gen(n, d)
  a <- rep(-Inf, n)
  b <- -runif(n, 0, 2)
  cov_parms <- c(1.0, 0.1, 0.01)
  cov_name <- "matern15_isotropic"
  cov_mat <- get(cov_name)(cov_parms, locs)
  odr_FIC <- FIC_reorder_univar(a, b, m, covMat = cov_mat)
  odr_FIC100 <- FIC_reorder_univar(a, b, 100, covMat = cov_mat)
  odr_Vecc <- Vecc_reorder(a, b, m, covMat = cov_mat)$order
  odr_univar <- TruncatedNormalBeta::cholperm(cov_mat, a, b)$perm
  return(list(
    a = a, b = b, locs = locs, cov_parms = cov_parms,
    cov_name = cov_name, cov_mat = cov_mat,
    odr_FIC = odr_FIC, odr_FIC100 = odr_FIC100, odr_Vecc = odr_Vecc,
    odr_univar = odr_univar
  ))
}
prob3_gen <- function(n, d, m) {
  locs <- grid_gen(n, d)$grid
  a <- rep(-1, n)
  b <- rep(1, n)
  cov_parms <- c(1.0, 0.1, 0.01)
  cov_name <- "matern15_isotropic"
  cov_mat <- get(cov_name)(cov_parms, locs)
  odr_FIC <- FIC_reorder_univar(a, b, m, covMat = cov_mat)
  odr_FIC100 <- FIC_reorder_univar(a, b, 100, covMat = cov_mat)
  odr_Vecc <- Vecc_reorder(a, b, m, covMat = cov_mat)$order
  odr_univar <- TruncatedNormalBeta::cholperm(cov_mat, a, b)$perm
  return(list(
    a = a, b = b, locs = locs, cov_parms = cov_parms,
    cov_name = cov_name, cov_mat = cov_mat,
    odr_FIC = odr_FIC, odr_FIC100 = odr_FIC100, odr_Vecc = odr_Vecc,
    odr_univar = odr_univar
  ))
}
## Experiment function -------------------------
exp_func <- function(pron_obj, odr, runTN = T) {
  time_Vecc <- system.time(est_Vecc <- VeccTMVN::pmvn(
    prob_obj$a[odr], prob_obj$b[odr], 0,
    locs = prob_obj$locs[odr, , drop = F], covName = prob_obj$cov_name,
    reorder = 0, covParms = prob_obj$cov_parms,
    m = m, verbose = T,
    NLevel1 = 10, NLevel2 = 1e3
  ))[[3]]
  if (runTN) {
    time_TN <- system.time(est_TN <- withTimeout(TruncatedNormalBeta::pmvnorm(
      rep(0, n), prob_obj$cov_mat[odr, odr],
      lb = prob_obj$a[odr], ub = prob_obj$b[odr], B = 1e4
    ), timeout = 60, onTimeout = "silent"))[[3]]
    if (is.null(est_TN)) {
      est_TN <- NA
      time_TN <- NA
    }
    return(list(
      est_Vecc = est_Vecc, time_Vecc = time_Vecc,
      est_TN = est_TN, time_TN = time_TN
    ))
  } else {
    return(list(
      est_Vecc = est_Vecc, time_Vecc = time_Vecc
    ))
  }
}
## Prob setups -----------------------
n <- 900
d <- 2
m <- 30
## Iteratively compute the same MVN prob -----------------------
nprob <- 3
niter <- 30
nmtd <- 2
nodr <- 5
rslt <- data.frame(matrix(NA, nprob * niter * nmtd * nodr, 6))
colnames(rslt) <- c("prob_ind", "iter_ind", "mtd", "order", "est", "time")
for (prob_ind in c(1:nprob)) {
  set.seed(123)
  prob_obj <- get(paste0("prob", prob_ind, "_gen"))(n, d, m)
  for (i in 1:niter) {
    ind_offset <- (prob_ind - 1) * niter * nodr * nmtd +
      (i - 1) * nodr * nmtd
    rslt_noodr <- exp_func(prob_obj, 1:n, T)
    rslt_FIC <- exp_func(prob_obj, prob_obj$odr_FIC, F)
    rslt_FIC100 <- exp_func(prob_obj, prob_obj$odr_FIC100, F)
    rslt_Vecc <- exp_func(prob_obj, prob_obj$odr_Vecc, F)
    rslt_univar <- exp_func(prob_obj, prob_obj$odr_univar, T)
    rslt[ind_offset + 1, ] <- c(
      prob_ind, i, "VMET", "no_order",
      rslt_noodr$est_Vecc, rslt_noodr$time_Vecc
    )
    rslt[ind_offset + 2, ] <- c(
      prob_ind, i, "VMET", "FIC",
      rslt_FIC$est_Vecc, rslt_FIC$time_Vecc
    )
    rslt[ind_offset + 3, ] <- c(
      prob_ind, i, "VMET", "FIC100",
      rslt_FIC100$est_Vecc, rslt_FIC100$time_Vecc
    )
    rslt[ind_offset + 4, ] <- c(
      prob_ind, i, "VMET", "Vecc",
      rslt_Vecc$est_Vecc, rslt_Vecc$time_Vecc
    )
    rslt[ind_offset + 5, ] <- c(
      prob_ind, i, "VMET", "univar",
      rslt_univar$est_Vecc, rslt_univar$time_Vecc
    )
    rslt[ind_offset + 6, ] <- c(
      prob_ind, i, "MET", "no_order",
      rslt_noodr$est_TN, rslt_noodr$time_TN
    )
    rslt[ind_offset + 7, ] <- c(prob_ind, i, "MET", "FIC", NA, NA)
    rslt[ind_offset + 8, ] <- c(prob_ind, i, "MET", "FIC100", NA, NA)
    rslt[ind_offset + 9, ] <- c(prob_ind, i, "MET", "Vecc", NA, NA)
    rslt[ind_offset + 10, ] <- c(
      prob_ind, i, "MET", "univar",
      rslt_univar$est_TN, rslt_univar$time_TN
    )
  }
}

if (!file.exists("results")) {
  dir.create("results")
}
save(rslt, file = paste0(
  "results/ordering_bias.RData"
))

# Plotting -----------------------------------
load(paste0(
  "results/ordering_bias.RData"
))
library(ggplot2)
library(scales)
library(tidyr)
library(dplyr)
if (!file.exists("plots")) {
  dir.create("plots")
}
rslt <- as_tibble(rslt)
rslt$est <- as.numeric(rslt$est)
rslt$time <- as.numeric(rslt$time)
for (i in 1:nprob) {
  rslt_tmp <- rslt %>%
    filter(order %in% c("no_order", "univar")) %>%
    filter(prob_ind == i) %>%
    filter(!is.na(est)) %>%
    unite("mtd_order", mtd, order, remove = F)
  rslt_tmp %>%
    ggplot(mapping = aes(x = mtd_order, y = est, fill = mtd)) +
    geom_boxplot() +
    ylab("Estimates of MVN probabilities") +
    scale_y_continuous(
      labels = signif(
        seq(
          from = min(rslt_tmp$est), to = max(rslt_tmp$est), length.out = 5
        ),
        digits = 2
      )
    ) +
    scale_x_discrete(labels = c(
      "MET", "MET with reorder",
      "VMET", "VMET with reorder"
    )) +
    ggtitle(paste("Scenario", i)) +
    theme(
      text = element_text(size = 14), legend.position = "none",
      axis.title.x = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  ggsave(paste0("plots/ordering_bias_MET_VMET", i, ".pdf"),
    width = 5,
    height = 5
  )
}

for (i in 1:nprob) {
  rslt_tmp <- rslt %>%
    filter(mtd == "VMET") %>%
    filter(prob_ind == i)
  rslt_tmp %>%
    ggplot(mapping = aes(x = order, y = est, fill = order)) +
    geom_boxplot() +
    ylab("Estimates of MVN probabilities") +
    scale_y_continuous(
      breaks = seq(
        from = min(rslt_tmp$est), to = max(rslt_tmp$est), length.out = 5
      ),
      labels = signif(
        seq(
          from = min(rslt_tmp$est), to = max(rslt_tmp$est), length.out = 5
        ),
        digits = 2
      )
    ) +
    scale_x_discrete(
      breaks = unique(rslt_tmp$order),
      labels = c(
        "No reorder", "FIC", "FIC100",
        "Vecchia", "Univariate"
      )
    ) +
    ggtitle(paste("Scenario", i)) +
    theme(
      text = element_text(size = 14), legend.position = "none",
      axis.title.x = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  ggsave(paste0("plots/ordering_bias_FIC_Vecc", i, ".pdf"),
    width = 5,
    height = 5
  )
}
