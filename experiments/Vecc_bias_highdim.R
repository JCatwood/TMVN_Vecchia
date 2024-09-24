rm(list = ls())

# Simulation -------------------------------
library(VeccTMVN)
library(tlrmvnmvt)
library(GpGp)
source("../funcs/CDFNormalAproxPackNoC_pmvn.R")
source("../funcs/utils.R")
## MVN prob gen funcs ----------------------------------------
prob1_gen <- function(n, d, ...) {
  locs <- grid_gen(n, d)$grid
  a <- rep(-Inf, n)
  b <- rep(0, n)
  cov_parms <- c(1.0, 0.1, 0.03)
  cov_name <- "matern15_isotropic"
  cov_mat <- get(cov_name)(cov_parms, locs)
  return(list(
    a = a, b = b, locs = locs, cov_parms = cov_parms,
    cov_name = cov_name, cov_mat = cov_mat
  ))
}
prob2_gen <- function(n, d, ...) {
  locs <- latin_gen(n, d)
  a <- rep(-Inf, n)
  b <- -runif(n, 0, 2)
  cov_parms <- c(1.0, 0.1, 0.03)
  cov_name <- "matern15_isotropic"
  cov_mat <- get(cov_name)(cov_parms, locs)
  return(list(
    a = a, b = b, locs = locs, cov_parms = cov_parms,
    cov_name = cov_name, cov_mat = cov_mat
  ))
}
prob3_gen <- function(n, d, ...) {
  locs <- grid_gen(n, d)$grid
  a <- rep(-1, n)
  b <- rep(1, n)
  cov_parms <- c(1.0, 0.1, 0.03)
  cov_name <- "matern15_isotropic"
  cov_mat <- get(cov_name)(cov_parms, locs)
  return(list(
    a = a, b = b, locs = locs, cov_parms = cov_parms,
    cov_name = cov_name, cov_mat = cov_mat
  ))
}
## Prob setups -----------------------
myargs <- commandArgs(trailingOnly = TRUE)
set.seed(321)
n <- 6400
d <- 2
m_vec <- seq(from = 30, to = 50, by = 10)
N_SOV <- c(5e3, 1e4, 2e4)
N_TLR <- c(2e4, 5e4, 10e4)
if (length(myargs) > 0) {
  prob_ind <- as.numeric(myargs[1])  
} else {
  prob_ind <- 1  
}
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
time_df <- data.frame(matrix(NA, niter, length(N_TLR) + length(N_SOV) + length(m_vec)))
prob_df <- data.frame(matrix(NA, niter, length(N_TLR) + length(N_SOV) + length(m_vec)))
for (i in 1:niter) {
  est_Vecc <- rep(NA, length(m_vec))
  time_Vecc <- rep(NA, length(m_vec))
  for (j in 1:length(m_vec)) {
    cat("VeccTMVN", j, "\n")
    ### Compute MVN prob with idea V -----------------------
    m <- m_vec[j]
    time_Vecc[j] <- system.time(est_Vecc[j] <- VeccTMVN::pmvn(a, b, 0,
      locs = locs, covName = cov_name,
      reorder = 3, covParms = cov_parms,
      m = m, verbose = T,
      NLevel1 = 10, NLevel2 = 1e4
    ))[[3]]
  }
  ### Compute MVN prob with other methods -----------------------
  est_TLR <- rep(NA, length(N_TLR))
  time_TLR <- rep(NA, length(N_TLR))
  for (j in 1:length(N_TLR)) {
    cat("TLR", j, "\n")
    N <- N_TLR[j]
    err_obj <- try(
      time_TLR[j] <- system.time(
        est_TLR[j] <- tlrmvnmvt::pmvn(a[z_order], b[z_order],
          sigma = cov_mat[z_order, z_order],
          algorithm = tlrmvnmvt::TLRQMC(N = round(N / 20), m = sqrt(n), epsl = 1e-6)
        )
      )[[3]]
    )
    if (class(err_obj) == "try-error") {
      time_TLR[j] <- NA
      est_TLR[j] <- NA
    }
  }

  est_SOV <- rep(NA, length(N_SOV))
  time_SOV <- rep(NA, length(N_SOV))
  for (j in 1:length(N_SOV)) {
    cat("SOV", j, "\n")
    N <- N_SOV[j]
    time_SOV[j] <- system.time(
      est_SOV[j] <- tlrmvnmvt::pmvn(a, b,
        sigma = cov_mat,
        algorithm = tlrmvnmvt::GenzBretz(N = round(N / 20))
      )
    )[[3]]
  }
  ### save results ------------------------
  time_df[i, ] <- c(time_Vecc, time_TLR, time_SOV)
  prob_df[i, ] <- c(est_Vecc, est_TLR, est_SOV)
}
if (!file.exists("results")) {
  dir.create("results")
}
save(m_vec, N_SOV, N_TLR, time_df, prob_df, file = paste0(
  "results/Vecc_bias_highdim_exp",
  prob_ind, ".RData"
))

# Plotting -----------------------------------
load(paste0(
  "results/Vecc_bias_highdim_exp",
  prob_ind, ".RData"
))
library(ggplot2)
library(tidyr)
library(scales)
box_plt_low_dim <- function(mydf, yName = NULL,
                            yTrans = "identity") {
  mtd_names <- c(paste0("m", m_vec), "TLR", "SOV")
  colnames(mydf) <- mtd_names
  mydf_pivot <- pivot_longer(mydf, cols = 1:ncol(mydf), names_to = "method")
  hex <- hue_pal()(3)
  if (yTrans == "identity") {
    breaks <- seq(from = min(mydf_pivot$value), to = max(mydf_pivot$value), length.out = 5)
  } else {
    breaks <- exp(seq(
      from = log(min(mydf_pivot$value)),
      to = log(max(mydf_pivot$value)), length.out = 5
    ))
  }

  ggplot(mydf_pivot, aes(x = method, y = value)) +
    geom_boxplot(fill = c(rep(hex[1], length(m_vec)), hex[2:3])) +
    scale_x_discrete(limits = mtd_names) +
    scale_y_continuous(
      name = yName, trans = yTrans, breaks = breaks,
      labels = signif(breaks, digits = 2)
    ) +
    ggtitle(paste("Scenario", prob_ind)) +
    theme(
      text = element_text(size = 14), legend.position = "none",
      axis.title.x = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
}
if (!file.exists("plots")) {
  dir.create("plots")
}
box_plt_low_dim(prob_df, yName = "Log-probability estimates", yTrans = "log2")
ggsave(paste0("plots/logprob_highdim_exp", prob_ind, ".pdf"),
  width = 5,
  height = 5
)
tmp_plt <- box_plt_low_dim(time_df, yName = "Time (seconds)", yTrans = "log2") +
  ggtitle("Computation times")
tmp_plt
ggsave(paste0("plots/time_highdim_exp", prob_ind, ".pdf"),
  width = 5, height = 5
)
