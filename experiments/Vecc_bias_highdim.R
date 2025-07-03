rm(list = ls())

# Simulation -------------------------------
library(VeccTMVN)
library(tlrmvnmvt)
library(GpGp)
source("../funcs/CDFNormalAproxPackNoC_pmvn.R")
source("../funcs/utils.R")
## MVN prob gen funcs ----------------------------------------
#' Generate MVN probabilities of Scenarios 1/2/3 in Table 1 of \href{https://arxiv.org/abs/2311.09426}{the paper}
#' @param n number of locations, also the dimension of the MVN probability
#' @param d dimension of the spatial domain, where the Matern kernel is defined
#' @param ... not used, for aligning interface
#' @return a list containing:
#' \\item{a}{lower integration limits}
#' \\item{b}{upper integration limits}
#' \\item{locs}{spatial locations for the spatial covariance matrix}
#' \\item{cov_parms}{covariance parameters corresponding to `cov_name` in the `GpGp` package}
#' \\item{cov_name}{covariance name as defined in `GpGp` package}
#' \\item{cov_mat}{the generate covariance matrix for the MVN probability}
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
set.seed(123)
n <- 6400 # dimension of the MVN probability
d <- 2 # dimension of the spatial domain where the covariance matrix is defined
m_vec <- c(30, 50, 70) # conditioning set sizes to try for VMET method
N_SOV <- c(2e4, 5e4, 10e4) # MC sample size to try for SOV method
N_TLR <- c(2e4, 5e4, 10e4) # MC sample size to try for TLR method
# optional, one can use command line arguments to specify the scenario index
if (length(myargs) > 0) {
  prob_ind <- as.numeric(myargs[1])
} else {
  prob_ind <- 3 # default scenario ID
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
# MVN probability estimation ------------------------------
# Estimate MVN probabilities for `niter` times using VMET, SOV, and TLR SOV
# Please refer to the paper https://arxiv.org/abs/2311.09426 for the details of
# each method
for (i in 1:niter) {
  ## VMET ---------------------------
  est_Vecc <- rep(NA, length(m_vec))
  time_Vecc <- rep(NA, length(m_vec))
  for (j in 1:length(m_vec)) {
    cat("VeccTMVN", j, "\n")
    m <- m_vec[j]
    time_Vecc[j] <- system.time(est_Vecc[j] <- VeccTMVN::pmvn(
      a, b, 0,
      locs = locs, covName = cov_name,
      reorder = 3, covParms = cov_parms,
      m = m, verbose = T,
      NLevel1 = 10, NLevel2 = 1e4
    ))[[3]]
  }
  ## TLR SOV method -----------------------
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
  ## SOV method -----------------------
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
# save the results. Once you have the results, you can comment out the 
# MVN probability estimation section to obtain the plots faster

# Plotting -----------------------------------
load(paste0(
  "results/Vecc_bias_highdim_exp",
  prob_ind, ".RData"
))
library(ggplot2)
library(tidyr)
library(scales)
library(dplyr)
six_colors <- hue_pal()(6)
six_shapes <- c(1:6)
names(six_colors) <- c("VMET", "VMET-ML", "MET", "TLR", "SOV", "VCDF")
names(six_shapes) <- c("VMET", "VMET-ML", "MET", "TLR", "SOV", "VCDF")
time_vs_err_plt <- function(probDf, timeDf) {
  col_names <- c(
    paste("VMET", m_vec, sep = "_"),
    paste("TLR", N_TLR, sep = "_"),
    paste("SOV", N_SOV, sep = "_")
  )
  colnames(probDf) <- col_names
  colnames(timeDf) <- col_names
  if (prob_ind < 3) {
    colname_benchmark <- paste("VMET", max(m_vec), sep = "_")
  } else {
    colname_benchmark <- paste("SOV", max(N_SOV), sep = "_")
  }
  benchmark <- log(mean(probDf[[colname_benchmark]]))
  rmse <- sqrt(colMeans((log(probDf) - benchmark)^2))
  secs <- colMeans(timeDf)
  mydf <- data.frame(
    time = secs, rmse = rmse,
    method = gsub("[_0-9e+]*", "", col_names)
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
    scale_color_manual(values = six_colors) +
    scale_shape_manual(values = six_shapes) +
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
time_vs_err_plt(prob_df, time_df)
ggsave(
  paste0("plots/err_vs_time_highdim_exp", prob_ind, "_new", ".pdf"),
  width = 5,
  height = 5
)
