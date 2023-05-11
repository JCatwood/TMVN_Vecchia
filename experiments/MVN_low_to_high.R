rm(list = ls())

# Simulation -------------------------------
library(VeccTMVN)
library(tlrmvnmvt)
library(GpGp)
source("../funcs/utils.R")
## MVN prob gen funcs ----------------------------------------
prob_gen_fix_dom <- function(n, d, upper = 0) {
  grid_obj <- grid_gen(n, d)
  locs <- grid_obj$grid
  unit_dist <- grid_obj$unit_dist
  a <- rep(-Inf, n)
  b <- rep(upper, n)
  cov_parms <- c(1.0, 0.1, 0.03)
  cov_name <- "matern15_isotropic"
  cov_mat <- get(cov_name)(cov_parms, locs)
  smoothness <- 1.5
  return(list(
    a = a, b = b, locs = locs, cov_parms = cov_parms,
    cov_name = cov_name, smoothness = smoothness, dom_type = "fixed_grid"
  ))
}
prob_gen_exp_dom <- function(n, d, upper = 0) {
  grid_obj <- grid_gen(n, d)
  locs <- grid_obj$grid
  unit_dist <- grid_obj$unit_dist
  a <- rep(-Inf, n)
  b <- rep(upper, n)
  cov_parms <- c(1.0, 3 * unit_dist, 0.01)
  cov_name <- "matern15_isotropic"
  cov_mat <- get(cov_name)(cov_parms, locs)
  smoothness <- 1.5
  return(list(
    a = a, b = b, locs = locs, cov_parms = cov_parms,
    cov_name = cov_name, smoothness = smoothness, dom_type = "expand_grid"
  ))
}
## Prob setups -----------------------
set.seed(123)
n_vec <- c(900, 1600, 2500, 3600, 4900, 6400)
d <- 2
m <- 30
upper <- -2
## Iteratively compute the same MVN prob -----------------------
niter <- 30
time_df_Vecc <- data.frame(matrix(NA, niter, length(n_vec)))
err_df_Vecc <- data.frame(matrix(NA, niter, length(n_vec)))
prob_df_Vecc <- data.frame(matrix(NA, niter, length(n_vec)))
time_df_TLR <- data.frame(matrix(NA, niter, length(n_vec)))
err_df_TLR <- data.frame(matrix(NA, niter, length(n_vec)))
prob_df_TLR <- data.frame(matrix(NA, niter, length(n_vec)))
for (j in 1:length(n_vec)) {
  ### prob gen -----------------------
  n <- n_vec[j]
  prob_obj <- prob_gen_exp_dom(n, d, upper = upper)
  a <- prob_obj$a
  b <- prob_obj$b
  locs <- prob_obj$locs
  cov_name <- prob_obj$cov_name
  cov_parms <- prob_obj$cov_parms
  smoothness <- prob_obj$smoothness
  dom_type <- prob_obj$dom_type
  for (i in 1:niter) {
    ### Compute MVN prob with idea V -----------------------
    time_df_Vecc[i, j] <- system.time(est_Vecc <- VeccTMVN::pmvn(a, b, 0,
      locs = locs, covName = cov_name,
      reorder = 1, covParms = cov_parms,
      m = m, verbose = T,
      NLevel1 = 10, NLevel2 = 1e3
    ))[[3]]
    prob_df_Vecc[i, j] <- est_Vecc
    err_df_Vecc[i, j] <- attr(est_Vecc, "error")
    ### Compute MVN prob with tlrmvnmvt -----------------------
    cov_parms_TLR <- c(cov_parms[1], cov_parms[2], smoothness, cov_parms[3])
    time_df_TLR[i, j] <- system.time(
      est_TLR <- tlrmvnmvt::pmvn(a, b,
        geom = locs, kernelType = "matern", para = cov_parms_TLR,
        algorithm = tlrmvnmvt::GenzBretz(N = 500)
      )
    )[[3]]
    prob_df_TLR[i, j] <- est_TLR
    err_df_TLR[i, j] <- attr(est_TLR, "error")
  }
}

if (!file.exists("results")) {
  dir.create("results")
}
save(n_vec, time_df_Vecc, prob_df_Vecc, err_df_Vecc,
  time_df_TLR, prob_df_TLR, err_df_TLR,
  file = paste0("results/MVN_low_to_high_", dom_type, "_b", b[1], ".RData")
)

# Plotting -----------------------------------
load("results/MVN_low_to_high_expand_grid.RData")
library(ggplot2)
library(tidyr)
library(tibble)
box_plt_low_to_high_cmb <- function(..., yLim = NULL, yName = NULL,
                                    legends = NULL) {
  mtd_names <- c(paste0("n = ", n_vec))
  dfs <- list(...)
  for (i in 1:length(dfs)) {
    colnames(dfs[[i]]) <- mtd_names
    dfs[[i]] <- pivot_longer(dfs[[i]], cols = 1:ncol(dfs[[i]]), names_to = "method")
  }
  if (is.null(legends)) {
    legends <- as.character(c(1:length(dfs)))
  }
  df_cmb <- add_column(dfs[[1]], group = legends[1])
  if (length(dfs) > 1) {
    for (i in 2:length(dfs)) {
      df_cmb <- rbind(df_cmb, add_column(dfs[[i]], group = legends[i]))
    }
  }
  ggplot(df_cmb, aes(x = method, y = value, fill = group)) +
    scale_x_discrete(limits = mtd_names) +
    scale_y_continuous(limits = yLim, name = yName) +
    geom_boxplot()
}
box_plt_low_to_high_cmb(log(prob_df_Vecc), log(prob_df_TLR),
  yName = "Log-probability Estimates",
  legends = c("Vecc", "TLR")
)
box_plt_low_to_high_cmb(prob_df_Vecc, prob_df_TLR,
  yName = "Probability Estimates",
  legends = c("Vecc", "TLR")
)
box_plt_low_to_high_cmb(time_df_Vecc, time_df_TLR,
  yName = "Time", legends = c("Vecc", "TLR")
)
