rm(list = ls())

# Simulation -------------------------------
library(VeccTMVN)
library(GpGp)
source("../funcs/utils.R")
## MVN prob gen funcs ----------------------------------------
prob_gen_fix_dom <- function(n, d) {
  grid_obj <- grid_gen(n, d)
  locs <- grid_obj$grid
  unit_dist <- grid_obj$unit_dist
  a <- rep(-Inf, n)
  b <- rep(0, n)
  cov_parms <- c(1.0, 0.1, 0.03)
  cov_name <- "matern15_isotropic"
  cov_mat <- get(cov_name)(cov_parms, locs)
  smoothness <- 1.5
  return(list(
    a = a, b = b, locs = locs, cov_parms = cov_parms,
    cov_name = cov_name, smoothness = smoothness, dom_type = "fixed_grid"
  ))
}
prob_gen_exp_dom <- function(n, d) {
  grid_obj <- grid_gen(n, d)
  locs <- grid_obj$grid
  unit_dist <- grid_obj$unit_dist
  a <- rep(-Inf, n)
  b <- rep(0, n)
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
n <- 6400
d <- 2
m <- 30
N_level2_vec <- c(1:4)^2 * 1e3
prob_obj <- prob_gen_exp_dom(n, d)
a <- prob_obj$a
b <- prob_obj$b
locs <- prob_obj$locs
cov_name <- prob_obj$cov_name
cov_parms <- prob_obj$cov_parms
smoothness <- prob_obj$smoothness
dom_type <- prob_obj$dom_type
## Iteratively compute the same MVN prob -----------------------
niter <- 30
time_df_Vecc <- data.frame(matrix(NA, niter, length(N_level2_vec)))
err_df_Vecc <- data.frame(matrix(NA, niter, length(N_level2_vec)))
prob_df_Vecc <- data.frame(matrix(NA, niter, length(N_level2_vec)))
for (i in 1:niter) {
  for (j in 1:length(N_level2_vec)) {
    N_level2 <- N_level2_vec[j]
    ### Compute MVN prob with idea V -----------------------
    time_df_Vecc[i, j] <- system.time(est_Vecc <- VeccTMVN::pmvn(a, b, 0,
      locs = locs, covName = cov_name,
      reorder = 1, covParms = cov_parms,
      m = m, verbose = T,
      NLevel1 = 10, NLevel2 = N_level2
    ))[[3]]
    prob_df_Vecc[i, j] <- est_Vecc
    err_df_Vecc[i, j] <- attr(est_Vecc, "error")
  }
}

if (!file.exists("results")) {
  dir.create("results")
}
save(N_level2_vec, time_df_Vecc, prob_df_Vecc, err_df_Vecc,
  file = paste0("results/Vecc_MC_sample_", dom_type, ".RData")
)

# Plotting -----------------------------------
load("results/Vecc_MC_sample_expand_grid.RData")
library(ggplot2)
library(tidyr)
box_plt_MC_sample <- function(mydf, yLim = NULL, yName = NULL) {
  mtd_names <- c(paste0("NLevel2 = ", N_level2_vec))
  colnames(mydf) <- mtd_names
  mydf_pivot <- pivot_longer(mydf, cols = 1:ncol(mydf), names_to = "NLevel2")
  ggplot(mydf_pivot, aes(x = NLevel2, y = value)) +
    scale_x_discrete(limits = mtd_names) +
    scale_y_continuous(limits = yLim, name = yName) +
    geom_boxplot()
}
get_relerr <- function(df) {
  col_mean_mat <- matrix(colMeans(df),
    nrow = nrow(df),
    ncol = ncol(df), byrow = T
  )
  abs((df - col_mean_mat) / col_mean_mat)
}
prob_relerr_df_Vecc <- get_relerr(prob_df_Vecc)
box_plt_MC_sample(prob_df_Vecc, yName = "Vecc MVN prob")
box_plt_MC_sample(log(prob_df_Vecc), yName = "Vecc MVN log-prob")
box_plt_MC_sample(time_df_Vecc, yName = "Vecc time")
