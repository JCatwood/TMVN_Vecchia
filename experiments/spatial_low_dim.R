rm(list = ls())

# Simulation -------------------------------
library(VeccTMVN)
library(tlrmvnmvt)
library(TruncatedNormal)
library(mvtnorm)
library(lhs)
library(GpGp)
source("../funcs/CDFNormalAproxPackNoC_pmvn.R")
## Locs gen funcs ----------------------------------------
latin_gen <- function(n, d) {
  as.matrix(lhs::randomLHS(n, d))
}
grid_gen <- function(n, d, unitDist = NULL) {
  m <- floor(n^{
    1 / d
  }) + 1
  if (is.null(unitDist)) {
    unitDist <- 1 / m
  }
  edges <- list()
  for (i in 1:d) {
    edges[[i]] <- (1:m) * unitDist
  }
  grid <- as.matrix(expand.grid(edges))
  return(list(grid = grid[1:n, , drop = F], unit_dist = unitDist))
}
## MVN prob gen funcs ----------------------------------------
prob1_gen <- function(n, d, retDenseCov = F) {
  locs <- latin_gen(n, d)
  a <- rep(-Inf, n)
  b <- rep(-2, n)
  cov_parms <- c(1.0, 0.1, 0.0)
  cov_name <- "matern15_isotropic"
  if (retDenseCov) {
    cov_mat <- get(cov_name)(cov_parms, locs)
    return(list(
      cov_mat = cov_mat, a = a, b = b, locs = locs,
      cov_parms = cov_parms, cov_name = cov_name
    ))
  } else {
    return(list(a = a, b = b, locs = locs, cov_parms = cov_parms, cov_name = cov_name))
  }
}
prob2_gen <- function(n, d, retDenseCov = F) {
  locs <- grid_gen(n, d)$grid
  a <- rep(-Inf, n)
  b <- rep(0, n)
  cov_parms <- c(1.0, 0.1, 0.0)
  cov_name <- "matern15_isotropic"
  if (retDenseCov) {
    cov_mat <- get(cov_name)(cov_parms, locs)
    return(list(
      a = a, b = b, locs = locs, cov_parms = cov_parms,
      cov_name = cov_name, cov_mat = cov_mat
    ))
  } else {
    return(list(
      a = a, b = b, locs = locs, cov_parms = cov_parms,
      cov_name = cov_name
    ))
  }
}
## Prob setups -----------------------
set.seed(123)
n <- 400
d <- 2
m <- 30
prob_obj <- prob2_gen(n, d, retDenseCov = T)
a <- prob_obj$a
b <- prob_obj$b
locs <- prob_obj$locs
cov_mat <- prob_obj$cov_mat
cov_name <- prob_obj$cov_name
cov_parms <- prob_obj$cov_parms
## Iteratively compute the same MVN prob -----------------------
niter <- 30
time_df <- data.frame(matrix(NA, niter, 4))
prob_df <- data.frame(matrix(NA, niter, 4))
for (i in 1:niter) {
  ### Compute MVN prob with idea V -----------------------
  time_Vecc <- system.time(est_Vecc <- VeccTMVN::pmvn(a, b, 0,
    locs = locs, covName = cov_name,
    reorder = 1, covParms = cov_parms,
    m = m, verbose = T,
    NLevel1 = 10, NLevel2 = 1e3,
  ))[[3]]
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
save(time_df, prob_df, file = "results/low_dim_prob_time.RData")

# Plotting -----------------------------------
load("results/low_dim_prob_time.RData")
library(ggplot2)
library(tidyr)
box_plt_low_dim <- function(mydf) {
  colnames(mydf) <- c("VeccTMVN", "Botev", "TLRMVN", "Nascimento")
  mydf_pivot <- pivot_longer(mydf, cols = 1:4, names_to = "method")
  ggplot(mydf_pivot, aes(x = method, y = value)) +
    geom_boxplot()
}
box_plt_low_dim(prob_df)
box_plt_low_dim(time_df)
