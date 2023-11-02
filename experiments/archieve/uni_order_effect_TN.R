rm(list = ls())
library(TruncatedNormal)
library(VeccTMVN)
library(GpGp)

## example MVN probabilities --------------------------------
n1 <- 20
n2 <- 20
n <- n1 * n2
locs <- as.matrix(expand.grid((1:n1) / n1, (1:n2) / n2))
covparms <- c(2, 0.1, 0)
cov_mat <- matern15_isotropic(covparms, locs)
a_list <- list(rep(-Inf, n), rep(-Inf, n), rep(-Inf, n))
b_list <- list(rep(-2, n), rep(1, n), -runif(n) * 2)

## Compute MVN probs --------------------------------
for (i in 1:length(a_list)) {
  ### maxmin ordering ----------------------------------
  ord <- GpGp::order_maxmin(locs)
  ### reorder ----------------------------------
  locs_ord <- locs[ord, , drop = F]
  cov_mat_ord <- matern15_isotropic(covparms, locs_ord)
  a_ord <- a_list[[i]][ord]
  b_ord <- b_list[[i]][ord]
  ### compare ---------------------
  est_TN_no_reorder <- TruncatedNormalBeta::pmvnorm(
    rep(0, n), cov_mat_ord,
    lb = a_ord, ub = b_ord
  )
  est_TN_reorder <- TruncatedNormal::pmvnorm(
    rep(0, n), cov_mat_ord,
    lb = a_ord, ub = b_ord
  )
  est_Vecc <- VeccTMVN::pmvn(a_ord, b_ord, 0,
    locs = locs_ord,
    covName = "matern15_isotropic",
    covParms = covparms, m = 30
  )
  cat(
    "TN without reorder est", est_TN_no_reorder,
    "rel err", attributes(est_TN_no_reorder)$relerr, "\n"
  )
  cat(
    "TN with reorder est", est_TN_reorder,
    "rel err", attributes(est_TN_reorder)$relerr, "\n"
  )
  cat(
    "VeccTMVN est", est_Vecc,
    "rel err", attributes(est_Vecc)$error / est_Vecc, "\n"
  )
}
