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
  ### ordering based on integration limits --------------------------------
  pnorm_at_a <- pnorm(a_list[[i]], sd = sqrt(covparms[1]))
  pnorm_at_b <- pnorm(b_list[[i]], sd = sqrt(covparms[1]))
  ord_int_lim <- order(pnorm_at_b - pnorm_at_a, decreasing = F)
  ### reorder ----------------------------------
  locs_ord_int_lim <- locs[ord_int_lim, , drop = F]
  cov_mat_ord_int_lim <- matern15_isotropic(covparms, locs_ord_int_lim)
  a_ord_int_lim <- a_list[[i]][ord_int_lim]
  b_ord_int_lim <- b_list[[i]][ord_int_lim]
  ### maxmin ordering ----------------------------------
  ord_maxmin <- GpGp::order_maxmin(locs)
  ### reorder ----------------------------------
  locs_ord_maxmin <- locs[ord_maxmin, , drop = F]
  cov_mat_ord <- matern15_isotropic(covparms, locs_ord_maxmin)
  a_ord_maxmin <- a_list[[i]][ord_maxmin]
  b_ord_maxmin <- b_list[[i]][ord_maxmin]
  ### compare ---------------------
  est_TN_reorder <- TruncatedNormal::pmvnorm(
    rep(0, n), cov_mat,
    lb = a_list[[i]], ub = b_list[[i]]
  )
  est_Vecc_maxmin <- VeccTMVN::pmvn(a_ord_maxmin, b_ord_maxmin, 0,
    locs = locs_ord_maxmin,
    covName = "matern15_isotropic",
    covParms = covparms, m = 30
  )
  est_Vecc_int_lim <- VeccTMVN::pmvn(a_ord_int_lim, b_ord_int_lim, 0,
    locs = locs_ord_int_lim,
    covName = "matern15_isotropic",
    covParms = covparms, m = 30
  )
  cat(
    "TN with reorder est", est_TN_reorder,
    "rel err", attributes(est_TN_reorder)$relerr, "\n"
  )
  cat(
    "VeccTMVN est maxmin", est_Vecc_maxmin,
    "rel err", attributes(est_Vecc_maxmin)$error / est_Vecc_maxmin, "\n"
  )
  cat(
    "VeccTMVN est int lim", est_Vecc_int_lim,
    "rel err", attributes(est_Vecc_int_lim)$error / est_Vecc_int_lim, "\n"
  )
}
