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
## local pmvn func for testing --------------------------------
pmvn <- function(a, b, orderTest) {
  a_ord <- a[orderTest]
  b_ord <- b[orderTest]
  locs_ord <- locs[orderTest, , drop = F]
  VeccTMVN::pmvn(a_ord, b_ord, 0,
    locs = locs_ord,
    covName = "matern15_isotropic",
    covParms = covparms, m = 30
  )
}
## Compute MVN probs --------------------------------
for (i in 1:length(a_list)) {
  ### ordering based on integration limits --------------------------------
  pnorm_at_a <- pnorm(a_list[[i]], sd = sqrt(covparms[1]))
  pnorm_at_b <- pnorm(b_list[[i]], sd = sqrt(covparms[1]))
  ord_int_lim <- order(pnorm_at_b - pnorm_at_a, decreasing = F)
  est_Vecc_int_lim <- pmvn(a_list[[i]], b_list[[i]], ord_int_lim)
  ### maxmin ordering ----------------------------------
  ord_maxmin <- GpGp::order_maxmin(locs)
  est_Vecc_maxmin <- pmvn(a_list[[i]], b_list[[i]], ord_maxmin)
  ### univariate under FIC ------------------------------
  ord_FIC_univar <- VeccTMVN::FIC_univar_reorder(a_list[[i]], b_list[[i]],
    m = 30, locs,
    "matern15_isotropic", covparms,
  )
  est_Vecc_FIC_univar <- pmvn(a_list[[i]], b_list[[i]], ord_FIC_univar)
  ### compare ---------------------
  est_TN_reorder <- TruncatedNormal::pmvnorm(
    rep(0, n), cov_mat,
    lb = a_list[[i]], ub = b_list[[i]]
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
  cat(
    "VeccTMVN est FIC_univar", est_Vecc_FIC_univar,
    "rel err", attributes(est_Vecc_FIC_univar)$error / est_Vecc_FIC_univar, "\n"
  )
}
