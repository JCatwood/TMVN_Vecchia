library(GpGp)
library(truncnorm)

#' Univariate ordering under FIC approximation
#'
#' @example
#' n1 <- 5
#' n2 <- 5
#' n <- n1 * n2
#' m <- 5
#' locs <- as.matrix(expand.grid((1:n1) / n1, (1:n2) / n2))
#' covparms <- c(2, 0.1, 0)
#' cov_name <- "matern15_isotropic"
#' a <- rep(-Inf, n)
#' b <- seq(from = -3, to = 3, length.out = n)
#' cat("The output order should be roughly increasing after", m, "numbers \n")
#' cat(FIC_univar_reorder(a, b, locs, cov_name, covparms, m))
FIC_univar_reorder <- function(a, b, locs, covName, covParms, m) {
  n <- nrow(locs)
  m <- min(m, n - 1)
  ## standardize ---------------------
  sd_par <- sqrt(covParms[1])
  a <- a / sd_par
  b <- b / sd_par
  covParms[1] <- 1
  ## maxmin order --------------------------
  odr_maxmin <- GpGp::order_maxmin(locs)
  odr_FIC_univar <- odr_maxmin
  a <- a[odr_maxmin]
  b <- b[odr_maxmin]
  locs <- locs[odr_maxmin, , drop = F]
  ## cov func -------------------------------
  cov_func_GpGp <- getFromNamespace(covName, "GpGp")
  cov_func <- function(ind) {
    cov_func_GpGp(covParms, locs[ind, , drop = F])
  }
  ## TMVN expectation ------------------------------
  x_first_m <- truncnorm::etruncnorm(a[1:m], b[1:m])
  ## compute TMVN probs under FIC ----------------------------
  tmvn_prob_1D <- rep(-1.0, n)
  for (i in (m + 1):n) {
    ind_cond <- c(i, 1:m)
    cov_mat_sub <- cov_func(ind_cond)
    cov_mat_sub_inv <- solve(cov_mat_sub[-1, -1])
    cov_vec_sub <- cov_mat_sub[1, -1]
    cond_sd <- sqrt(cov_mat_sub[1, 1] - as.numeric(
      t(cov_vec_sub) %*% cov_mat_sub_inv %*% cov_vec_sub
    ))
    cond_mean <- as.numeric(t(cov_vec_sub) %*% cov_mat_sub_inv %*% x_first_m)
    tmvn_prob_1D[i] <- pnorm(b[i], mean = cond_mean, sd = cond_sd) -
      pnorm(a[i], mean = cond_mean, sd = cond_sd)
  }
  odr_maxmin[order(tmvn_prob_1D, decreasing = F)]
}
