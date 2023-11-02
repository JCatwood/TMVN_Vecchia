library(GpGp)
source("../funcs/inv_chol.R")
## example MVN probabilities --------------------------------
n1 <- 10
n2 <- 10
n <- n1 * n2
locs <- as.matrix(expand.grid((1:n1) / n1, (1:n2) / n2))
covparms <- c(2, 0.3, 0)
cov_mat <- matern15_isotropic(covparms, locs)

## random diag matrix with 1 and -1
n_rnd_D <- 10

## ordering --------------------------------
ord <- order_maxmin(locs)
locs_ord <- locs[ord, , drop = FALSE]
cov_mat_ord <- matern15_isotropic(covparms, locs_ord)

## KL function
KL_GP_zero_mean <- function(cov_1, inv_chol_2) {
  n <- nrow(cov_1)
  0.5 * (-2 * sum(log(diag(inv_chol_2))) -
    determinant(cov_1, logarithm = T)$modulus[1] -
    n + sum(diag(U %*% t(U) %*% cov_1)))
}

## find KL at different m --------------------------------
m_list <- round(seq(from = 1, to = 30, length.out = 5))
m_max <- max(m_list)
NNarray_max <- find_ordered_nn(locs_ord, m = m_max)
KL <- matrix(0, n_rnd_D, length(m_list))
for (k in 1:n_rnd_D) {
  D <- sample(c(-1, 1), n, replace = T)
  ### compute DKD --------------------------------
  D_cov_mat_ord_D <- D * t(t(cov_mat_ord) * D)
  for (j in 1:length(m_list)) {
    m <- m_list[j]
    NNarray <- NNarray_max[, 1:(m + 1)]
    #### Vecchia approx --------------------------------
    U <- get_sp_inv_chol(D_cov_mat_ord_D, NNarray)
    KL[k, j] <- KL_GP_zero_mean(D_cov_mat_ord_D, U)
  }
}
