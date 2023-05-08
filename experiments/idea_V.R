library(GpGp)
library(TruncatedNormal)
library(VeccTMVN)
## example MVN probabilities --------------------------------
n1 <- 10
n2 <- 10
n <- n1 * n2
locs <- as.matrix(expand.grid((1:n1) / n1, (1:n2) / n2))
covparms <- c(2, 0.3, 0)
cov_mat <- matern15_isotropic(covparms, locs)
a_list <- list(rep(-Inf, n), rep(-1, n), -runif(n) * 2)
b_list <- list(rep(-2, n), rep(1, n), runif(n) * 2)
# mu <- rep(0, n)
# a_list <- lapply(a_list, function(x){x - mu})
# b_list <- lapply(b_list, function(x){x - mu})

## ordering and NN --------------------------------
m <- 30
ord <- order_maxmin(locs)
locs_ord <- locs[ord, , drop = FALSE]
cov_mat_ord <- matern15_isotropic(covparms, locs_ord)
a_list_ord <- lapply(a_list, function(x) {
  x[ord]
})
b_list_ord <- lapply(b_list, function(x) {
  x[ord]
})
NNarray <- find_ordered_nn(locs_ord, m = m)

## Vecchia approx --------------------------------
U <- get_sp_inv_chol(cov_mat_ord, NNarray)
cov_mat_Vecc <- solve(U %*% t(U))
sqrt(mean((cov_mat_ord - cov_mat_Vecc)^2))

## Compare MVN prob errors caused by Vecchia approx ---------------
cat(
  "Comparing ", length(a_list_ord),
  " MVN probabilities computed with/out Vecchia\n"
)
for (i in 1:length(a_list_ord)) {
  a_ord <- a_list_ord[[i]]
  b_ord <- b_list_ord[[i]]
  tmp_rslt <-
    pmvnorm(mu = rep(0, n), sigma = cov_mat_ord, lb = a_ord, ub = b_ord)
  cat(
    "Without Vecchia: Est ", tmp_rslt[1], "Relerr ",
    attributes(tmp_rslt)$relerr, "\n"
  )
  tmp_rslt <-
    pmvnorm(mu = rep(0, n), sigma = cov_mat_Vecc, lb = a_ord, ub = b_ord)
  cat(
    "With Vecchia: Est ", tmp_rslt[1], "Relerr ",
    attributes(tmp_rslt)$relerr, "\n"
  )
}
cat("Done")
