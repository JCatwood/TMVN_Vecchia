rm(list = ls())
library(GpGp)
library(nleqslv)
source("../funcs/nleqsv.R")
source("../funcs/gradpsi_sp.R")
source("../funcs/gradpsi.R")
source("../funcs/vecc_cond_mean_var.R")

## example MVN probabilities --------------------------------
set.seed(123)
n1 <- 10
n2 <- 10
n <- n1 * n2
locs <- as.matrix(expand.grid((1:n1) / n1, (1:n2) / n2))
covparms <- c(2, 0.3, 0)
cov_mat <- matern15_isotropic(covparms, locs)
a_list <- list(rep(-Inf, n), rep(-1, n), -runif(n) * 2)
b_list <- list(rep(-2, n), rep(1, n), runif(n) * 2)

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
vecc_cond_mean_var_obj <- vecc_cond_mean_var(cov_mat_ord, NNarray)
vecc_cond_mean_var_obj_sp <- vecc_cond_mean_var_sp(cov_mat_ord, NNarray)

## Compare system solution --------------------------------
for (i in 3:length(a_list_ord)) {
  a_ord <- a_list_ord[[i]]
  b_ord <- b_list_ord[[i]]
  x0 <- rep(0, 2 * length(a_ord) - 2)
  solv_idea_5 <- nleqslv(x0,
    fn = grad_idea5,
    jac = jac_idea5,
    veccCondMeanVarObj = vecc_cond_mean_var_obj,
    a = a_ord, b = b_ord,
    global = "pwldog",
    method = "Newton",
    control = list(maxit = 500L)
  )
  x0 <- rep(0, 2 * length(a_ord))
  solv_idea_5_sp <- my_nleqslv(x0,
    fn = grad_idea5_sp,
    gradNewtonFn = function(x, ...) {
      ret <- grad_jacprod_jacsolv_idea5(x, ...)
      list(
        grad = ret$jac_grad,
        Newton_step = -ret$jac_inv_grad
      )
    },
    veccCondMeanVarObj = vecc_cond_mean_var_obj_sp,
    a = a_ord, b = b_ord,
    control = list(maxit = 500L)
  )
  ### Check solution consistency -------------------------------
  cat("Terminal codes of dense optimization is", solv_idea_5$termcd, "\n")
  cat("Terminal codes of sparse optimization is", solv_idea_5_sp$code, "\n")
  cat(
    "Difference dense optimization and sparse optimization is",
    sum(abs(solv_idea_5$x - solv_idea_5_sp$x[-c(n, 2 * n)])), "\n"
  )
}
