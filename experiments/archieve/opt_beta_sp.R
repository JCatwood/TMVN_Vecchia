rm(list = ls())
library(GpGp)
library(nleqslv)
library(truncnorm)
library(VeccTMVN)

## example MVN probabilities --------------------------------
set.seed(123)
n1 <- 10
n2 <- 10
n <- n1 * n2
locs <- as.matrix(expand.grid((1:n1) / n1, (1:n2) / n2))
covparms <- c(2, 0.3, 0)
cov_mat <- matern15_isotropic(covparms, locs)
a_list <- list(rep(-4, n), rep(-1, n), -runif(n) * 2)
b_list <- list(-runif(n) * 2, rep(1, n), runif(n) * 2)

## ordering and NN --------------------------------
m <- 10
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
vecc_cond_mean_var_obj_sp <- vecc_cond_mean_var_sp(NNarray, cov_mat_ord)

## Compare system solution --------------------------------
for (i in 1:length(a_list_ord)) {
  a_ord <- a_list_ord[[i]]
  b_ord <- b_list_ord[[i]]
  trunc_expect <- etruncnorm(a_ord, b_ord)
  x0 <- c(trunc_expect[-n], rep(0, n - 1))
  x0_all <- c(trunc_expect, rep(0, n))
  cat("Dense nleqslv used", system.time(
    solv_idea_5 <- nleqslv(x0,
      fn = grad_idea5,
      jac = jac_idea5,
      veccCondMeanVarObj = vecc_cond_mean_var_obj,
      a = a_ord, b = b_ord,
      global = "pwldog",
      method = "Newton",
      control = list(maxit = 500L)
    )
  )[[3]], "seconds\n")
  cat("Sparse nleqslv used", system.time(
    solv_idea_5_sp <- my_nleqslv(
      x0_all,
      fn = function(x, ...) {
        ret <- grad_jacprod_jacsolv_idea5(x, ...,
          retJac = F,
          retProd = F, retSolv = F
        )
        ret$grad
      },
      jacTransFn = function(x, ...) {
        ret <- grad_jacprod_jacsolv_idea5(x, ...,
          retJac = F,
          retProd = T, retSolv = F
        )
        ret$jac_grad
      },
      jacInvFn = function(x, ...) {
        ret <- grad_jacprod_jacsolv_idea5(x, ...,
          retJac = F,
          retProd = F, retSolv = T
        )
        ret$jac_inv_grad
      },
      veccCondMeanVarObj = vecc_cond_mean_var_obj_sp,
      a = a_ord, b = b_ord, verbose = F,
      control = list(maxit = 1000L)
    )
  )[[3]], "seconds\n")
  cat("Sparse optim used", system.time(
    solv_idea_5_sp_optim <- optim(
      x0_all,
      fn = function(x, ...) {
        ret <- grad_jacprod_jacsolv_idea5(x, ...,
          retJac = F,
          retProd = F, retSolv = F
        )
        0.5 * sum((ret$grad)^2)
      },
      gr = function(x, ...) {
        ret <- grad_jacprod_jacsolv_idea5(x, ...,
          retJac = F,
          retProd = T, retSolv = F
        )
        ret$jac_grad
      },
      method = "L-BFGS-B",
      veccCondMeanVarObj = vecc_cond_mean_var_obj_sp,
      a = a_ord, b = b_ord, verbose = F,
      lower = c(a_ord, rep(-Inf, n)), upper = c(b_ord, rep(Inf, n)),
      control = list(maxit = 1000L)
    )
  )[[3]], "seconds\n")


  ### Check solution consistency -------------------------------
  cat("Terminal codes of dense optimization is", solv_idea_5$termcd, "\n")
  cat("Terminal codes of sparse optimization is", solv_idea_5_sp$code, "\n")
  cat(
    "Difference dense nleqslv and sparse nleqslv is",
    sum(abs(solv_idea_5$x - solv_idea_5_sp$x[-c(n, 2 * n)])), "\n"
  )
  cat(
    "Difference dense nleqslv and sparse optim is",
    sum(abs(solv_idea_5$x - solv_idea_5_sp_optim$par[-c(n, 2 * n)])), "\n"
  )
}
