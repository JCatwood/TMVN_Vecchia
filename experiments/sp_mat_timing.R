library(Matrix)
library(bench)
library(VeccTMVN)

rm(list = ls())

set.seed(1)
n <- 1e4
m <- 30
row_inds <- c(sapply(1:n, function(x) {
  rep(x, m)
}))
col_inds <- c(sapply(1:n, function(x) {
  sample(1:n, m)
}))
vals <- rnorm(n * m)
mat_sp <- sparseMatrix(
  i = row_inds,
  j = col_inds,
  x = vals
)
# matrix multiplication -----------------------------------------
mat_ds <- as.matrix(mat_sp)
bench::mark(
  mat_sp %*% t(mat_sp),
  mat_ds %*% t(mat_ds),
  check = F
)
# query coeffs in matrix self-conjugate multiplication -----------------------
nnz_inds <- cbind(row_inds, col_inds)[sample(1:length(row_inds), 2e3, F), ]
bench::mark(
  sum(sp_mat_mul_query(
    nnz_inds[, 1], nnz_inds[, 2], mat_sp@i, mat_sp@p,
    mat_sp@x
  )),
  sum(apply(nnz_inds, 1, function(ind) {
    sum(mat_sp[, ind[1]] * mat_sp[, ind[2]])
  })),
  sum(apply(nnz_inds, 1, function(ind) {
    sum(mat_ds[, ind[1]] * mat_ds[, ind[2]])
  })),
  check = T
)
# sparse triangular matrix solve -----------------------
diag(mat_ds) <- 0.1 + runif(n)
y <- runif(n)
mat_ds <- triu(mat_ds)
mat_sp <- as(mat_ds, "CsparseMatrix")
bench::mark(
  solve(mat_ds, y),
  solve(mat_sp, y),
  check = T
)
# sparse matrix vector multiplication -----------------------
bench::mark(
  as.vector(t(mat_sp) %*% y),
  as.vector(t(y) %*% mat_sp),
  check = T
)
# extract diag of sparse upper tri -----------------------
sp_diag_mat <- sparseMatrix(
  i = rep(1:(10 * n)),
  j = rep(1:(10 * n)),
  x = 1
)
bench::mark(
  diag(sp_diag_mat),
  sp_diag_mat@x[sp_diag_mat@p[-1]],
  check = T
)
