library(Matrix)
library(bench)

rm(list = ls())

set.seed(1)
n <- 1e4
m <- 30
row_inds <- c(sapply(1 : n, function(x){rep(x, m)}))
col_inds <- c(sapply(1 : n, function(x){sample(1 : n, m)}))
vals <- rnorm(n * m)
mat_sp <- sparseMatrix(
  i = row_inds,
  j = col_inds,
  x = vals
)

mat_ds <- as.matrix(mat_sp)
bench::mark(
  mat_sp %*% t(mat_sp),
  mat_ds %*% t(mat_ds),
  check = F
)