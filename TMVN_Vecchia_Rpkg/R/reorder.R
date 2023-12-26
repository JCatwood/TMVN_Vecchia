library(GpGp)
library(truncnorm)

# Univariate ordering under FIC approximation, first m chosen by maxmin ordering
#
# @param a lower bound vector for TMVN
# @param b upper bound vector for TMVN
# @param m Vecchia conditioning set size
# @param locs location (feature) matrix n X d
# @param covName covariance function name from the `GpGp` package
# @param covParms parameters for `covName`
# @param covMat dense covariance matrix, not needed when `locs` is not null
# @return a vector of new order based on FIC assumption and maxmin ordering
#
# example of a legacy function
# n1 <- 5
# n2 <- 5
# n <- n1 * n2
# m <- 5
# locs <- as.matrix(expand.grid((1:n1) / n1, (1:n2) / n2))
# covparms <- c(2, 0.1, 0)
# cov_name <- "matern15_isotropic"
# a <- rep(-Inf, n)
# b <- seq(from = -3, to = 3, length.out = n)
# cat("The output order should be roughly increasing after", m, "numbers")
# cat(FIC_reorder_maxmin(a, b, m, locs, cov_name, covparms))
FIC_reorder_maxmin <- function(a, b, m, locs = NULL, covName = NULL,
                               covParms = NULL, covMat = NULL) {
  if (!is.null(covMat)) {
    n <- nrow(covMat)
  } else {
    n <- nrow(locs)
  }
  m <- min(m, n - 1)
  ## standardize ---------------------
  if (!is.null(covMat)) {
    sd_par <- sqrt(diag(covMat))
    a <- a / sd_par
    b <- b / sd_par
    covMat <- t(t(covMat / sd_par) / sd_par)
  } else {
    sd_par <- sqrt(covParms[1])
    a <- a / sd_par
    b <- b / sd_par
    covParms[1] <- 1
  }
  ## maxmin order --------------------------
  if (!is.null(covMat)) {
    odr_maxmin <- sample(1:n, n, F)
  } else {
    odr_maxmin <- GpGp::order_maxmin(locs)
  }
  odr_FIC_univar <- odr_maxmin
  a <- a[odr_maxmin]
  b <- b[odr_maxmin]
  locs <- locs[odr_maxmin, , drop = F]
  ## cov func -------------------------------
  if (!is.null(covMat)) {
    cov_func <- function(ind) {
      covMat[ind, ind, drop = F]
    }
  } else {
    cov_func_GpGp <- utils::getFromNamespace(covName, "GpGp")
    cov_func <- function(ind) {
      cov_func_GpGp(covParms, locs[ind, , drop = F])
    }
  }
  ## TMVN expectation ------------------------------
  x_first_m <- truncnorm::etruncnorm(a[1:m], b[1:m])
  x_first_m[is.nan(x_first_m)] <- a[1:m]
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
    tmvn_prob_1D[i] <- stats::pnorm(b[i], mean = cond_mean, sd = cond_sd) -
      stats::pnorm(a[i], mean = cond_mean, sd = cond_sd)
  }
  odr_maxmin[order(tmvn_prob_1D, decreasing = F)]
}


#' Univariate ordering under FIC approximation, first m chosen by m iter of
#'   dense univariate reordering
#'
#' @param a lower bound vector for TMVN
#' @param b upper bound vector for TMVN
#' @param m Vecchia conditioning set size
#' @param locs location (feature) matrix n X d
#' @param covName covariance function name from the `GpGp` package
#' @param covParms parameters for `covName`
#' @param covMat dense covariance matrix, not needed when `locs` is not null
#' @return a vector of new order based on FIC assumption and maxmin ordering
#' 
#' @examples
#' library(VeccTMVN)
#' n1 <- 5
#' n2 <- 5
#' n <- n1 * n2
#' m <- 5
#' locs <- as.matrix(expand.grid((1:n1) / n1, (1:n2) / n2))
#' covparms <- c(2, 0.1, 0)
#' cov_name <- "matern15_isotropic"
#' a <- rep(-Inf, n)
#' b <- seq(from = -3, to = 3, length.out = n)
#' cat("The output order should be roughly 1 to ", n)
#' cat(FIC_reorder_univar(a, b, m, locs, cov_name, covparms))
#' 
#' @export
FIC_reorder_univar <- function(a, b, m, locs = NULL, covName = NULL,
                               covParms = NULL, covMat = NULL) {
  if (!is.null(covMat)) {
    n <- nrow(covMat)
  } else {
    n <- nrow(locs)
  }
  m <- min(m, n - 1)
  if (m > 100) {
    warning("Univariate reordering complexity grows at m^4.\n")
  }
  ## standardize ---------------------
  if (!is.null(covMat)) {
    sd_par <- sqrt(diag(covMat))
    a <- a / sd_par
    b <- b / sd_par
    covMat <- t(t(covMat / sd_par) / sd_par)
  } else {
    sd_par <- sqrt(covParms[1])
    a <- a / sd_par
    b <- b / sd_par
    covParms[1] <- 1
  }
  ## cov func -------------------------------
  if (is.null(covMat)) {
    cov_func_GpGp <- utils::getFromNamespace(covName, "GpGp")
    covMat <- cov_func_GpGp(covParms, locs)
  }
  cov_func <- function(indRow, indCol) {
    covMat[indRow, indCol, drop = F]
  }
  ## TMVN expectation ------------------------------
  x_first_m <- rep(0, m)
  ## Initial order ---------------------------
  odr <- 1:n
  ## m iterations of univar reordering ---------------------
  for (i in 1:m) {
    if (i > 1) {
      cov_mat_rows <- cov_func(odr[1:(i - 1)], odr)
      cov_mat_sub <- cov_mat_rows[, 1:(i - 1), drop = F]
      cov_mat_rows_sub <- cov_mat_rows[, i:n, drop = F]
      cov_mat_sub_inv <- solve(cov_mat_sub)
      mu_cond <- as.numeric(t(cov_mat_rows_sub) %*%
        (cov_mat_sub_inv %*% x_first_m[1:(i - 1)]))
      sd_cond <- sqrt(1 - apply(cov_mat_rows_sub, 2, function(x) {
        as.numeric(t(x) %*% cov_mat_sub_inv %*% x)
      }))
    } else {
      mu_cond <- rep(0, n + 1 - i)
      sd_cond <- rep(1, n + 1 - i)
    }
    a_tilde <- (a[i:n] - mu_cond) / sd_cond
    b_tilde <- (b[i:n] - mu_cond) / sd_cond
    tmvn_prob_1D <- stats::pnorm(b_tilde) - stats::pnorm(a_tilde)
    j <- which.min(tmvn_prob_1D)
    j_hat <- j + i - 1
    x_first_m[i] <- truncnorm::etruncnorm(a_tilde[j], b_tilde[j]) * sd_cond[j] +
      mu_cond[j]
    if (is.nan(x_first_m[i])) {
      x_first_m[i] <- a[j_hat]
    }
    tmp <- a[j_hat]
    a[j_hat] <- a[i]
    a[i] <- tmp
    tmp <- b[j_hat]
    b[j_hat] <- b[i]
    b[i] <- tmp
    tmp <- odr[j_hat]
    odr[j_hat] <- odr[i]
    odr[i] <- tmp
  }
  odr
}


#' Univariate ordering under Vecchia approximation
#'
#' @param a lower bound vector for TMVN
#' @param b upper bound vector for TMVN
#' @param m Vecchia conditioning set size
#' @param locs location (feature) matrix n X d
#' @param covName covariance function name from the `GpGp` package
#' @param covParms parameters for `covName`
#' @param covMat dense covariance matrix, not needed when `locs` is not null
#' @return a vector of new order based on FIC assumption and maxmin ordering
#'
#' @examples
#' library(lhs)
#' library(GpGp)
#' library(VeccTMVN)
#' set.seed(123)
#' n <- 100
#' m <- 5
#' locs <- lhs::geneticLHS(n, 2)
#' covparms <- c(1, 0.1, 0)
#' cov_name <- "matern15_isotropic"
#' cov_mat <- get(cov_name)(covparms, locs)
#' a <- rep(-Inf, n)
#' b <- runif(n)
#' odr_TN <- TruncatedNormal::cholperm(cov_mat, a, b)$perm
#' rslt <- Vecc_reorder(a, b, m,
#'   locs = locs, covName = cov_name,
#'   covParms = covparms
#' )
#' # compare order
#' cat(rslt$order)
#' cat(odr_TN)
# # check NN array
# locs_odr <- locs[rslt$order, , drop = FALSE]
# cov_mat_odr <- get(cov_name)(covparms, locs_odr)
# NN <- GpGp::find_ordered_nn(locs_odr, m)
# cat("nn difference is", sum(rslt$nn - NN, na.rm = TRUE))
# # check A matrix and conditional var
# rslt_check <- vecc_cond_mean_var_sp(rslt$nn, covMat = cov_mat_odr)
# cat("A difference is", sum(rslt$A - rslt_check$A))
# cat(
#   "Conditional variance difference is",
#   sum(rslt$cond_var - rslt_check$cond_var)
# )
#
#' @export
Vecc_reorder <- function(a, b, m, locs = NULL, covName = NULL,
                         covParms = NULL, covMat = NULL) {
  n <- length(a)
  if (is.null(covMat)) {
    covMat <- get(covName)(covParms, locs)
  }
  margin_sd <- sqrt(diag(covMat))
  b <- b / margin_sd
  a <- a / margin_sd
  corr_mat <- t(t(covMat / margin_sd) / margin_sd)
  rslt <- univar_order_vecc(a, b, corr_mat, m)
  # adjust for indexing
  rslt$order <- as.numeric(rslt$order + 1)
  rslt$nn <- cbind(1:n, rslt$nn + 1)
  for (i in 1:m) {
    rslt$nn[i, (i + 1):(m + 1)] <- NA
  }
  # create sparse matrix A
  nnz_A <- (m + 1) * m / 2 + (n - m) * (m + 1) # num of non-zero
  A_row_inds <- rep(0, nnz_A)
  A_col_inds <- rep(0, nnz_A)
  A_vals <- rep(0, nnz_A)
  ind <- 1
  # iteration through rows
  for (i in 1:n) {
    nnz_A_i <- min(m + 1, i)
    A_row_inds[ind:(ind + nnz_A_i - 1)] <- i
    # first col ind is i
    A_col_inds[ind:(ind + nnz_A_i - 1)] <- rslt$nn[i, 1:nnz_A_i]
    # first A[i, i] should be zero
    A_vals[ind] <- 0
    if (i > 1) {
      A_vals[(ind + 1):(ind + nnz_A_i - 1)] <-
        rslt$cond_mean_coeff[i, 1:min(i - 1, m)]
    }
    ind <- ind + nnz_A_i
  }
  # create sparse A
  rslt$A <- Matrix::sparseMatrix(
    i = A_row_inds, j = A_col_inds, x = A_vals,
    dims = c(n, n)
  )
  return(rslt)
}
