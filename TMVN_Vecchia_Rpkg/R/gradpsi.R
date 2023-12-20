# Gradient of the psi function in Idea V
#
# @param xAndBeta a vector whose first (n - 1) coeffs are x and second (n - 1)
# coeffs are beta. (n - 1) instead of n because beta_n = 0 and x_n does
# not matter. n is the dim of the MVN prob
# @param veccCondMeanVarObj contains information of the conditional mean
# coefficient, the conditional variance, and the NN array of the Vecchia
# approximation
# @param a lower bound vector for TMVN
# @param b upper bound vector for TMVN
# @return
#   a vector of length 2n - 2, representing the gradient of psi w.r.t. x[-n]
#   and beta[-n]
#
grad_idea5 <- function(xAndBeta, veccCondMeanVarObj, a, b) {
  n <- length(a)
  x <- rep(0, n)
  beta <- rep(0, n)
  x[1:(n - 1)] <- xAndBeta[1:(n - 1)]
  beta[1:(n - 1)] <- xAndBeta[n:(2 * n - 2)]
  D <- sqrt(veccCondMeanVarObj$cond_var)
  mu_c <- as.vector(veccCondMeanVarObj$A %*% x)
  a_tilde_shift <- (a - mu_c) / D - beta
  b_tilde_shift <- (b - mu_c) / D - beta
  log_diff_cdf <- TruncatedNormal::lnNpr(a_tilde_shift, b_tilde_shift)
  pl <- exp(-0.5 * a_tilde_shift^2 - log_diff_cdf) / sqrt(2 * pi)
  pu <- exp(-0.5 * b_tilde_shift^2 - log_diff_cdf) / sqrt(2 * pi)
  Psi <- pl - pu

  dpsi_dx <- as.vector(
    -t(diag(rep(1, n - 1)) - veccCondMeanVarObj$A[-n, -n]) %*%
      ((beta / D)[1:(n - 1)]) +
      t(veccCondMeanVarObj$A[, -n]) %*% (Psi / D)
  )
  dpsi_dbeta <- (beta - (x - mu_c) / D + Psi)[-n]
  c(dpsi_dx, dpsi_dbeta)
}


# Jacobian of the psi function in Idea V
#
# @param xAndBeta a vector whose first (n - 1) coeffs are x and second (n - 1)
# coeffs are beta. (n - 1) instead of n because beta_n = 0 and x_n does
# not matter. n is the dim of the MVN prob
# @param veccCondMeanVarObj contains information of the conditional mean
# coefficient, the conditional variance, and the NN array of the Vecchia
# approximation
# @param a lower bound vector for TMVN
# @param b upper bound vector for TMVN
# @return
#   a square matrix of dim 2n - 2
#
jac_idea5 <- function(xAndBeta, veccCondMeanVarObj, a, b) {
  n <- length(a)
  x <- rep(0, n)
  beta <- rep(0, n)
  x[1:(n - 1)] <- xAndBeta[1:(n - 1)]
  beta[1:(n - 1)] <- xAndBeta[n:(2 * n - 2)]
  D <- sqrt(veccCondMeanVarObj$cond_var)
  mu_c <- as.vector(veccCondMeanVarObj$A %*% x)
  a_tilde_shift <- (a - mu_c) / D - beta
  b_tilde_shift <- (b - mu_c) / D - beta
  log_diff_cdf <- TruncatedNormal::lnNpr(a_tilde_shift, b_tilde_shift)
  pl <- exp(-0.5 * a_tilde_shift^2 - log_diff_cdf) / sqrt(2 * pi)
  pu <- exp(-0.5 * b_tilde_shift^2 - log_diff_cdf) / sqrt(2 * pi)
  Psi <- pl - pu

  a_tilde_shift[is.infinite(a_tilde_shift)] <- 0
  b_tilde_shift[is.infinite(b_tilde_shift)] <- 0
  dPsi <- (-Psi^2) + a_tilde_shift * pl - b_tilde_shift * pu
  dpsi_dx_dx <- t(veccCondMeanVarObj$A[, -n]) %*%
    (dPsi / D / D * veccCondMeanVarObj$A[, -n])
  dpsi_dx_dbeta <- t(dPsi / D * veccCondMeanVarObj$A)[-n, -n] -
    t((diag(rep(1, n - 1)) - veccCondMeanVarObj$A[-n, -n]) / D[-n])
  dpsi_dbeta_dbeta <- diag(1 + dPsi[-n])
  rbind(
    cbind(dpsi_dx_dx, dpsi_dx_dbeta),
    cbind(t(dpsi_dx_dbeta), dpsi_dbeta_dbeta)
  )
}


# TEST-------------------------------

# ## copy from TruncatedNormal V2.2.2-------------------------------
# gradpsi_TN <- function(y, L, l, u){
#   # implements grad_psi(x) to find optimal exponential twisting;
#   # assumes scaled 'L' with zero diagonal;
#   d <- length(u);
#   cv <- rep(0,d);
#   x <- cv;
#   mu <- cv
#   x[1:(d-1)] <- y[1:(d-1)];
#   mu[1:(d-1)] <- y[d:(2*d-2)]
#   # compute now ~l and ~u
#   cv[-1] <- L[-1,] %*% x;
#   lt <- l - mu - cv;
#   ut <- u - mu - cv;
#   # compute gradients avoiding catastrophic cancellation
#   w <- TruncatedNormal::lnNpr(lt, ut);
#   pl <- exp(-0.5*lt^2-w)/sqrt(2*pi);
#   pu <- exp(-0.5*ut^2-w)/sqrt(2*pi)
#   P <- pl-pu;
#   # output the gradient
#   dfdx <- -mu[-d] + as.vector(crossprod(P, L[,-d]))
#   dfdm <- mu - x + P
#   grad <- c(dfdx, dfdm[-d])
#   # here compute Jacobian matrix
#   grad
# }
#
# jacpsi_TN <-  function(y, L, l, u){
#   # implements grad_psi(x) to find optimal exponential twisting;
#   # assume scaled 'L' with zero diagonal;
#   d <- length(u);
#   cv <- rep(0,d);
#   x <- cv;
#   mu <- cv
#   x[1:(d-1)] <- y[1:(d-1)];
#   mu[1:(d-1)] <- y[d:(2*d-2)]
#   # compute now ~l and ~u
#   cv[-1] <- L[-1,] %*% x;
#   lt <- l - mu - cv;
#   ut <- u - mu - cv;
#   # compute gradients avoiding catastrophic cancellation
#   w <- TruncatedNormal::lnNpr(lt, ut);
#   pl <- exp(-0.5*lt^2-w)/sqrt(2*pi);
#   pu <- exp(-0.5*ut^2-w)/sqrt(2*pi)
#   P <- pl-pu;
#
#   # here compute Jacobian matrix
#   lt[is.infinite(lt)] <- 0
#   ut[is.infinite(ut)] <- 0
#   dP <- (-P^2) + lt * pl - ut * pu # dPdm
#   DL <- rep(dP,1,d) * L
#   mx <- -diag(d) + DL
#   xx <- crossprod(L, DL)
#   mx <- mx[1:(d-1),1:(d-1)]
#   xx <- xx[1:(d-1),1:(d-1)]
#   if (d>2){
#     Jac <- rbind(cbind(xx, t(mx)), cbind(mx, diag(1+dP[1:(d-1)])))
#   } else {
#     Jac <- rbind(cbind(xx, t(mx)), cbind(mx, 1 + dP[1:(d-1)]))
#   }
#   Jac
# }
#
# ## example MVN probabilities --------------------------------
# library(GpGp)
# library(TruncatedNormal)
# library(nleqslv)
# source("inv_chol.R")
# source("vecc_cond_mean_var.R")
# set.seed(123)
# n1 <- 10
# n2 <- 10
# n <- n1*n2
# locs <- as.matrix(expand.grid((1 : n1) / n1, (1 : n2) / n2))
# covparms <- c(2, 0.3, 0)
# cov_mat <- matern15_isotropic(covparms, locs)
# a_list <- list(rep(-Inf, n), rep(-1, n), -runif(n) * 2)
# b_list <- list(rep(-2, n), rep(1, n), runif(n) * 2)
#
# ## ordering and NN --------------------------------
# m <- 30
# ord <- order_maxmin(locs)
# locs_ord <- locs[ord, , drop = FALSE]
# cov_mat_ord <- matern15_isotropic(covparms, locs_ord)
# a_list_ord <- lapply(a_list, function(x){x[ord]})
# b_list_ord <- lapply(b_list, function(x){x[ord]})
# NNarray <- find_ordered_nn(locs_ord, m = m)
#
# ## Vecchia approx --------------------------------
# U <- get_sp_inv_chol(cov_mat_ord, NNarray)
# cov_mat_Vecc <- solve(U %*% t(U))
# L_Vecc <- t(chol(cov_mat_Vecc))
# vecc_cond_mean_var_obj <- vecc_cond_mean_var(cov_mat_ord, NNarray)
# D_Vecc <- diag(L_Vecc)
# L_Vecc_scaled <- L_Vecc / D_Vecc
# diag(L_Vecc_scaled) <- 0
#
# ## Compare dpsi ------------------------------
# for(i in 1 : length(a_list_ord)){
#   a_ord <- a_list_ord[[i]]
#   b_ord <- b_list_ord[[i]]
#   x0 <- runif(2 * n - 2)
#   y0 <- x0
#   # y0[1 : (n - 1)] <- (x0[1 : (n - 1)] -
#   #                       vecc_cond_mean_var_obj$A[-n, -n] %*% x0[1 : (n - 1)]) /
#   #   D_Vecc[-n]
#   x0[1 : (n - 1)] <- L_Vecc[-n, -n] %*% y0[1 : (n - 1)]
#   a_ord_scaled <- a_ord / D_Vecc
#   b_ord_scaled <- b_ord / D_Vecc
#   grad_TN <- gradpsi_TN(y0, L_Vecc_scaled, a_ord_scaled, b_ord_scaled)
#   grad_idea_5 <- grad_idea5(x0, vecc_cond_mean_var_obj, a_ord, b_ord)
#   err_beta_grad <- max(abs(
#     grad_TN[n : (2 * n - 2)] - grad_idea_5[n : (2 * n - 2)]))
#   err_x_grad <- max(abs(grad_idea_5[1 : (n - 1)] -
#     t((diag(rep(1, n - 1)) - vecc_cond_mean_var_obj$A[-n, -n]) / D_Vecc[-n]) %*%
#     grad_TN[1 : (n - 1)]))
#   cat("err_beta_grad is ", err_beta_grad, "\n")
#   cat("err_beta_grad is ", err_x_grad, "\n")
# }
#
# ## Compare ddpsi/ddbeta ------------------------------
# for(i in 1 : length(a_list_ord)){
#   a_ord <- a_list_ord[[i]]
#   b_ord <- b_list_ord[[i]]
#   x0 <- runif(2 * n - 2)
#   y0 <- x0
#   y0[1 : (n - 1)] <- (x0[1 : (n - 1)] -
#                         vecc_cond_mean_var_obj$A[-n, -n] %*% x0[1 : (n - 1)]) /
#     D_Vecc[-n]
#   a_ord_scaled <- a_ord / D_Vecc
#   b_ord_scaled <- b_ord / D_Vecc
#   jac_TN <- jacpsi_TN(y0, L_Vecc_scaled, a_ord_scaled, b_ord_scaled)
#   jac_idea_5 <- jac_idea5(x0, vecc_cond_mean_var_obj, a_ord, b_ord)
#   err_beta_jac <- max(jac_TN[n : (2 * n - 2), n : (2 * n - 2)] -
#                         jac_idea_5[n : (2 * n - 2), n : (2 * n - 2)])
#   cat("err_beta_jac is ", err_beta_jac, "\n")
# }
#
# ## Compare system solution --------------------------------
# for(i in 3 : length(a_list_ord)){
#   a_ord <- a_list_ord[[i]]
#   b_ord <- b_list_ord[[i]]
#   x0 <- rep(0, 2 * length(a_ord) - 2)
#   D_Vecc <- diag(L_Vecc)
#   L_Vecc_scaled <- L_Vecc / D_Vecc
#   a_ord_scaled <- a_ord / D_Vecc
#   b_ord_scaled <- b_ord / D_Vecc
#   diag(L_Vecc_scaled) <- 0
#   solv_TN <- nleqslv(x0,
#                      fn = gradpsi_TN,
#                      jac = jacpsi_TN,
#                      L = L_Vecc_scaled, l = a_ord_scaled, u = b_ord_scaled,
#                      global = "pwldog",
#                      method = "Newton",
#                      control = list(maxit = 500L))
#   solv_idea_5 <- nleqslv(x0,
#                          fn = grad_idea5,
#                          jac = jac_idea5,
#                          veccCondMeanVarObj = vecc_cond_mean_var_obj,
#                          a = a_ord, b = b_ord,
#                          global = "pwldog",
#                          method = "Newton",
#                          control = list(maxit = 500L))
#   ### Check gradient consistency -------------------------------
#   y0 <- solv_TN$x
#   x0 <- y0
#   x0[1 : (n - 1)] <- L_Vecc[-n, -n] %*% y0[1 : (n - 1)]
#   grad_idea_5 <- grad_idea5(x0, vecc_cond_mean_var_obj, a_ord, b_ord)
#   # grad_TN <- gradpsi_TN(y0, L_Vecc_scaled, a_ord_scaled, b_ord_scaled)
#   cat("Rnage of the gradient error of idea 5 at the solution ",
#       "(after transformation) of TN is", range(grad_idea_5), "\n")
#
#   ### Check jacobian numerically -----------------------------------
#   y0 <- solv_TN$x
#   x0 <- y0
#   x0[1 : (n - 1)] <- L_Vecc[-n, -n] %*% y0[1 : (n - 1)]
#   epsl <- 1e-5
#   x0_epsl <- x0
#   ind_epsl <- 120
#   x0_epsl[ind_epsl] <- x0[ind_epsl] + epsl
#   grad_idea_5 <- grad_idea5(x0, vecc_cond_mean_var_obj, a_ord, b_ord)
#   grad_idea_5_epsl <- grad_idea5(x0_epsl, vecc_cond_mean_var_obj, a_ord, b_ord)
#   jac_row_numerical <- (grad_idea_5_epsl - grad_idea_5) / epsl
#   jac_idea_5 <- jac_idea5(x0, vecc_cond_mean_var_obj, a_ord, b_ord)
#   # jac_TN <- jacpsi_TN(y0, L_Vecc_scaled, a_ord_scaled, b_ord_scaled)
#   cat("Range of error of the ", ind_epsl, "column in Jacobian is",
#       range(jac_row_numerical - jac_idea_5[, ind_epsl]), "\n")
#
#   ### Check solution consistency -------------------------------
#   err_beta_hat <- max(abs(solv_TN$x[n : (2 * n - 2)] -
#                         solv_idea_5$x[n : (2 * n - 2)]))
#   cat("Terminal codes of TN is", solv_TN$termcd, "\n")
#   cat("Terminal codes of idea_5 is", solv_idea_5$termcd, "\n")
#   cat("err_beta_hat is", err_beta_hat, "\n")
# }
