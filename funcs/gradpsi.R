source("utils.R")

#' Gradient of the psi function in Idea V
#' Input:
#'   xAndBeta - a vector whose first (n - 1) coeffs are x and second (n - 1)
#'     coeffs are beta. (n - 1) instead of n because beta_n = 0 and x_n does
#'     not matter. n is the dim of the MVN prob
#'   veccCondMeanVarObj - contains information of the conditional mean
#'     coefficient, the conditional variance, and the NN array of the Vecchia
#'     approximation
#'   a - lower bound vector for TMVN
#'   b - upper bound vector for TMVN
#' Return:
#'   a vector of length 2n - 2, representing the gradient of psi w.r.t. beta[-1]
#'     and x[-1]
#'
grad_idea5 <- function(xAndBeta, veccCondMeanVarObj, a, b) {
  n <- length(a)
  x <- xAndBeta[1:(n - 1)]
  beta <- xAndBeta[n:(2 * n - 2)]
  D <- sqrt(veccCondMeanVarObj$cond_var[-1])  # vector of length n - 1
  mu_c <- as.vector(veccCondMeanVarObj$A[-1,-1] %*% x)
  a_tilde_shift <- (a[-1] - mu_c) / D - beta
  b_tilde_shift <- (b[-1] - mu_c) / D - beta
  log_diff_cdf <- lnNpr(a_tilde_shift, b_tilde_shift)
  pl <- exp(-0.5 * a_tilde_shift ^ 2 - w) / sqrt(2 * pi)
  pu <- exp(-0.5 * b_tilde_shift ^ 2 - w) / sqrt(2 * pi)
  Psi <- pl - pu
  dpsi_dx = -D * beta + as.vector(t(A) %*% (Psi / D))
  dpsi_dbeta = beta - (x - mu_c) / D + Psi
  c(dpsi_dx, dpsi_dbeta)
}


#' Jacobian of the psi function in Idea V
#' Input:
#'   xAndBeta - a vector whose first (n - 1) coeffs are x and second (n - 1)
#'     coeffs are beta. (n - 1) instead of n because beta_n = 0 and x_n does
#'     not matter. n is the dim of the MVN prob
#'   veccCondMeanVarObj - contains information of the conditional mean
#'     coefficient, the conditional variance, and the NN array of the Vecchia
#'     approximation
#'   a - lower bound vector for TMVN
#'   b - upper bound vector for TMVN
#' Return:
#'   a square matrix of dim 2n - 2
#'
jac_idea5 <- function(xAndBeta, veccCondMeanVarObj, a, b) {
  n <- length(a)
  x <- xAndBeta[1:(n - 1)]
  beta <- xAndBeta[n:(2 * n - 2)]
  D <- sqrt(veccCondMeanVarObj$cond_var[-1])  # vector of length n - 1
  mu_c <- as.vector(veccCondMeanVarObj$A[-1,-1] %*% x)
  a_tilde_shift <- (a[-1] - mu_c) / D - beta
  b_tilde_shift <- (b[-1] - mu_c) / D - beta
  log_diff_cdf <- lnNpr(a_tilde_shift, b_tilde_shift)
  pl <- exp(-0.5 * a_tilde_shift ^ 2 - w) / sqrt(2 * pi)
  pu <- exp(-0.5 * b_tilde_shift ^ 2 - w) / sqrt(2 * pi)
  Psi <- pl - pu
  
  a_tilde_shift[is.infinite(a_tilde_shift)] <- 0
  b_tilde_shift[is.infinite(b_tilde_shift)] <- 0
  dPsi <- (-Psi ^ 2) + a_tilde_shift * pl - b_tilde_shift * pu
  dpsi_dx_dx <- t(veccCondMeanVarObj$A[-1, -1]) %*%
    (dPsi / D / D * veccCondMeanVarObj$A[-1, -1])
  dpsi_dx_dbeta <- t(dPsi / D * veccCondMeanVarObj$A[-1, -1]) - diag(1 / D)
  dpsi_dbeta_dbeta <- diag(1 + dPsi)
  rbind(cbind(dpsi_dx_dx, dpsi_dx_dbeta), 
        cbind(t(dpsi_dx_dbeta), dpsi_dbeta_dbeta))
}










