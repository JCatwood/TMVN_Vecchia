library(truncnorm)

mvnrnd_wrap <- function(a, b, mu, NN, veccObj, N, verbose = 0) {
  n <- length(a)
  tmp <- cbind(a, b, mu)
  a <- tmp[, 1]
  b <- tmp[, 2]
  mu <- tmp[, 3]
  # find tilting parameter beta -----------------------------------
  trunc_expect <- truncnorm::etruncnorm(a, b, mean = mu)
  x0 <- c(trunc_expect, rep(0, n))
  solv_idea_5_sp <- stats::optim(
    x0,
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
    veccCondMeanVarObj = veccObj,
    a = a, b = b, mu = mu, verbose = verbose,
    lower = c(a, rep(-Inf, n)), upper = c(b, rep(Inf, n)),
    control = list(maxit = 500)
  )
  if (verbose) {
    cat(
      "Gradient norm at the optimal beta is", sqrt(2 * solv_idea_5_sp$value),
      "\n"
    )
  }
  if (any(solv_idea_5_sp$par[1:n] < a) ||
    any(solv_idea_5_sp$par[1:n] > b)) {
    warning("Optimal x is outside the integration region during minmax tilting\n")
  }
  x_star <- solv_idea_5_sp$par[1:n]
  beta <- solv_idea_5_sp$par[(n + 1):(2 * n)]
  # 2nd opt for finding x_star ---------------------------
  solv_xstar <- stats::optim(
    x_star,
    fn = function(x, ...) {
      -psi_wrapper(x, ...)
    },
    gr = function(x, ...) {
      -dpsi_dx(x, ...)
    },
    method = "L-BFGS-B",
    beta = beta,
    veccCondMeanVarObj = veccObj,
    a = a, b = b, mu = mu, NN = NN,
    lower = c(a, rep(-Inf, n)), upper = c(b, rep(Inf, n)),
    control = list(maxit = 500)
  )
  if (verbose) {
    cat("Psi value is", -solv_xstar$value, "\n")
  }
  if (any(solv_xstar$par < a) ||
    any(solv_xstar$par > b)) {
    warning("Optimal x is outside the integration region during minmax tilting\n")
  }
  x_star <- solv_xstar$par
  # psi_star --------------------------------------
  psi_star <- psi_wrapper(x_star,
    beta = beta,
    veccCondMeanVarObj = veccObj,
    a = a, b = b, mu = mu, NN = NN
  )
  # sample until N samples are collected ------------------
  X <- matrix(0, nrow = n, ncol = N)
  accept <- 0L
  iter <- 0L
  n_sim <- floor(N / 2) * 2
  ntotsim <- 0L
  while (accept < N) { # while # of accepted is less than N
    call <- mvnrnd(
      a, b, NN, mu, veccObj$cond_mean_coeff,
      sqrt(veccObj$cond_var), beta, n_sim
    )
    ntotsim <- ntotsim + n_sim
    idx <- stats::rexp(n_sim) > (psi_star - call$logpr) # acceptance tests
    m <- sum(idx)
    if (m > N - accept) {
      m <- N - accept
      idx <- which(idx)[1:m]
    }
    if (m > 0) {
      X[, (accept + 1):(accept + m)] <-
        t(call$X_trans[idx, , drop = F]) # accumulate accepted
    }
    accept <- accept + m # keep track of # of accepted
    iter <- iter + 1L # keep track of while loop iterations
    if (verbose) {
      cat("Iteration", iter, "n_sim", n_sim, "n_accept", m, "\n")
    }
    n_sim <- min(c(1e6, N, ceiling(n_sim / m * (N - accept))))
    n_sim <- ceiling(n_sim / 2) * 2
    if ((ntotsim > 1e4) && (accept / ntotsim < 1e-3)) { # if iterations are getting large, give warning
      warning("Acceptance probability smaller than 0.001")
    } else if (iter > 1e5) { # if iterations too large, seek approximation only
      if (accept == 0) {
        stop("Could not sample from truncated Normal - check input")
      } else if (accept > 1) {
        X <- X[, 1:accept]
        warning("Sample of size smaller than N returned.")
      }
    }
  }
  return(X)
}


# TEST -------------------------------------------------------
# library(GpGp)
# library(VeccTMVN)
# library(truncnorm)
# library(TruncatedNormal)
# ## example MVN probabilities --------------------------------
# set.seed(123)
# n1 <- 10
# n2 <- 10
# n <- n1 * n2
# locs <- as.matrix(expand.grid((1:n1) / n1, (1:n2) / n2))
# covparms <- c(2, 0.3, 0.01)
# cov_name <- "matern15_isotropic"
# cov_mat <- get(cov_name)(covparms, locs)
# a <- rep(-Inf, n)
# b <- rep(-2, n)
# mu <- rep(0, n)
# N <- 1e3
# ## Sample with TruncatedNormal -----------------------
# samp_TN <- TruncatedNormal::mvrandn(a, b, cov_mat, n = N, mu = mu)
# ## Vecc approx objs --------------------------------
# m <- 30 # num of nearest neighbors
# ord <- VeccTMVN::Vecc_reorder(
#   a, b, m, locs, cov_name,
#   covparms
# )$order
# a_vecc_ord <- a[ord]
# b_vecc_ord <- b[ord]
# locs_vecc_ord <- locs[ord, , drop = F]
# NN <- GpGp::find_ordered_nn(locs_vecc_ord, m)
# vecc_obj <- vecc_cond_mean_var_sp(NN,
#                                   locs = locs_vecc_ord, covName = cov_name,
#                                   covParms = covparms
# )
# ## Sample with mvnrnd_wrap -----------------------
# samp_Vecc_ord <- mvnrnd_wrap(
#   a_vecc_ord, b_vecc_ord, mu,
#   NN = NN, veccObj = vecc_obj, N = N, verbose = 1
# )
# ord_rev <- integer(n)
# ord_rev[ord] <- 1 : n
# samp_Vecc <- samp_Vecc_ord[ord_rev, , drop = F]
# ## Visual comparison of two TMVN samples -------------------------------
# par(mfrow = c(1, 2))
# hist(samp_Vecc)
# hist(samp_TN)
