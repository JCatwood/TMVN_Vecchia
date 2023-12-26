# Line search along Newton step
#
# @param x current point
# @param grad gradient of `objFn` at `x`
# @param objFn objective function
# @param NewtonStep the Newton step, $- H^{-1}g$. Practically, any direction
# with a negative inner product with grad will work
# @param alpha sufficient decrease coefficient
# @param maxStep maximum step length
# @param stepTol absolute convergence parameter
# @return a list of `code` and `x_new`
# code 1: qualifying `x_new` found
# code 2: qualifying `x_new` not found, return the original input `x`
# code 3: `NewtonStep` is not a descending direction
# @examples
# Fn <- function(x) {
#   c(x[1]^2 + x[2]^2 - 2, exp(x[1] - 1) + x[2]^3 - 2)
# }
# J <- function(x) {
#   matrix(c(
#     2 * x[1], 2 * x[2],
#     exp(x[1] - 1), 3 * x[2]^2
#   ), 2, 2, byrow = T)
# }
# obj_fn <- function(x) {
#   0.5 * sum(Fn(x)^2)
# }
# x <- c(2, 0.5)
# for (i in 1:10) {
#   Newton_step <- -as.vector(solve(J(x)) %*% Fn(x))
#   grad <- as.vector(t(J(x)) %*% Fn(x))
#   ret <- line_search(x, grad, obj_fn, Newton_step)
#   if (ret$code == 0) {
#     x <- ret$x_new
#   } else {
#     break
#   }
# }
# cat("Solution is", x, "where f(x) is", obj_fn(x), "\n")
line_search <- function(x, grad, objFn, NewtonStep, alpha = 1e-4,
                        maxStep = Inf, stepTol = 1e-4) {
  lambda <- 1
  Newton_length <- sqrt(sum(NewtonStep^2))
  if (lambda * Newton_length < stepTol) {
    return(list(code = 2, x_new = x))
  }
  if (Newton_length > maxStep) {
    NewtonStep <- maxStep / Newton_length * NewtonStep
    Newton_length <- maxStep
  }
  init_slope <- sum(grad * NewtonStep)
  if (init_slope >= 0) {
    warning("NewtonStep isn't a descending direction\n")
    return(list(code = 3, x_new = x))
  }
  f <- objFn(x)
  while (T) {
    x_new <- x + lambda * NewtonStep
    f_new <- objFn(x_new)
    if (f_new <= f + alpha * lambda * init_slope) {
      return(list(code = 1, x_new = x_new))
    } else if (lambda * Newton_length < stepTol) {
      return(list(code = 2, x_new = x))
    } else {
      if (lambda == 1) {
        lambda_tmp <- (-init_slope) / 2 / (f_new - f - init_slope)
      } else {
        a_b <- 1 / (lambda - lambda_prev) * as.vector(
          matrix(
            c(
              1 / lambda^2, -1 / (lambda_prev^2),
              -lambda_prev / (lambda^2), lambda / (lambda_prev^2)
            ),
            2, 2,
            byrow = T
          ) %*% c(
            f_new - f - lambda * init_slope,
            f_prev - f - lambda_prev * init_slope
          )
        )
        a <- a_b[1]
        b <- a_b[2]
        disc <- b^2 - 3 * a * init_slope
        if (a == 0) {
          lambda_tmp <- -init_slope / (2 * b)
        } else {
          lambda_tmp <- (-b + sqrt(disc)) / (3 * a)
        }
        if (lambda_tmp > 0.5 * lambda) {
          lambda_tmp <- 0.5 * lambda
        }
      }
      lambda_prev <- lambda
      f_prev <- f_new
      if (lambda_tmp <= 0.1 * lambda) {
        lambda <- 0.1 * lambda
      } else {
        lambda <- lambda_tmp
      }
    }
  }
}


# Locally defined function for solving non-linear systems
#
# @param x0 initial guess
# @param fn non-linear system, fn: $R^{n} -> R^{n}$
# @param jacTransFn function that computes Jacobian^{T} cdot fn
# @param jacInvFn function that computes Jacobian^{T} cdot fn
# @param ... arguments to `fn`, `jacTransFn`, and `jacInvFn`
# @param control a list of tuning parameters for optimization
# @return a list of `code` and `x`
# code 1: qualifying `x` found
# code 2: qualifying `x` not found after reaching `maxit`
# code 3: qualifying `x` not found before reaching `maxit`
# @examples
# Fn <- function(x) {
#   c(x[1]^2 + x[2]^2 - 2, exp(x[1] - 1) + x[2]^3 - 2)
# }
# J <- function(x) {
#   matrix(c(
#     2 * x[1], 2 * x[2],
#     exp(x[1] - 1), 3 * x[2]^2
#   ), 2, 2, byrow = T)
# }
# jac_trans_fn <- function(x) {
#   as.vector(t(J(x)) %*% Fn(x))
# }
# jac_inv_fn <- function(x) {
#   as.vector(solve(J(x)) %*% Fn(x))
# }
# x <- c(2, 0.5)
# ret <- my_nleqslv(x, Fn, jac_trans_fn, jac_inv_fn)
# cat(
#   "Solution is", ret$x, "where f(x) is", ret$obj, "and code is", ret$code,
#   "\n"
# )
my_nleqslv <- function(x0, fn, jacTransFn, jacInvFn, ..., control = list()) {
  x <- x0
  maxit <- control[["maxit"]]
  if (is.null(maxit)) {
    maxit <- 500
  }
  step_tol <- control[["step_tol"]]
  if (is.null(step_tol)) {
    step_tol <- 1e-4
  }
  obj_tol <- control[["obj_tol"]]
  if (is.null(obj_tol)) {
    obj_tol <- 1e-4
  }
  obj_fn <- function(x) {
    0.5 * sum(fn(x, ...)^2)
  }
  for (i in 1:maxit) {
    grad <- jacTransFn(x, ...)
    Newton_step <- -jacInvFn(x, ...)
    ret <- line_search(x, grad, obj_fn, Newton_step, stepTol = step_tol)
    if (ret$code == 1) {
      x <- ret$x_new
    } else if (obj_fn(x) > obj_tol) {
      ret <- line_search(x, grad, obj_fn, -grad, stepTol = step_tol)
      if (ret$code == 1) {
        x <- ret$x_new
      } else {
        break
      }
    }
  }
  obj_ret <- obj_fn(x)
  if (obj_ret <= obj_tol) {
    code_ret <- 1
  } else if (i == maxit) {
    code_ret <- 2
  } else {
    code_ret <- 3
  }
  return(list(x = x, obj = obj_ret, code = code_ret))
}
