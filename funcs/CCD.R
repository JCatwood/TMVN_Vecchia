#' Constrained coordinate descent for a quadratic function under rectangular
#' bounds
#' $$ \frac{1}{2} x^\top A x + b^\top x $$ s.t. lb <= x <= ub
#' A should be positive definite
#'
cons_coord_desc <- function(A, b = rep(0, nrow(A)), lb = rep(-Inf, nrow(A)),
                            ub = rep(Inf, nrow(A)), silent = T, convtol = 1e-3,
                            maxEpoch = 20) {
  parms_old <- ub
  parms_old[is.na(ub)] <- lb[is.na(ub)]
  parms_old[is.na(parms_old)] <- 0
  n <- nrow(A)
  if (any(diag(A) <= 0)) {
    stop("Input matrix A should be positive definite\n")
  }
  obj_old <- as.numeric(0.5 * t(parms_old) %*% A %*% parms_old +
    t(b) %*% parms_old)
  if (!silent) {
    cat("Epoch 0, objective is", obj_old, "\n")
  }
  parms_new <- parms_old
  for (k in 1:maxEpoch)
  {
    for (j in 1:n)
    {
      opt_j <- (-sum(A[j, -j] * parms_new[-j]) - b[j]) / A[j, j]
      if (opt_j > ub[j]) {
        parms_new[j] <- ub[j]
      } else if (opt_j < lb[j]) {
        parms_new[j] <- lb[j]
      } else {
        parms_new[j] <- opt_j
      }
    }
    obj_new <- as.numeric(0.5 * t(parms_new) %*% A %*% parms_new +
      t(b) %*% parms_new)
    if (!silent) {
      cat("Epoch", k, "objective is", obj_new, "\n")
    }
    if ((obj_old - obj_new) / abs(obj_old) < convtol) {
      cat("Convergence condition met\n")
      return(parms_new)
    }
    parms_old <- parms_new
    obj_old <- obj_new
  }
  return(parms_new)
}
