line_search <- function(x, grad, objFn, gradFn, NewtonStep, alpha = 1e-4, 
                        maxStep = Inf, minStep = 0){
  lambda <- 1
  Newton_length <- sqrt(sum(NewtonStep^2))
  if(Newton_length > maxStep){
    NewtonStep <- maxStep / Newton_length * NewtonStep
    Newton_length <- maxStep
  }
  init_slope <- sum(grad * NewtonStep)
  if(init_slope >= 0)
    stop("NewtonStep isn't a descending direction\n")
  f <- objFn(x)
  while (T) {
    x_new <- x + lambda * p
    f_new <- objFn(x_new)
    if(f_new <= f + alpha * lambda * init_slope)
      return(list(code = 0, ))
  }
}
