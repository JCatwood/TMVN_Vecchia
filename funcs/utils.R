#' Compute log(pnorm(b) - pnorm(a)) with the input mu and sd
lnNpr <- function(a, b, mu = 0, sd = 1){
  logcdf_a <- pnorm(a, mean = mu, sd = sd, log.p = T)
  logcdf_b <- pnorm(b, mean = mu, sd = sd, log.p = T)
  if(logcdf_b - logcdf_a < 1e-10)
    return(logcdf_a + log(b - a))
  else
    return(logcdf_a + log(exp(b - a) - 1))
}