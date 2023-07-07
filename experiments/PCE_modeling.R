library(VeccTMVN)
rm(list = ls())

load("PCE.RData")
summary(data.PCE.censored)
unique(data.PCE.censored$dl_units)
unique(data.PCE.censored$Units)
hist(log(data.PCE.censored$result_va)) # it looks like censored normal indeed!
# extract raw data -------------------------------
y <- log(data.PCE.censored$result_va + 1e-8)
b_censor <- log(data.PCE.censored$detection_level + 1e-8)
ind_censor <- which(data.PCE.censored$left_censored)
locs <- cbind(
  data.PCE.censored$lon, data.PCE.censored$lat,
  data.PCE.censored$startDate
)
n <- nrow(locs)
d <- ncol(locs)
# standardize -------------------------------
b_scaled <- (b_censor - mean(y, na.rm = T)) / sd(y, na.rm = T)
y_scaled <- (y - mean(y, na.rm = T)) / sd(y, na.rm = T)
locs_scaled <- locs
for (i in 1:d) {
  locs_scaled[, i] <- (locs[, i] - min(locs[, i])) /
    (max(locs[, i]) - min(locs[, i]))
}
# model fitting -------------------------------
cov_name <- "matern15_isotropic"
covparms_init <- c(1, 0.1, 0.01)
lk_func <- function(covparms, ...) {
  set.seed(123)
  loglk_censor_MVN(
    locs_scaled, ind_censor, y_scaled, b_scaled, cov_name,
    covparms, ...
  )
}
range_search <- seq(from = 0.01, to = 0.1, by = 0.01)
lk_search <- rep(0, length(range_search))
ind <- 1
for (range in range_search) {
  covparms <- covparms_init
  covparms[2] <- range
  lk_search[ind] <- lk_func(covparms, m = 30, NLevel2 = 1e3)
  ind <- ind + 1
}
# save results ------------------------------------
save(range_search, lk_search, file = paste0(
  "results/PCE_modeling.RData"
))
