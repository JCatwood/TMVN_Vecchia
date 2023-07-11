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
locs_names <- c("lon", "lat", "date")
colnames(locs) <- locs_names
# standardize -------------------------------
b_scaled <- (b_censor - mean(y, na.rm = T)) / sd(y, na.rm = T)
y_scaled <- (y - mean(y, na.rm = T)) / sd(y, na.rm = T)
locs_scaled <- locs
for (i in 1:d) {
  locs_scaled[, i] <- (locs[, i] - min(locs[, i])) /
    (max(locs[, i]) - min(locs[, i]))
}
# # variogram analysis ----------------------
# library(sp)
# library(gstat)
# mydf <- data.frame(cbind(locs_scaled, y_scaled))
# mydf <- na.omit(mydf)
# colnames(mydf) <- c(locs_names, "log_PCE")
# # mydf$lon <- 1
# # mydf$lat <- 1
# coordinates(mydf) = ~lon + lat + date
# myvariog <- variogram(log_PCE~1, mydf, cutoff=2.0, width = 0.01)
# plot(myvariog)
# # analysis on temporal dependence --------------------------
# library(tidyr)
# library(dplyr)
# mytib <- as_tibble(mydf)
# avg_by_date <- mytib %>%
#   filter(lon < 0.55) %>%
#   filter(lat < 0.55) %>%
#   filter(lon > 0.45) %>%
#   filter(lat > 0.45) %>%
#   group_by(date) %>% summarize(mean = mean(log_PCE))
# plot(avg_by_date)  # it seems that the dependence between avg_y and date is weak
# library(ggplot2)
# library(patchwork)
# date_unique <- unique(mytib$date)
# for(i in date_unique) {
#   date_i <- date_unique[i]
#   mytib %>% filter(date < 0.1) %>%
#     # filter(date > 0.01) %>%
#     ggplot(aes(x = lon, y = lat, color=log_PCE)) +
#     geom_point()
# }
# model fitting -------------------------------
cov_name <- "matern_spacetime"
covparms_init <- c(1, 0.1, 0.1, 0.5) # var, range1, range2, nugget
neglk_func <- function(covparms, ...) {
  if (any(covparms < c(0.01, 1e-6, 1e-6, 0.01))) {
    return(Inf)
  }
  set.seed(123)
  smoothness <- 1.5
  covparms_with_smooth <- c(covparms_init[1:3], smoothness, covparms_init[4])
  -loglk_censor_MVN(
    locs_scaled, ind_censor, y_scaled, b_scaled, cov_name,
    covparms_with_smooth, ...
  )
}
opt_obj <- optim(
  par = covparms_init, fn = neglk_func,
  control = list(trace = 1), m = 50, NLevel2 = 1e4
)
# save results ------------------------------------
save(opt_obj, file = paste0(
  "results/PCE_modeling.RData"
))
# plots -------------------------------------------
