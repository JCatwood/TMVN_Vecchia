library(VeccTMVN)
library(sf)
library(spData)
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
smoothness <- 1.5
neglk_func <- function(covparms, ...) {
  if (any(covparms < c(0.01, 1e-6, 1e-6, 0.01))) {
    return(Inf)
  }
  set.seed(123)
  covparms_with_smooth <- c(covparms[1:3], smoothness, covparms[4])
  negloglk <- -loglk_censor_MVN(
    locs_scaled, ind_censor, y_scaled, b_scaled, cov_name,
    covparms_with_smooth, ...
  )
  cat("covparms is", covparms, "\n")
  cat("Neg loglk is", negloglk, "\n")
  return(negloglk)
}
# opt_obj <- optim(
#   par = covparms_init, fn = neglk_func,
#   control = list(trace = 1), m = 50, NLevel2 = 1e3
# )
# if (!file.exists("results")) {
#   dir.create("results")
# }
# save(opt_obj, file = paste0(
#   "results/PCE_modeling.RData"
# ))
# Find State information for locs -----------------------------------
lonlat_to_state <- function(locs) {
  ## State DF
  states <- spData::us_states
  ## Convert points data.frame to an sf POINTS object
  pts <- st_as_sf(locs, coords = 1:2, crs = 4326)
  ## Transform spatial data to some planar coordinate system
  ## (e.g. Web Mercator) as required for geometric operations
  states <- st_transform(states, crs = 3857)
  pts <- st_transform(pts, crs = 3857)
  ## Find names of state (if any) intersected by each point
  state_names <- states[["NAME"]]
  ii <- as.integer(st_intersects(pts, states))
  state_names[ii]
}
state_names <- lonlat_to_state(data.frame(locs))
ind_Texas <- which(state_names == "Texas")
ind_Texas_neighbor <- which(
  state_names %in%
    c("Texas", "Oklahoma", "New Mexico", "Arkansas", "Louisiana")
)
# Define a enveloping box for Texas -------------------------------------------
## TX boundary lat: 25.83333 to 36.5 lon: -93.51667 to -106.6333
ind_Texas_box <- which((locs[, 1] < -93.02) & (locs[, 1] > -107.13) &
  (locs[, 2] > 25.33) & (locs[, 2] < 37))
# Sample at locations given by `ind_Texas_big` --------------------------------
load("results/PCE_modeling.RData")
ind_Texas_big <- ind_Texas
ind_obs <- which(!is.na(y))
ind_Texas_big <- union(ind_Texas_big, ind_obs)
ind_censor_Texas_big <- which(is.na(y[ind_Texas_big]))
locs_scaled_Texas_big <- locs_scaled[ind_Texas_big, , drop = F]
y_scaled_Texas_big <- y_scaled[ind_Texas_big]
b_scaled_Texas_big <- b_scaled[ind_Texas_big]
N <- 1000
time_sim_Texas_big <- system.time(
  samp_Texas_big <- ptmvrandn(
    locs_scaled_Texas_big,
    ind_censor_Texas_big,
    y_scaled_Texas_big,
    b_scaled_Texas_big, cov_name, c(opt_obj$par[1:3], smoothness, opt_obj$par[4]),
    m = 50, N = N
  )
)[[3]]
if (!file.exists("results")) {
  dir.create("results")
}
save(samp_Texas_big, time_sim_Texas_big,
  file = "results/PCE_modeling_sample.RData"
)
# plots -------------------------------------------
