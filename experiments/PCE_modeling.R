library(VeccTMVN)
library(sf)
library(spData)
library(GpGp)
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
locs_scaled[, 1:(d - 1)] <-
  (locs_scaled[, 1:(d - 1)] - min(locs_scaled[, 1:(d - 1)])) /
    (max(locs_scaled[, 1:(d - 1)]) - min(locs_scaled[, 1:(d - 1)]))
locs_scaled[, d] <- (locs_scaled[, d] - min(locs_scaled[, d])) /
  (max(locs_scaled[, d]) - min(locs_scaled[, d]))
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
opt_obj <- optim(
  par = covparms_init, fn = neglk_func,
  control = list(trace = 1), m = 50, NLevel2 = 1e3
)
if (!file.exists("results")) {
  dir.create("results")
}
save(opt_obj, file = paste0(
  "results/PCE_modeling.RData"
))
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
covparms <- c(opt_obj$par[1:3], smoothness, opt_obj$par[4])
time_sim_Texas_big <- system.time(
  samp_Texas_big <- ptmvrandn(
    locs_scaled_Texas_big,
    ind_censor_Texas_big,
    y_scaled_Texas_big,
    b_scaled_Texas_big, cov_name, covparms,
    m = 50, N = N
  )
)[[3]]
if (!file.exists("results")) {
  dir.create("results")
}
save(samp_Texas_big, time_sim_Texas_big, covparms, cov_name, ind_Texas_big,
  file = "results/PCE_modeling_sample.RData"
)
# plots -------------------------------------------
library(fields)
library(RColorBrewer)
load("results/PCE_modeling_sample.RData")
if (!file.exists("plots")) {
  dir.create("plots")
}
lon_grid <- seq(from = -106.6, to = -93.5, length.out = 100)
lat_grid <- seq(from = 25.8, to = 36.6, length.out = 100)
TX_grid <- expand.grid(lon_grid, lat_grid)
colnames(TX_grid) <- c("lon", "lat")
TX_grid_scaled <- TX_grid
TX_grid_scaled[, 1:(d - 1)] <-
  (TX_grid_scaled[, 1:(d - 1)] - min(locs[, 1:(d - 1)])) /
    (max(locs[, 1:(d - 1)]) - min(locs[, 1:(d - 1)]))
TX_grid_scaled <- cbind(TX_grid_scaled, rep(1, nrow(TX_grid)))
colnames(TX_grid_scaled) <- c("lon", "lat", "time")
## pred_VMET --------------------------------------
cov_mat <- get(cov_name)(covparms, rbind(
  locs_scaled_Texas_big,
  as.matrix(TX_grid_scaled)
))
n_obs <- nrow(locs_scaled_Texas_big)
n_grid <- nrow(TX_grid)
pred_VMET <- rowMeans(cov_mat[(n_obs + 1):(n_obs + n_grid), 1:n_obs] %*%
  solve(cov_mat[1:n_obs, 1:n_obs], samp_Texas_big))
pred_VMET <- pred_VMET * sd(y, na.rm = T) + mean(y, na.rm = T)
pred_VMET[!(lonlat_to_state(TX_grid) == "Texas") |
  (is.na(lonlat_to_state(TX_grid)))] <- NA
## pred_GP --------------------------------------
locs_obs_scaled <- locs_scaled[-ind_censor, , drop = F]
y_obs_scaled <- y_scaled[-ind_censor]
cov_mat <- get(cov_name)(covparms, rbind(
  locs_obs_scaled,
  as.matrix(TX_grid_scaled)
))
n_obs <- nrow(locs_obs_scaled)
n_grid <- nrow(TX_grid)
pred_GP <- rowMeans(cov_mat[(n_obs + 1):(n_obs + n_grid), 1:n_obs] %*%
  solve(cov_mat[1:n_obs, 1:n_obs], y_obs_scaled))
pred_GP <- pred_GP * sd(y, na.rm = T) + mean(y, na.rm = T)
pred_GP[!(lonlat_to_state(TX_grid) == "Texas") |
  (is.na(lonlat_to_state(TX_grid)))] <- NA
## actual plot ----------------------------------------
pdf(file = "plots/PCE_modeling.pdf", width = 12, height = 5)
par(mfrow = c(1, 2), mar = c(4, 4, 1, 6))
image.plot(lon_grid, lat_grid, matrix(pred_VMET, 100, 100),
  col = colorRampPalette(brewer.pal(11, "RdBu")[11:1])(30),
  xlab = "longitude", ylab = "latitude", cex.lab = 1.3,
  cex.axis = 1.3, legend.shrink = 0.8, legend.cex = 2.5, legend.width = 2,
  mgp = c(2, 1, 0), zlim = c(-13.1, 0.07)
)
points(
  locs[ind_Texas_big[ind_censor_Texas_big], 1],
  locs[ind_Texas_big[ind_censor_Texas_big], 2],
  col = "grey",
  cex = 0.6, pch = 1,
)
points(
  x = locs[ind_obs, 1], y = locs[ind_obs, 2], col = "black",
  cex = 0.6, pch = 4,
)
# fields::US(xlim = c(-106.6, -93.5), ylim = c(25.8, 36.6), add = T)
image.plot(lon_grid, lat_grid, matrix(pred_GP, 100, 100),
  col = colorRampPalette(brewer.pal(11, "RdBu")[11:1])(30),
  xlab = "longitude", ylab = "latitude", cex.lab = 1.3,
  cex.axis = 1.3, legend.shrink = 0.8, legend.cex = 2.5, legend.width = 2,
  mgp = c(2, 1, 0), zlim = c(-13.1, 0.07)
)
points(
  x = locs[ind_obs, 1], y = locs[ind_obs, 2], col = "black",
  cex = 0.6, pch = 4,
)
dev.off()
