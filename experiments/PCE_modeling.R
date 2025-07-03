library(VeccTMVN)
library(sf)
library(spData)
library(GpGp)
library(scoringRules)
rm(list = ls())

# preliminary data analysis -------------------------------
load("PCE.RData")
summary(data.PCE.censored)
unique(data.PCE.censored$dl_units)
unique(data.PCE.censored$Units)
hist(log(data.PCE.censored$result_va)) # it looks like censored normal indeed!
# extract raw data -------------------------------
y <- log(data.PCE.censored$result_va + 1e-8) # log precipitation with a small nugget
b_censor <- log(data.PCE.censored$detection_level + 1e-8) # log censoring thresholds with a small nugget
ind_censor <- which(data.PCE.censored$left_censored) # indices for censored responses
locs <- cbind(
  data.PCE.censored$lon, data.PCE.censored$lat,
  data.PCE.censored$startDate
)  # spatial locations for the spatial covariance matrix
n <- nrow(locs)  # number of locations, also the dimension of the MVN probability
d <- ncol(locs)  # dimension of the spatial domain, where the Matern kernel is defined
locs_names <- c("lon", "lat", "date")
colnames(locs) <- locs_names
# standardize -------------------------------
# scale y and b by the mean and standard deivation of the observed log-precipitation
# scale locs to the unit hypercube
b_scaled <- (b_censor - mean(y, na.rm = T)) / sd(y, na.rm = T)
y_scaled <- (y - mean(y, na.rm = T)) / sd(y, na.rm = T)
locs_scaled <- locs
locs_scaled[, 1:(d - 1)] <-
  (locs_scaled[, 1:(d - 1)] - min(locs_scaled[, 1:(d - 1)])) /
    (max(locs_scaled[, 1:(d - 1)]) - min(locs_scaled[, 1:(d - 1)]))
locs_scaled[, d] <- (locs_scaled[, d] - min(locs_scaled[, d])) /
  (max(locs_scaled[, d]) - min(locs_scaled[, d]))
# data splitting ---------------------------------------
mask_censor <- rep(FALSE, n)
mask_censor[ind_censor] <- TRUE
n_train <- round(n * 0.8)
set.seed(123)
# randomly generated indices for training locs
ind_train <- sample(1:n, n_train, replace = FALSE)
y_train <- y[ind_train]
y_scaled_train <- y_scaled[ind_train]
y_scaled_test <- y_scaled[-ind_train]
b_scaled_train <- b_scaled[ind_train]
b_scaled_test <- b_scaled[-ind_train]
locs_scaled_train <- locs_scaled[ind_train, ]
locs_scaled_test <- locs_scaled[-ind_train, ]
locs_train <- locs[ind_train, ]
locs_test <- locs[-ind_train, ]
# model specification -------------------------------
cov_name <- "matern_spacetime"
covparms_init <- c(1, 0.1, 0.1, 0.5) # var, range1, range2, nugget
smoothness <- 1.5
# model fitting -------------------------------
# negative log-likelihood function that optimizes the covariance parameters
neglk_func <- function(covparms, ...) {
  if (any(covparms < c(0.01, 1e-6, 1e-6, 0.01))) {
    return(Inf)  # prevent near-singular situations
  }
  set.seed(123)
  # augment covparms with smoothness that is fixed (not optimized)
  covparms_with_smooth <- c(covparms[1:3], smoothness, covparms[4])
  # use the loglk_censor_MVN function from VeccTMVN to estimate the log-likelihood
  # of a partially censored GP
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
# save the results. Once you have the results, you can comment out the 
# model fitting section to obtain the results faster

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
state_names_train <- lonlat_to_state(data.frame(locs_train))
# find indices for training locations that are in Texas
ind_Texas_train <- which(state_names_train == "Texas")

# Define a enveloping box for Texas -------------------------------------------
## TX boundary lat: 25.83333 to 36.5 lon: -93.51667 to -106.6333
## find indices for training locations that envelops Texas
ind_Texas_box_train <- which((locs_train[, 1] < -93.02) & (locs_train[, 1] > -107.13) &
  (locs_train[, 2] > 25.33) & (locs_train[, 2] < 37))
# Sample at locations given by `ind_Texas_big_train` --------------------------------
load("results/PCE_modeling.RData")  # load previously trained covariance parameters
ind_Texas_big_train <- ind_Texas_train
ind_obs_train <- which(!is.na(y_train))
# conditioning all observed responses, including those outside Texas, to sample
# the censored responses within Texas
ind_Texas_big_train <- union(ind_Texas_big_train, ind_obs_train)
ind_censor_Texas_big_train <- which(is.na(y_train[ind_Texas_big_train]))
# subset the training data to those indexed by `ind_Texas_big_train`
locs_scaled_Texas_big_train <- locs_scaled_train[ind_Texas_big_train, , drop = F]
y_scaled_Texas_big_train <- y_scaled_train[ind_Texas_big_train]
b_scaled_Texas_big_train <- b_scaled_train[ind_Texas_big_train]
N <- 1000
covparms <- c(opt_obj$par[1:3], smoothness, opt_obj$par[4])
# sample the censored responses ----------------------
time_sim_Texas_big_train <- system.time(
  samp_Texas_big_train <- ptmvrandn(
    locs_scaled_Texas_big_train,
    ind_censor_Texas_big_train,
    y_scaled_Texas_big_train,
    b_scaled_Texas_big_train, cov_name, covparms,
    m = 50, N = N, reorder = FALSE
  )
)[[3]]  # use ptmvrandn function from VeccTMVN to sample a partially censored GP
if (!file.exists("results")) {
  dir.create("results")
}
save(samp_Texas_big_train, time_sim_Texas_big_train, covparms,
     cov_name, ind_Texas_big_train,
  file = "results/PCE_modeling_sample.RData"
)
# save the results. Once you have the results, you can comment out the 
# sample the censored responses section to obtain the results faster

# MSE at testing ----------------------------------
load("results/PCE_modeling_sample.RData")
state_names_test <- lonlat_to_state(data.frame(locs_test))
# find the indices for those located in Texas among the testing dataset
ind_Texas_test <- which(state_names_test == "Texas")
# extract subsets correspond to those in Texas 
locs_scaled_Texas_test <- locs_scaled_test[ind_Texas_test, ]
y_scaled_Texas_test <- y_scaled_test[ind_Texas_test]
b_scaled_Texas_test <- b_scaled_test[ind_Texas_test]
mask_censored_Texas_test <- is.na(y_scaled_Texas_test)
# for censored responses in the testing dataset, they are converted to binary variables
# e.g., below the threshold or not to enable to computation of RMSE
score_benchmark <- y_scaled_Texas_test
score_benchmark[mask_censored_Texas_test] <- 1
# construct the covariance matrix for training and testing locations that 
# correspond to Texas
cov_mat <- get(cov_name)(covparms, rbind(
  locs_scaled_Texas_big_train,
  as.matrix(locs_scaled_Texas_test)
))
n_obs <- nrow(locs_scaled_Texas_big_train)
n_test_Texas <- nrow(locs_scaled_Texas_test)
## pred_VMET --------------------------------------
# ordinary kriging based on each sample of the partially censored GP
pred_VMET <- cov_mat[(n_obs + 1):(n_obs + n_test_Texas), 1:n_obs] %*%
  solve(cov_mat[1:n_obs, 1:n_obs], samp_Texas_big_train)
# convert the predicted value at locations with censored responses to binary
# e.g., below the threshold or not
pred_VMET[mask_censored_Texas_test, ] <-
  pred_VMET[mask_censored_Texas_test, ] <=
    b_scaled_Texas_test[mask_censored_Texas_test]
# CRPS_VMET <- mean(crps_sample(score_benchmark, pred_VMET))
pred_VMET <- rowMeans(pred_VMET)
RMSE_VMET <- sqrt(mean((pred_VMET - score_benchmark)^2))
## pred_LOD_GP --------------------------------------
# do the same for level-of-detection supplemented GP 
# substitute the censored responses by their level-of-detection
y_aug <- samp_Texas_big_train[, 1]
y_aug[ind_censor_Texas_big_train] <-
  b_scaled_Texas_big_train[ind_censor_Texas_big_train]
pred_GP_aug <- as.vector(cov_mat[(n_obs + 1):(n_obs + n_test_Texas), 1:n_obs] %*%
  solve(cov_mat[1:n_obs, 1:n_obs], y_aug))  # ordinary kriging
pred_GP_aug[mask_censored_Texas_test] <-
  pred_GP_aug[mask_censored_Texas_test] <=
    b_scaled_Texas_test[mask_censored_Texas_test]  # convert to binary
RMSE_GP_aug <- sqrt(mean((pred_GP_aug - score_benchmark)^2))
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
  locs_scaled_Texas_big_train,
  as.matrix(TX_grid_scaled)
))
n_obs <- nrow(locs_scaled_Texas_big_train)
n_grid <- nrow(TX_grid)
pred_VMET <- rowMeans(cov_mat[(n_obs + 1):(n_obs + n_grid), 1:n_obs] %*%
  solve(cov_mat[1:n_obs, 1:n_obs], samp_Texas_big_train))
pred_VMET <- pred_VMET * sd(y, na.rm = T) + mean(y, na.rm = T)
pred_VMET[!(lonlat_to_state(TX_grid) == "Texas") |
  (is.na(lonlat_to_state(TX_grid)))] <- NA
## pred_LOD_GP --------------------------------------
y_aug <- samp_Texas_big_train[, 1]
y_aug[ind_censor_Texas_big_train] <-
  b_scaled_Texas_big_train[ind_censor_Texas_big_train]
cov_mat <- get(cov_name)(covparms, rbind(
  locs_scaled_Texas_big_train,
  as.matrix(TX_grid_scaled)
))
n_obs <- nrow(locs_scaled_Texas_big_train)
n_grid <- nrow(TX_grid)
pred_GP_aug <- cov_mat[(n_obs + 1):(n_obs + n_grid), 1:n_obs] %*%
  solve(cov_mat[1:n_obs, 1:n_obs], y_aug)
pred_GP_aug <- pred_GP_aug * sd(y, na.rm = T) + mean(y, na.rm = T)
pred_GP_aug[!(lonlat_to_state(TX_grid) == "Texas") |
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
  locs_train[ind_Texas_big_train[ind_censor_Texas_big_train], 1],
  locs[ind_Texas_big_train[ind_censor_Texas_big_train], 2],
  col = "grey",
  cex = 0.6, pch = 1,
)
points(
  x = locs_train[ind_obs_train, 1], y = locs_train[ind_obs_train, 2], col = "black",
  cex = 0.6, pch = 4,
)
# fields::US(xlim = c(-106.6, -93.5), ylim = c(25.8, 36.6), add = T)
image.plot(lon_grid, lat_grid, matrix(pred_GP_aug, 100, 100),
  col = colorRampPalette(brewer.pal(11, "RdBu")[11:1])(30),
  xlab = "longitude", ylab = "latitude", cex.lab = 1.3,
  cex.axis = 1.3, legend.shrink = 0.8, legend.cex = 2.5, legend.width = 2,
  mgp = c(2, 1, 0), zlim = c(-13.1, 0.07)
)
points(
  locs_train[ind_Texas_big_train[ind_censor_Texas_big_train], 1],
  locs[ind_Texas_big_train[ind_censor_Texas_big_train], 2],
  col = "grey",
  cex = 0.6, pch = 1,
)
points(
  x = locs_train[ind_obs_train, 1], y = locs_train[ind_obs_train, 2], col = "black",
  cex = 0.6, pch = 4,
)
dev.off()
