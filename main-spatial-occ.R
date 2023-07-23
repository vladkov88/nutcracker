# main-spatial-occ.R: script to fit a multi-season spatial occupancy model 
#                     to assess Clark's Nutcracker use of whitebark pine
#                     forests.
# Authors: Vladimir Kovalenko and Jeffrey W. Doser
rm(list = ls())
# Load libraries ----------------------------------------------------------
library(tidyverse)
library(viridis)
library(spOccupancy)
library(coda)
library(MCMCvis)
library(sf)
# devtools::install_github("iamamutt/ggdistribute")
library(ggdistribute)

# Set working directory as needed
# setwd()

# Read in the count and acoustic data -------------------------------------
# c: point count data with sites as rows and repeat visits as columns
c <- as.matrix(read.csv("point-count-data.csv")[, c(2:4)])

# Occupancy covariates ----------------------------------------------------
# Cones
cones <- c(read.csv("covariate-point-count-data.csv")[, c(2)])

# Live basal area
livebasalarea <- c(read.csv("covariate-point-count-data.csv")[, c(3)])
totalbasalarea <- c(read.csv("covariate-point-count-data.csv")[, c(4)])
meanlivedbh <- c(read.csv("covariate-point-count-data.csv")[, c(5)])

propinfected <- c(read.csv("covariate-point-count-data.csv")[, c(6)])

# Site group
site.group.ind <- c(read.csv("point-count-data.csv")[,5])
n.sitegroups <- length(unique(site.group.ind))

# Extract the year associated with each data point.
X.lambda <- as.matrix(read.csv("year-point-count.csv")[,-1])
X.lambda[, 2] <- ifelse(X.lambda[, 2] == 2020, 0,
                        ifelse(X.lambda[, 2] == 2021,1,2))
head(X.lambda)

# Detection covariates ----------------------------------------------------
# R: total number of sites
R <- nrow(c)
# n.count: number of point count visits for each of the J possible sites.
n.count <- apply(c, 1, function(a) sum(!is.na(a)))
# Date --------------------------------
date <- as.matrix(read.csv("point-count-data.csv")[,c(6,7,8)])
# Convert date to Julian date
date <- as.POSIXct(date, format = "%m/%d/%Y")
date <- as.numeric(format(date, "%j"))
str(date) # This is the Julian day of the year, or the specific day of the year (1 = Jan 1)
# Put it in matrix format
date <- matrix(date, R, max(n.count))

# Wind --------------------------------
wind <- read.csv("point-count-data.csv")[, c(9,10,11)]
wind <- unlist(c(wind))
wind <- matrix(c(wind), R, max(n.count))

# Format data for spOccupancy ---------------------------------------------
# Get detection-nondetection data from count data
y.mat <- ifelse(c > 0, 1, 0)
# Reformat data into a 3-D array with dimensions site, year, replicate survey.
n.years <- 3
n.sites <- R / n.years
y <- array(NA, dim = c(n.sites, n.years, ncol(y.mat)))
for (i in 1:ncol(y.mat)) {
  y[, , i] <- matrix(y.mat[, i], n.sites, n.years, byrow = TRUE)
}
# Occupancy covariates
occ.covs <- list(year = matrix(X.lambda[, 2], n.sites, n.years, byrow = TRUE), 
		 cones = matrix(cones, n.sites, n.years, byrow = TRUE),
		 live.basal.area = matrix(livebasalarea, n.sites, n.years, byrow = TRUE),
		 site.group.ind = matrix(site.group.ind, n.sites, n.years, byrow = TRUE),
		 mean.live.dbh = matrix(meanlivedbh, n.sites, n.years, byrow = TRUE),
		 prop.infected = matrix(propinfected, n.sites, n.years, byrow = TRUE))

# Detection covariates
date.array <- array(NA, dim = dim(y))
wind.array <- array(NA, dim = dim(y))
for (i in 1:ncol(y.mat)) {
  date.array[, , i] <- matrix(date[, i], n.sites, n.years, byrow = TRUE)
  wind.array[, , i] <- matrix(wind[, i], n.sites, n.years, byrow = TRUE)
}
det.covs <- list(date = date.array,
		 wind = wind.array, 
                 total.basal.area = matrix(totalbasalarea, n.sites, n.years, byrow = TRUE),
                 year = matrix(X.lambda[, 2], n.sites, n.years, byrow = TRUE))

# Load in plot coordinates
coords <- read.csv("coordinates.csv")
site.names <- coords$site
coords <- coords[, c('Lon', 'Lat')]
# Convert to projected coordinate system
coords <- st_as_sf(coords, 
		   coords = c('Lon', 'Lat'), 
		   crs = "+proj=longlat +datum=WGS84")
coords.sf <- coords %>%
  st_transform("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs")
coords.aea <- st_coordinates(coords.sf)
data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs, coords = coords.aea)

# Set up for occupancy model ----------------------------------------------
n.samples <- 60000
batch.length <- 25
n.batch <- n.samples / batch.length 
n.burn <- 30000
n.thin <- 15
n.chains <- 3
# Generate informative prior for the spatial range ------------------------
# Restrict the lower bound of the spatial autocorrelation to be the average
# home range size obtained from the GPS data.
dist.matrix <- dist(coords.aea)
# Average home range area from GPS data: 32.0km2. Use the diameter of this 
# area as the lower bound on the spatial autocorrelation between sites. Then
# use the maximum inter-site distance as the upper bound.
diameter.hr <- sqrt(32 / pi) * 2
priors <- list(phi.unif = c(3 / max(dist.matrix), 3 / diameter.hr))

# Run the spatial multi-season occupancy model ----------------------------
out <- stPGOcc(occ.formula = ~ scale(year) + scale(cones) + scale(live.basal.area) + 
                             scale(mean.live.dbh) + scale(prop.infected), 
               det.formula = ~ scale(date) + I(scale(date)^2) + scale(wind) + 
                               scale(total.basal.area) + scale(year), 
               priors = priors,
               data = data.list, 
	       n.batch = n.batch,
	       batch.length = batch.length,
	       n.burn = n.burn, 
	       n.thin = n.thin, 
	       n.chains = n.chains, 
	       NNGP = TRUE,
	       n.neighbors = 15,
               n.report = 50) 

summary(out)

# Calculate the probability the covariate effects are positive
# (The probability the covariate effects is negative is just 1 - the prob they are positive.)
# A value of 0.5 indicates equal probability of a positive or negative effect, aka
# not a lot of support for an effect. The farther away from .5, the further support for 
# either a positive (close to 1) or negative (close to 0) effect.
# Occurrence
apply(out$beta.samples, 2, function(a) mean(a > 0))
# Detection
apply(out$alpha.samples, 2, function(a) mean(a > 0))

# Probability of a negative effect
# Occurrence
apply(out$beta.samples, 2, function(a) mean(a < 0))

# Detection
apply(out$alpha.samples, 2, function(a) mean(a < 0))

# Generate Figure 3 -------------------------------------------------------
# Plot the posterior densities for the covariate effects
beta.names <- colnames(out$beta.samples)[-1]
plot.df <- data.frame(val = c(out$beta.samples[, -1]), 
		      parameter = rep(beta.names, each = out$n.post * out$n.chains))

ggplot(data = plot.df, aes(x = val, fill = parameter), color = 'black') + 
  geom_vline(xintercept = 0, col = 'black', lty = 2) + 
  geom_posterior(draw_ci = TRUE, ci_width = 0.95, draw_sd = FALSE, brighten = FALSE, midline = 'black', 
		 color = 'black', alpha = 0.85) + 
  scale_fill_viridis_d() + 
  scale_y_continuous(breaks = c(1, 5,8,11.5,14.5), 
		     labels = c('Year','Prop Infected','Live BA','Cone Density','Live DBH')) + 
  theme_classic(base_size = 12) + 
  labs(x = 'Effect Size', y = 'Parameter') +
  guides(fill = 'none') +
  scale_color_viridis_d()

# Detection covariates
# Plot the posterior densities for the covariate effects
alpha.names <- colnames(out$alpha.samples)[-1]
plot.df <- data.frame(val = c(out$alpha.samples[, -1]), 
                      parameter = rep(alpha.names, each = out$n.post * out$n.chains))

ggplot(data = plot.df, aes(x = val, fill = parameter), color = 'black') + 
  geom_vline(xintercept = 0, col = 'black', lty = 2) + 
  geom_posterior(draw_ci = TRUE, ci_width = 0.95, draw_sd = FALSE, brighten = FALSE, midline = 'black', 
                 color = 'black', alpha = 0.85) + 
  scale_fill_viridis_d() + 
  scale_y_continuous(breaks = c(1.5,4.5,7,9.5,12), 
                     labels = c('Quadratic Date','Date','Wind','Year','Total Basal Area')) + 
  theme_classic(base_size = 12) + 
  labs(x = 'Effect Size', y = 'Parameter') +
  guides(fill = 'none') +
  scale_color_viridis_d()

# Generate Figure S3 (plot of detection prob vs. date) --------------------
n.pred <- 500
det.pred.vals <- seq(from = min(c(data.list$det.covs$date), na.rm = TRUE), 
                     to = max(c(data.list$det.covs$date), na.rm = TRUE), length.out = n.pred)
X.p.pred <- matrix(0, n.pred, ncol(out$X.p))
n.post <- out$n.post * out$n.chains
p.pred <- matrix(NA, n.post, n.pred)
# Intercept
X.p.pred[, 1] <- 1
# Linear effect of date
X.p.pred[, 2] <- (det.pred.vals - mean(c(data.list$det.covs$date), na.rm = TRUE)) / 
  sd(c(data.list$det.covs$date), na.rm = TRUE)
# Quadratic effect of date
X.p.pred[, 3] <- X.p.pred[, 2]^2
# Setting all other values to 0 (the mean), which means we're predicting detection 
# probabili
for (i in 1:n.post) {
  p.pred[i, ] <- plogis(X.p.pred %*% as.matrix(out$alpha.samples[i, ]))
}
# You can fine tune this as you like. 
plot.df <- data.frame(mean.val = apply(p.pred, 2, mean),
                      low = apply(p.pred, 2, quantile, 0.025),
                      high = apply(p.pred, 2, quantile, 0.975), 
                      day = as.Date(det.pred.vals, origin = as.Date('2021-01-01')))
ggplot(plot.df, aes(x = day, y = mean.val)) +
  geom_ribbon(aes(ymin = low, ymax = high), fill = 'gray') +
  geom_line() +
  theme_bw() +
  labs(x = 'Date', y = 'Detection Probability')
