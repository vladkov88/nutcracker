# main-spatial-occ.R: script to fit a multi-season spatial occupancy model 
#                     to assess Clark's Nutcracker use of whitebark pine
#                     forests.
# Authors: Vladimir Kovalenko and Jeffrey W. Doser
rm(list = ls())
# Load libraries ----------------------------------------------------------
library(tidyverse)
library(spOccupancy)
library(coda)
library(sf)
# Set seed to get same exact results as presented in the paper
set.seed(100)

# Set working directory as needed
# setwd()

# Read in the count and acoustic data -------------------------------------
# c: point count data with sites as rows and repeat visits as columns
c <- as.matrix(read.csv("data/point-count-data.csv")[, c(2:4)])

# Occupancy covariates ----------------------------------------------------
# Cones
cones <- c(read.csv("data/covariate-point-count-data.csv")[, c(2)])

# Live basal area
livebasalarea <- c(read.csv("data/covariate-point-count-data.csv")[, c(3)])
totalbasalarea <- c(read.csv("data/covariate-point-count-data.csv")[, c(4)])
meanlivedbh <- c(read.csv("data/covariate-point-count-data.csv")[, c(5)])

propinfected <- c(read.csv("data/covariate-point-count-data.csv")[, c(6)])

# Site group
site.group.ind <- c(read.csv("data/point-count-data.csv")[,5])
n.sitegroups <- length(unique(site.group.ind))

# Extract the year associated with each data point.
X.lambda <- as.matrix(read.csv("data/year-point-count.csv")[,-1])
X.lambda[, 2] <- ifelse(X.lambda[, 2] == 2020, 0,
                        ifelse(X.lambda[, 2] == 2021,1,2))
head(X.lambda)

# Detection covariates ----------------------------------------------------
# R: total number of sites
R <- nrow(c)
# n.count: number of point count visits for each of the J possible sites.
n.count <- apply(c, 1, function(a) sum(!is.na(a)))
# Date --------------------------------
date <- as.matrix(read.csv("data/point-count-data.csv")[,c(6,7,8)])
# Convert date to Julian date
date <- as.POSIXct(date, format = "%m/%d/%Y")
date <- as.numeric(format(date, "%j"))
str(date) # This is the Julian day of the year, or the specific day of the year (1 = Jan 1)
# Put it in matrix format
date <- matrix(date, R, max(n.count))

# Wind --------------------------------
wind <- read.csv("data/point-count-data.csv")[, c(9,10,11)]
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
coords <- read.csv("data/coordinates.csv")
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
save(data.list, site.names, file = 'data/spOccupancy-data.rda')

# Get size of study area (i.e., convex hull around 100m buffer of the points)
coords.buffered <- st_buffer(coords.sf, dist = 0.1)
st_area(st_convex_hull(st_union(coords.buffered)))

# Set up for occupancy model ----------------------------------------------
n.samples <- 60000
batch.length <- 25
n.batch <- n.samples / batch.length 
n.burn <- 30000
n.thin <- 15
n.chains <- 3
# Generate informative prior for the spatial range ------------------------
# Restrict the lower bound of the spatial autocorrelation to be the 
# diameter of the maximum core use area from 3 birds that were tagged with 
# GPS. 
dist.matrix <- dist(coords.aea)
# Radius of the maximum core use area for the three birds is 2.18082km, so will
# use 2 * 2.18082 for the lower bound of the spatial autocorrelation. Then
# use the maximum inter-site distance as the upper bound.
diameter.core <- 2.18082 * 2
priors <- list(phi.unif = c(3 / max(dist.matrix), 3 / diameter.core))

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

# Save model results to hard drive ----------------------------------------
save(out, file = 'results/spatial-occ-results.rda')
