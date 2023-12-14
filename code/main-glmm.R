# main-glmm.R: script to fit a generalized linear mixed model to the 
#              acoustic vocalization data. 
# Authors: Vladimir Kovalenko and Jeffrey W. Doser
rm(list = ls())
# Load libraries ----------------------------------------------------------
library(tidyverse)
library(spAbundance)
library(coda)
# Set seed to get the same results as presented in the manuscript.
set.seed(111)

# Set working directory as needed
# setwd()

# Read in the count and acoustic data -------------------------------------
c <- as.matrix(read.csv("data/acoustic-vocalization-data.csv")[-1,-1])
# Get site names
site.names <- read.csv("data/acoustic-vocalization-data.csv", header = TRUE)[-1, 1]
site.names <- str_extract(site.names, "[^_]+")

# Occupancy covariates ----------------------------------------------------
# Cones
cones <- c(read.csv("data/acoustic-covariates-2.csv")[, c(2)])
cones_bin <- ifelse(cones==0,0,1)

# Live basal area

livebasalarea <- c(read.csv("data/acoustic-covariates-1.csv")[, c(3)])
totalbasalarea <- c(read.csv("data/acoustic-covariates-1.csv")[, c(4)])
meanlivedbh <- c(read.csv("data/acoustic-covariates-1.csv")[, c(5)])
propinfected <- c(read.csv("data/acoustic-covariates-2.csv")[, c(6)])
dougfir <- as.integer((read.csv("data/acoustic-covariates-2.csv")[, c(7)]))

# Site group
site.group.ind <- c(read.csv("data/site-indx-acoustic.csv")[,2])
n.sitegroups <- length(unique(site.group.ind))

X.lambda <- as.matrix(read.csv("data/year-acoustic.csv")[,-1])
X.lambda[, 2] <- ifelse(X.lambda[, 2] == 2020, 0,
                        ifelse(X.lambda[, 2] == 2021,1,2))
head(X.lambda)

# Detection covariates ----------------------------------------------------
# R: total number of sites
R <- nrow(c)
n.count <- apply(c, 1, function(a) sum(!is.na(a)))
# Date --------------------------------
date <- as.matrix(read.csv("data/acoustic-det-covs.csv"))
# Convert date to Julian date
date <- as.POSIXct(date, format = "%m/%d/%Y")
date <- as.numeric(format(date, "%j"))
str(date) # This is the Julian day of the year, or the specific day of the year (1 = Jan 1)
# Put it in matrix format
date <- matrix(date, R, max(n.count))

# Format data for spAbundance ---------------------------------------------
# Covariates
covs <- list(date = date, 
             total.basal.area = totalbasalarea, 
             year = X.lambda[, 2], 
             cones = cones, 
             live.basal.area = livebasalarea, 
             site.group.ind = site.group.ind,
             mean.live.dbh = meanlivedbh,
             prop.infected = propinfected,
             doug.fir = dougfir,
             cones.bin = cones_bin)

# Remove sites without.glmm any data
site.miss.indx <- c(7,10,18,28,30)
#site.miss.indx <- which(apply(c, 1, function(a) sum(is.na(a)) == ncol(c)))
y <- c[-site.miss.indx, ]
covs$date <- matrix(covs$date[-site.miss.indx, ], nrow(y), ncol(y))
covs$total.basal.area <- matrix(covs$total.basal.area[-site.miss.indx], nrow(y), ncol(y))
covs$year <- matrix(covs$year[-site.miss.indx], nrow(y), ncol(y))
covs$cones <- matrix(covs$cones[-site.miss.indx], nrow(y), ncol(y))
covs$live.basal.area <- matrix(covs$live.basal.area[-site.miss.indx], nrow(y), ncol(y))
covs$site.group.ind <- matrix(covs$site.group.ind[-site.miss.indx], nrow(y), ncol(y))
covs$mean.live.dbh <- matrix(covs$mean.live.dbh[-site.miss.indx], nrow(y), ncol(y))
covs$prop.infected <- matrix(covs$prop.infected[-site.miss.indx], nrow(y), ncol(y))
covs$doug.fir <- matrix(covs$doug.fir[-site.miss.indx], nrow(y), ncol(y))
covs$cones.bin <- matrix(covs$cones.bin[-site.miss.indx], nrow(y), ncol(y))
site.names <- site.names[-site.miss.indx]

plot(propinfected,meanlivedbh)

# Data list required for spAbundance
data.list <- list(y = matrix(as.numeric(y), nrow(y), ncol(y)), covs = covs)
save(data.list, site.names, file = 'data/spAbundance-data.rda')

# Get summary statistics for recording days at each site during each year
tmp <- apply(data.list$y, 1, function(a) sum(!is.na(a)))
mean(tmp)
range(tmp)
sd(tmp)



# Tuning values
tuning.list <- list(beta = 0.1, kappa = 0.3, beta.star = 0.5)

# Fit a GLMM with the acoustic data ---------------------------------------
n.samples <- 60000
batch.length <- 25
n.batch <- n.samples / batch.length
n.burn <- 30000
n.thin <- 20
n.chains <- 3

out.glmm <- abund(formula = ~ scale(year) + scale(cones) + 
                              scale(mean.live.dbh) + scale(prop.infected) + scale(date) + 
			      I(scale(date)^2) + (1 | site.group.ind) + 
			      (scale(date) | site.group.ind) + 
			      (I(scale(date)^2) | site.group.ind), 
             data = data.list, 
	     n.batch = n.batch,
	     family = 'NB',
	     tuning = tuning.list,
	     batch.length = batch.length,
	     n.burn = n.burn, 
	     n.thin = n.thin, 
	     n.chains = n.chains, 
             n.report = 10) 

summary(out.glmm)

# This is the full model object, which is too large for GitHub, so it is not
# in the associated repository.
save(out.glmm, 'results/full-glmm-vocalization-results.rda')

# Bayesian p-value
ppc.out.glmm.ft <- ppcAbund(out.glmm, fit.stat = 'freeman-tukey', group = 0)
summary(ppc.out.glmm.ft)
ppc.out.glmm.chi <- ppcAbund(out.glmm, fit.stat = 'chi-square', group = 0)
summary(ppc.out.glmm.chi)
bpv.glmm.ft <- mean(ppc.out.glmm.ft$fit.y < ppc.out.glmm.ft$fit.y.rep)
bpv.glmm.chi <- mean(ppc.out.glmm.chi$fit.y < ppc.out.glmm.chi$fit.y.rep)

# Calculate 95% coverage rate
y.rep.samples <- fitted(out.glmm)
y.rep.means <- apply(y.rep.samples, c(2, 3), mean, na.rm = TRUE)
y.rep.low <- apply(y.rep.samples, c(2, 3), quantile, 0.025, na.rm = TRUE)
y.rep.high <- apply(y.rep.samples, c(2, 3), quantile, 0.975, na.rm = TRUE)
y.true <- out.glmm$y
y.coverage <- ifelse(c(y.true) <= c(y.rep.high) & c(y.true) >= c(y.rep.low), 1, 0)
mean(y.coverage, na.rm = TRUE)

# Compare average site-level vocalizations with the true and replicate data
y.site.rep.means <- apply(y.rep.means, 1, mean, na.rm = TRUE)
y.site.true.means <- apply(y.true, 1, mean, na.rm = TRUE)
plot(y.site.true.means, y.site.rep.means, pch = 19)
abline(0, 1)


# Save a smaller version of the model results to include on GitHub
out.glmm$like.samples <- NULL
out.glmm$y.rep.samples <- NULL
# Mean, 2.5% quantile, and 97.5% quantile of the abundance values at each site/year 
# combination.
abund.means <- apply(out.glmm$mu.samples, 2, mean, na.rm = TRUE)
abund.low <- apply(out.glmm$mu.samples, 2, quantile, 0.025, na.rm = TRUE)
abund.high <- apply(out.glmm$mu.samples, 2, quantile, 0.975, na.rm = TRUE)
out.glmm$mu.samples <- NULL
# Save resulting object to hard drive -------------------------------------
save(out.glmm, abund.means, abund.low, abund.high, y.coverage,
     y.site.rep.means, y.site.true.means, bpv.glmm.ft, bpv.glmm.chi,
     file = 'results/glmm-vocalization-results.rda')
