# main-count-glmm.R: fits a generalized linear mixed model to the count 
#                    data. 

setwd("~/R/Occupancy/Bayesian/DoserAbundance/main-count-glmm")

rm(list = ls())
# Load libraries ----------------------------------------------------------
library(tidyverse)
# Can install using devtools::install_github("doserjef/spAbundance")
library(spAbundance)
library(coda)
# A useful package for summarizing MCMC results
library(MCMCvis)
library(viridis)

# Read in the count and acoustic data -------------------------------------
# Note that a "site" is the combination of a sampling location and year
# c: point count data with sites as rows and repeat visits as columns
c <- as.matrix(read.csv("v_20-21-22_narm.csv")[-1,-1])
# Get site names
site.names <- read.csv("v_20-21-22_narm.csv", header = TRUE)[-1, 1]
site.names <- str_extract(site.names, "[^_]+")


# Occupancy covariates ----------------------------------------------------
# Cones
cones <- c(read.csv("x_acoustic_narm_conemod.csv")[, c(2)])
cones_bin <- ifelse(cones==0,0,1)

# Live basal area

livebasalarea <- c(read.csv("x_acoustic_narm.csv")[, c(3)])
totalbasalarea <- c(read.csv("x_acoustic_narm.csv")[, c(4)])
meanlivedbh <- c(read.csv("x_acoustic_narm.csv")[, c(5)])
propinfected <- c(read.csv("x_acoustic_narm.csv")[, c(6)])
dougfir <- as.integer((read.csv("x_acoustic_narm_conemod.csv")[, c(7)]))

# Site group
# ditto with this. Converted it to a vector
# site.group.ind <- as.matrix(read.csv("c_modC.csv")[,-c(1:4,6:8)])
# NOTE: double check that these are correct. Just figured its worth a check
#       since they are nearly all in order, except the 2nd sitegroup.
site.group.ind <- c(read.csv("site_groups_narm.csv")[,2])
site.group.ind
n.sitegroups <- length(unique(site.group.ind))

# X.lambda: design matrix consisting of an intercept in column 1 and an 
#           indicator variable for year in column 2.
X.lambda <- as.matrix(read.csv("X.lambda_acoustic_narm.csv")[,-1])
X.lambda[, 2] <- ifelse(X.lambda[, 2] == 2020, 0,
                        ifelse(X.lambda[, 2] == 2021,1,2))
head(X.lambda)

# Detection covariates ----------------------------------------------------
# R: total number of sites
R <- nrow(c)
# n.count: number of point count visits for each of the J possible sites.
n.count <- apply(c, 1, function(a) sum(!is.na(a)))
# Date --------------------------------
date <- as.matrix(read.csv("v_20-21-22_dates_narm.csv"))
# Convert date to Julian date
date <- as.POSIXct(date, format = "%m/%d/%Y")
date <- as.numeric(format(date, "%j"))
str(date) # This is the Julian day of the year, or the specific day of the year (1 = Jan 1)
# Put it in matrix format
date <- matrix(date, R, max(n.count))

# Format data for spAbundance ---------------------------------------------
# Get detection-nondetection data from count data
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

chart.Correlation(cbind(cones,cones_bin,meanlivedbh,livebasalarea,propinfected,as.factor(dougfir)))

plot(propinfected,meanlivedbh)

# JWD: the problem was occurring because the "y" that you were putting into
#      data.list was a character matrix, not a numeric matrix, and that led to 
#      an error in ppcAbund. Below I convert it to a numeric matrix, and then 
#      everything works as desired. 
data.list <- list(y = matrix(as.numeric(y), nrow(y), ncol(y)), covs = covs)

# JWD: you can increase speed to convergence by including the "tuning" 
#      values below. Details on these are in the spOccupancy vignettes. 
#      Essentially, they are values that help the algorithm "find" the 
#      correct values.
tuning.list <- list(beta = 0.1, kappa = 0.3, beta.star = 0.5)

# Fit a GLMM with the point count data ------------------------
# JWD: I just messed around with this for testing just to make stuff faster. 
# n.samples <- 50000
# batch.length <- 25
# n.batch <- n.samples / batch.length
# n.burn <- 25000
# n.thin <- 25
# n.samples <- 20000
# batch.length <- 25
# n.batch <- n.samples / batch.length
# n.burn <- 2500
# n.thin <- 20
# JWD: I just changed for testing
n.samples <- 60000
batch.length <- 25
n.batch <- n.samples / batch.length
n.burn <- 30000
n.thin <- 20
n.chains <- 3

out.glmm <- abund(formula = ~ scale(year) + scale(cones) +                          scale(mean.live.dbh) + scale(prop.infected) + scale(date) + 
			 I(scale(date)^2) + 
			 (1 | site.group.ind) + (scale(date) | site.group.ind) + 
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

# Plot the posterior densities for the covariate effects
beta.names <- colnames(out.glmm$beta.samples)[-1]
plot.df <- data.frame(val = c(out.glmm$beta.samples[, -1]), 
                      parameter = rep(beta.names, each = out.glmm$n.post * out.glmm$n.chains))

ggplot(data = plot.df, aes(x = val, fill = parameter), color = 'black') + 
  geom_vline(xintercept = 0, col = 'black', lty = 2) + 
  geom_posterior(draw_ci = TRUE, ci_width = 0.95, draw_sd = FALSE, brighten = FALSE, midline = 'black', 
                 color = 'black', alpha = 0.85) + 
  scale_fill_viridis_d() + 
  scale_y_continuous(breaks = c(1, 3, 5, 7, 9, 11), 
    labels = c('Date','Quadratic Date','Year','Prop Infected','Live DBH','Cone Density')) + 
  theme_classic(base_size = 12) + 
  labs(x = 'Effect Size', y = 'Parameter') +
  guides(fill = 'none') +
  scale_color_viridis_d()

####################

apply(out.glmm$beta.samples, 2, function(a) mean(a > 0))
apply(out.glmm$beta.samples, 2, function(a) mean(a < 0))


# Get estimated relative abundance at each site/year combination
abund.means <- apply(out.glmm$mu.samples, 2, mean, na.rm = TRUE)
abund.low <- apply(out.glmm$mu.samples, 2, quantile, 0.025, na.rm = TRUE)
abund.high <- apply(out.glmm$mu.samples, 2, quantile, 0.975, na.rm = TRUE)
tmp <- data.list$covs$year[, 1]
tmp <- ifelse(tmp == 0, 2020, ifelse(tmp == 1, 2021, 2022))
plot.df <- data.frame(abund.mean = abund.means, 
		      abund.low = abund.low,
		      abund.high = abund.high,
		      year = factor(tmp),
		      site = site.names)

plot.df <- arrange(plot.df,(plot.df$abund.mean))

plot.df.wide <- plot.df %>%
  pivot_wider(names_from = year, values_from = c(abund.mean, abund.low, abund.high), 
              names_sep = '.')
site.indx <- 1:nrow(plot.df.wide)
plot.df.wide$index.2020 <- site.indx - 0.25
plot.df.wide$index.2021 <- site.indx 
plot.df.wide$index.2022 <- site.indx + 0.25

# tmp <- plot.df %>%
#   group_by(site) %>%
#   summarize(val = mean(abund.mean))
# site.indx <- order(tmp$val)

my.cols <- viridis(3)
names(my.cols) <- c('2020', '2021', '2022')
p2 <- ggplot(data = plot.df.wide) + 
  geom_point(aes(x = index.2020, y = abund.mean.2020, col = '2020'), size = 2) +
  geom_segment(aes(x = index.2020, y = abund.low.2020, xend = index.2020, 
		   yend = abund.high.2020, col = '2020'), lineend = 'butt', size = 1) +
  geom_point(aes(x = index.2021, y = abund.mean.2021, col = '2021'), size = 2) +
  geom_segment(aes(x = index.2021, y = abund.low.2021, xend = index.2021, 
		   yend = abund.high.2021, col = '2021'), lineend = 'butt', size = 1) +
  geom_point(aes(x = index.2022, y = abund.mean.2022, col = '2022'), size = 2) +
  geom_segment(aes(x = index.2022, y = abund.low.2022, xend = index.2022, 
		   yend = abund.high.2022, col = '2022'), lineend = 'butt', size = 1) +
  theme_classic(base_size = 12) + 
  scale_color_manual(values = my.cols) + 
  scale_x_continuous(breaks = site.indx, labels = plot.df.wide$site, name = 'Site') + 
  labs(y = 'Expected number of vocalizations', color = 'Year') +
  coord_flip() + theme(axis.title.y = element_blank(),axis.text.y = element_text(size= 10, color = "black", hjust = 0.5), plot.margin = unit(c(1,1,1,1),"mm"))


# Generate a plot of average number of vocalizations across all sites per year
mu.samples.big <- array(NA, dim = c(nrow(out.glmm$mu.samples), length(site.indx) * 3, 
				    ncol(out.glmm$y)))
mu.samples.big[, -site.miss.indx, ] <- out.glmm$mu.samples
mu.samples.big <- aperm(mu.samples.big, c(1, 3, 2))
mu.samples.big <- array(mu.samples.big, c(nrow(out.glmm$mu.samples), 
					  ncol(out.glmm$y), 3, length(site.indx)))

plot.df <- data.frame(abund.median = apply(mu.samples.big, 3, median, na.rm = TRUE),
		      abund.lowest = apply(mu.samples.big, 3, quantile, 
					   0.025, na.rm = TRUE),
		      abund.low = apply(mu.samples.big, 3, quantile, 
					0.25, na.rm = TRUE),
		      abund.high = apply(mu.samples.big, 3, quantile, 
					 0.75, na.rm = TRUE),
		      abund.highest = apply(mu.samples.big, 3, quantile, 
					    0.975, na.rm = TRUE),
		      year = factor(2020:2022))
# The box is the 50% credible interval. The lines are the 95% credible interval.
p2 <- ggplot(data = plot.df, aes(x = year, fill=year)) +
  geom_boxplot(aes(ymin = abund.lowest, lower = abund.low, middle = abund.median,
		   upper = abund.high, ymax = abund.highest),
	       stat = 'identity', size = 0.75, alpha = 0.75, width = 0.5) +
  theme_bw(base_size = 12) +
  scale_fill_viridis_d() + 
  labs(x = '', y = 'Average number of vocalizations') #+ ylim(0,60)

# Boxplot of covariate effects side-by-side -------------------------------
# Covariate effects are stored in "out$beta.samples"
beta.quants <- apply(out.glmm$beta.samples, 2, quantile, c(0.025, 0.25, 0.5, 0.75, 0.975))
beta.quants
# Probability each effect is positive
prob.pos.beta <- apply(out.glmm$beta.samples, 2, function(a) mean(a > 0))
# Note the -1 just removes the intercept, I also removed the date (linear and 
# quadratic) so that its the same variables that are shown for occupancy, but
# you can put those back in if you want.
plot.df <- data.frame(med = beta.quants[3, -c(1)], 
		      lowest = beta.quants[1, -c(1)],
		      low = beta.quants[2, -c(1)],
		      high = beta.quants[4, -c(1)],
		      highest = beta.quants[5, -c(1)],
		      param = c('Year', 'Cone Density','Live Basal Area', 'Live Tree Diameter', 'Proportion Infected','Date','Date (Quadratic)'), 
                      prob.positive = prob.pos.beta[-c(1)])
# The boxplots are colored by the probability the effect of that covariate is positive
# such that values with a high probability of a positive effect are blue, and with a low
# probability of being positive (aka a high prob of being negative) are red.
ggplot(data = plot.df, aes(x = param)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_boxplot(aes(ymin = lowest, lower = low, middle = med, 
		   upper = high, ymax = highest), 
	       stat = 'identity', size = 0.75, width = 0.5) + 
  theme_bw(base_size = 12) + 
  scale_fill_gradient2(midpoint = 0.5, low = '#B2182B', mid = 'white', high = '#2166AC',
    	               na.value = NA) +
  guides(fill = 'none') +
  labs(x = 'Vocalization Covariate', y = 'Effect Size') + theme(axis.text.x = element_text(angle = -45, vjust = 0, hjust=0.2)) 

# Marginal effects plot of covariate effect -------------------------------
# CHANGE HERE
pred.vals <- seq(min(data.list$covs$mean.live.dbh),
                 max(data.list$covs$mean.live.dbh),
                 length.out = 100)

# Scale predicted values by mean and standard deviation used to fit the model
# CHANGE HERE
pred.vals.scale <- (pred.vals - mean(c(data.list$covs$mean.live.dbh))) /
  sd(c(data.list$covs$mean.live.dbh))
# Predict occupancy across live.basal.area values at mean values of all other variables
# JWD: note that these need to be in the same order as specified in the "formula" argument
#      when running the model.
# CHANGE HERE
pred.df <- as.matrix(data.frame(intercept = 1, year = 0, cones = 0,
                                mean.live.dbh = pred.vals.scale, prop.infected = 0,
                                date = 0, date2 = 0))
pred.array <- array(NA, dim = c(nrow(pred.df), 1, ncol(pred.df)))
pred.array[, 1, ] <- pred.df
# Make a copy for doing prediction with the goal of producing a marginal probability plot
# like below. This essentially equates to generating the prediction for a site that
# has the average spatial random effect value. 
# JWD: get rid of this, this was only necessary when including spatial random effects.
# out.copy <- out.glmm
# class(out.copy) <- 'tPGOcc'
# JWD: also need to assign names to the third dimension in the array for the 
#      specific variables
dimnames(pred.array)[[3]] <- c('(Intercept)', 'scale(year)', 'scale(cones)', 
			       'scale(mean.live.dbh)', 'scale(prop.infected)', 
			       'scale(date)', 'I(scale(date)^2)')
out.pred <- predict(out.glmm, pred.array, ignore.RE = TRUE)
str(out.pred)
# JWD: note all the "psi" values are now "mu"
mu.0.quants <- apply(out.pred$mu.0.samples, 2, quantile,
                      prob = c(0.025, 0.5, 0.975))
mu.plot.dat <- data.frame(mu.med = mu.0.quants[2, ],
                           mu.low = mu.0.quants[1, ],
                           mu.high = mu.0.quants[3, ],
                           mean.live.dbh = pred.vals)
ggplot(mu.plot.dat, aes(x = 2.5*mean.live.dbh, y = mu.med)) +
  geom_ribbon(aes(ymin = mu.low, ymax = mu.high), fill = 'grey70') +
  geom_line() +
  theme_bw() +
  labs(x = 'Mean Live DBH (cm)', y = 'Expected Number of Vocalizations', size = 12)


