# main-count-glmm.R: script to fit a generalized linear mixed model to the 
#                    acoustic vocalization data. 
# Authors: Vladimir Kovalenko and Jeffrey W. Doser
rm(list = ls())
# Load libraries ----------------------------------------------------------
library(tidyverse)
# Can install using devtools::install_github("doserjef/spAbundance")
library(spAbundance)
library(coda)
library(MCMCvis)
library(viridis)
# devtools::install_github("iamamutt/ggdistribute")
library(ggdistribute)

# Set working directory as needed
# setwd()

# Read in the count and acoustic data -------------------------------------
c <- as.matrix(read.csv("acoustic-vocalization-data.csv")[-1,-1])
# Get site names
site.names <- read.csv("acoustic-vocalization-data.csv", header = TRUE)[-1, 1]
site.names <- str_extract(site.names, "[^_]+")

# Occupancy covariates ----------------------------------------------------
# Cones
cones <- c(read.csv("acoustic-covariates-2.csv")[, c(2)])
cones_bin <- ifelse(cones==0,0,1)

# Live basal area

livebasalarea <- c(read.csv("acoustic-covariates-1.csv")[, c(3)])
totalbasalarea <- c(read.csv("acoustic-covariates-1.csv")[, c(4)])
meanlivedbh <- c(read.csv("acoustic-covariates-1.csv")[, c(5)])
propinfected <- c(read.csv("acoustic-covariates-2.csv")[, c(6)])
dougfir <- as.integer((read.csv("acoustic-covariates-2.csv")[, c(7)]))

# Site group
site.group.ind <- c(read.csv("site-indx-acoustic.csv")[,2])
n.sitegroups <- length(unique(site.group.ind))

X.lambda <- as.matrix(read.csv("year-acoustic.csv")[,-1])
X.lambda[, 2] <- ifelse(X.lambda[, 2] == 2020, 0,
                        ifelse(X.lambda[, 2] == 2021,1,2))
head(X.lambda)

# Detection covariates ----------------------------------------------------
# R: total number of sites
R <- nrow(c)
n.count <- apply(c, 1, function(a) sum(!is.na(a)))
# Date --------------------------------
date <- as.matrix(read.csv("acoustic-det-covs.csv"))
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

# Tuning values
tuning.list <- list(beta = 0.1, kappa = 0.3, beta.star = 0.5)

# Fit a GLMM with the point count data ------------------------
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

# Generate Figure 4 -------------------------------------------------------
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


# Generate part of Figure 2 -----------------------------------------------
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


# Marginal effects plot of covariate effect -------------------------------
# Example code to create one of the plots in Figure 5
pred.vals <- seq(min(data.list$covs$mean.live.dbh),
                 max(data.list$covs$mean.live.dbh),
                 length.out = 100)

# Scale predicted values by mean and standard deviation used to fit the model
pred.vals.scale <- (pred.vals - mean(c(data.list$covs$mean.live.dbh))) /
  sd(c(data.list$covs$mean.live.dbh))
# Predict occupancy across live.basal.area values at mean values of all other variables
pred.df <- as.matrix(data.frame(intercept = 1, year = 0, cones = 0,
                                mean.live.dbh = pred.vals.scale, prop.infected = 0,
                                date = 0, date2 = 0))
pred.array <- array(NA, dim = c(nrow(pred.df), 1, ncol(pred.df)))
pred.array[, 1, ] <- pred.df
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
