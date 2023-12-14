# summary.R: script to summarize the occupancy model and GLMM results and 
#            generate the resulting figures in the manuscript.
# Authors: Vladimir Kovalenko and Jeffrey W. Doser

rm(list = ls())
library(tidyverse)
library(spOccupancy)
library(spAbundance)
library(coda)
library(MCMCvis)
library(viridis)
# devtools::install_github("iamamutt/ggdistribute")
library(ggdistribute)
library(ggpubr)
library(sf)
library(pals)

# Load model results ------------------------------------------------------
# Occupancy model results (stored in an object called out) 
load("results/spatial-occ-results.rda")
# Relative activity GLMM results (out.glmm, abund.means, abund.low, abund.high)
load("results/glmm-vocalization-results.rda")

# Assess convergence ------------------------------------------------------
# Summary of model with convergence assessment
summary(out)
summary(out.glmm)
# Visual assessment of convergence ----
# Occupancy model occupancy regression coefficients
plot(out, 'beta', density = FALSE)
# Occupancy model detection regression coefficients
plot(out, 'alpha', density = FALSE)
# Occupancy model spatial parameters
plot(out, 'theta')
# GLMM regression coefficients
plot(out.glmm, 'beta', density = FALSE)
# GLMM random effect variances
plot(out.glmm, 'sigma.sq.mu', density = FALSE)

# Summarize point count occupancy model results ---------------------------
# Calculate the probability the covariate effects are positive
# A value of 0.5 indicates equal probability of a positive or negative effect, aka
# not a lot of support for an effect. The farther away from .5, the further support for 
# either a positive (close to 1) or negative (close to 0) effect.
# Occurrence
apply(out$beta.samples, 2, function(a) mean(a > 0))
# Detection
apply(out$alpha.samples, 2, function(a) mean(a > 0))

# Posterior predictive checks ---------------------------------------------
# Occupancy model with point count data
out.ppc.occ.ft.1 <- ppcOcc(out, fit.stat = 'freeman-tukey', group = 1)
out.ppc.occ.ft.2 <- ppcOcc(out, fit.stat = 'freeman-tukey', group = 2)
out.ppc.occ.chi.1 <- ppcOcc(out, fit.stat = 'chi-squared', group = 1)
out.ppc.occ.chi.2 <- ppcOcc(out, fit.stat = 'chi-squared', group = 2)
summary(out.ppc.occ.ft.1)
summary(out.ppc.occ.ft.2)
summary(out.ppc.occ.chi.1)
summary(out.ppc.occ.chi.2)
# GLMM with acoustic vocalization data
bpv.glmm.ft
bpv.glmm.chi
# Abundance coverage
mean(y.coverage, na.rm = TRUE)

# Generate figures --------------------------------------------------------
# Figure 3 ----------------------------
# Plot the posterior densities for the covariate effects
beta.names <- colnames(out$beta.samples)[-1]
plot.df <- data.frame(val = c(out$beta.samples[, -1]),
		      parameter = rep(beta.names, each = out$n.post * out$n.chains))

panel.a <- ggplot(data = plot.df, aes(x = val, group = parameter), color = 'black') +
  geom_vline(xintercept = 0, col = 'black', lty = 2) +
  geom_posterior(draw_ci = TRUE, ci_width = 0.95, draw_sd = FALSE, brighten = FALSE, midline = 'black',
		 color = 'black', alpha = 1, fill = 'grey90') +
  # scale_fill_viridis_d() +
  scale_y_continuous(breaks = c(1, 5,8,11.5,14.5),
		     labels = c('Year','Prop Infected','Live BA','Cone Density','Live DBH')) +
  theme_bw(base_size = 14) +
  labs(x = 'Effect Size (logit scale)', y = 'Parameter', title = '(a) Occupancy') +
  theme(text = element_text(family="LM Roman 10"), 
	# plot.title = element_text(size = 14), 
	panel.grid = element_blank())

# Detection covariates
# Plot the posterior densities for the covariate effects
alpha.names <- colnames(out$alpha.samples)[-1]
plot.df <- data.frame(val = c(out$alpha.samples[, -1]),
                      parameter = rep(alpha.names, each = out$n.post * out$n.chains))

panel.b <- ggplot(data = plot.df, aes(x = val, group = parameter), color = 'black') +
  geom_vline(xintercept = 0, col = 'black', lty = 2) +
  geom_posterior(draw_ci = TRUE, ci_width = 0.95, draw_sd = FALSE, brighten = FALSE, midline = 'black',
                 color = 'black', alpha = 1, fill = 'grey90') +
  scale_fill_viridis_d() +
  scale_y_continuous(breaks = c(1.5,4.5,7,9.5,12),
                     labels = c('Quadratic Date','Date','Wind','Year','Total Basal Area')) +
  theme_bw(base_size = 14) +
  labs(x = 'Effect Size (logit scale)', y = 'Parameter', title = '(b) Detection') +
  theme(text = element_text(family="LM Roman 10"), 
	#plot.title = element_text(size = 14), 
	panel.grid = element_blank())
plot.fig.3 <- ggarrange(panel.a, panel.b)

ggsave(plot.fig.3, file = 'figures/Figure-3.png', device = 'png', units = 'in', 
       width = 10, height = 5)

# Figure S3 ---------------------------
load('data/spOccupancy-data.rda')
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
plot.fig.s3 <- ggplot(plot.df, aes(x = day, y = mean.val)) +
  geom_ribbon(aes(ymin = low, ymax = high), fill = 'gray') +
  geom_line() +
  theme_bw(base_size = 14) +
  labs(x = 'Date', y = 'Detection Probability') + 
  theme(text = element_text(family = 'LM Roman 10'))
ggsave(plot.fig.s3, file = 'figures/Figure-S3.png', device = 'png', units = 'in', 
       width = 5, height = 5)

# Figure 4 ----------------------------
beta.names <- colnames(out.glmm$beta.samples)[-1]
plot.df <- data.frame(val = c(out.glmm$beta.samples[, -1]),
                      parameter = rep(beta.names, each = out.glmm$n.post * out.glmm$n.chains))

plot.fig.4 <- ggplot(data = plot.df, aes(x = val, group = parameter), color = 'black') +
  geom_vline(xintercept = 0, col = 'black', lty = 2) +
  geom_posterior(draw_ci = TRUE, ci_width = 0.95, draw_sd = FALSE, brighten = FALSE, midline = 'black',
                 color = 'black', alpha = 1, fill = 'grey90') +
  scale_y_continuous(breaks = c(1, 3, 5, 7, 9, 11),
    labels = c('Date','Quadratic Date','Year','Prop Infected','Live DBH','Cone Density')) +
  theme_bw(base_size = 14) +
  labs(x = 'Effect Size', y = 'Parameter') +
  theme(text = element_text(family="LM Roman 10"), 
	panel.grid = element_blank())

ggsave(plot.fig.4, file = 'figures/Figure-4.png', device = 'png', units = 'in', 
       width = 5, height = 5)

# Figure 2 ----------------------------
load("data/spAbundance-data.rda")
# Save the GLMM site names separately
site.names.glmm <- site.names
# Vocalization GLMM (Panel b)
# Year associated with each site
tmp <- data.list$covs$year[, 1]
# Replace indicator with the actual year
tmp <- ifelse(tmp == 0, 2020, ifelse(tmp == 1, 2021, 2022))
plot.df <- data.frame(abund.mean = abund.means, 
		      abund.low = abund.low,
		      abund.high = abund.high,
		      year = factor(tmp),
		      site = site.names.glmm)

plot.df <- arrange(plot.df,(plot.df$abund.mean))
# Get site order for help in making the occupancy plot line up
site.order <- unique(plot.df$site)

plot.df.wide <- plot.df %>%
  pivot_wider(names_from = year, values_from = c(abund.mean, abund.low, abund.high), 
              names_sep = '.')
site.indx <- 1:nrow(plot.df.wide)
plot.df.wide$index.2020 <- site.indx - 0.25
plot.df.wide$index.2021 <- site.indx 
plot.df.wide$index.2022 <- site.indx + 0.25

my.cols <- viridis(3)
names(my.cols) <- c('2020', '2021', '2022')
plot.fig.2.b <- ggplot(data = plot.df.wide) + 
  geom_point(aes(x = index.2020, y = abund.mean.2020, col = '2020'), size = 2) +
  geom_segment(aes(x = index.2020, y = abund.low.2020, xend = index.2020, 
		   yend = abund.high.2020, col = '2020'), lineend = 'butt', size = 1) +
  geom_point(aes(x = index.2021, y = abund.mean.2021, col = '2021'), size = 2) +
  geom_segment(aes(x = index.2021, y = abund.low.2021, xend = index.2021, 
		   yend = abund.high.2021, col = '2021'), lineend = 'butt', size = 1) +
  geom_point(aes(x = index.2022, y = abund.mean.2022, col = '2022'), size = 2) +
  geom_segment(aes(x = index.2022, y = abund.low.2022, xend = index.2022, 
		   yend = abund.high.2022, col = '2022'), lineend = 'butt', size = 1) +
  theme_bw(base_size = 14) + 
  scale_color_manual(values = my.cols) + 
  scale_x_continuous(breaks = site.indx, labels = plot.df.wide$site, name = 'Site') + 
  labs(y = 'Expected number of vocalizations', color = 'Year', title = '(b)') +
  coord_flip() + 
  theme(axis.title.y = element_blank(),
	axis.text.y = element_text(size= 10, color = "black", hjust = 0.5), 
	text = element_text(family = 'LM Roman 10'), 
	legend.position = c(0.75, 0.5),
	panel.grid = element_blank(),
	plot.margin = unit(c(1,1,1,1),"mm"))
# Occpuancy panel
load("data/spOccupancy-data.rda")
site.names.occ <- site.names
psi.means <- apply(out$psi.samples, c(2, 3), mean)
psi.low <- apply(out$psi.samples, c(2, 3), quantile, 0.025)
psi.high <- apply(out$psi.samples, c(2, 3), quantile, 0.975)
# Restrict to only sites in GLMM
site.indx <- rep(NA, length(unique(site.names.glmm)))
for (i in 1:length(site.indx)) {
  site.indx[i] <- which(site.names.occ == unique(site.names.glmm)[i])
}
psi.means <- psi.means[site.indx, ]
psi.low <- psi.low[site.indx, ]
psi.high <- psi.high[site.indx, ]
plot.df <- data.frame(psi.mean = c(psi.means), 
		      psi.low = c(psi.low),
		      psi.high = c(psi.high), 
		      year = rep(c(2020, 2021, 2022), each = nrow(psi.means)),
		      site = rep(site.names.occ[site.indx], times = ncol(psi.means)))
plot.df$site <- factor(plot.df$site, levels = site.order, ordered = TRUE)
plot.df <- plot.df %>%
  arrange(site)

plot.df.wide <- plot.df %>%
  pivot_wider(names_from = year, values_from = c(psi.mean, psi.low, psi.high), 
              names_sep = '.')
site.indx <- 1:nrow(plot.df.wide)
plot.df.wide$index.2020 <- site.indx - 0.25
plot.df.wide$index.2021 <- site.indx 
plot.df.wide$index.2022 <- site.indx + 0.25

my.cols <- viridis(3)
names(my.cols) <- c('2020', '2021', '2022')
plot.fig.2.a <- ggplot(data = plot.df.wide) + 
  geom_point(aes(x = index.2020, y = psi.mean.2020, col = '2020'), size = 2) +
  geom_segment(aes(x = index.2020, y = psi.low.2020, xend = index.2020, 
		   yend = psi.high.2020, col = '2020'), lineend = 'butt', size = 1) +
  geom_point(aes(x = index.2021, y = psi.mean.2021, col = '2021'), size = 2) +
  geom_segment(aes(x = index.2021, y = psi.low.2021, xend = index.2021, 
		   yend = psi.high.2021, col = '2021'), lineend = 'butt', size = 1) +
  geom_point(aes(x = index.2022, y = psi.mean.2022, col = '2022'), size = 2) +
  geom_segment(aes(x = index.2022, y = psi.low.2022, xend = index.2022, 
		   yend = psi.high.2022, col = '2022'), lineend = 'butt', size = 1) +
  theme_bw(base_size = 14) + 
  scale_color_manual(values = my.cols) + 
  scale_x_continuous(breaks = site.indx, labels = plot.df.wide$site, name = 'Site') + 
  labs(y = 'Occupancy probability', color = 'Year', title = '(a)') +
  coord_flip() + 
  guides(color = 'none') +
  theme(axis.title.y = element_blank(),
	axis.text.y = element_text(size= 10, color = "black", hjust = 0.5), 
	text = element_text(family = 'LM Roman 10'), 
	panel.grid = element_blank(),
	plot.margin = unit(c(1,1,1,1),"mm"))
plot.fig.2 <- ggarrange(plot.fig.2.a, plot.fig.2.b, ncol = 2)
ggsave(plot.fig.2, file = 'figures/Figure-2.png', device = 'png', units = 'in', 
       width = 10, height = 7)

# Generate Figure 5 -------------------------------------------------------
# Occupancy portions in 5a
load('data/spOccupancy-data.rda')
n.pred <- 500
occ.pred.vals <- seq(from = min(c(data.list$occ.covs$cones), na.rm = TRUE),
                     to = max(c(data.list$occ.covs$cones), na.rm = TRUE), length.out = n.pred)
X.pred <- matrix(0, n.pred, dim(out$X)[3])
n.post <- out$n.post * out$n.chains
psi.pred <- matrix(NA, n.post, n.pred)
# Intercept
X.pred[, 1] <- 1
# Cones
X.pred[, 3] <- (occ.pred.vals - mean(c(data.list$occ.covs$cones), na.rm = TRUE)) / 
  sd(c(data.list$occ.covs$cones), na.rm = TRUE)
# Setting all other values to 0 (the mean)
for (i in 1:n.post) {
  psi.pred[i, ] <- plogis(X.pred %*% as.matrix(out$beta.samples[i, ]))
}
# Make the plot
plot.df <- data.frame(mean.val = apply(psi.pred, 2, mean),
                      low = apply(psi.pred, 2, quantile, 0.025),
                      high = apply(psi.pred, 2, quantile, 0.975),
		      cones = occ.pred.vals / 0.4)
plot.fig.5.a <- ggplot(plot.df, aes(x = cones, y = mean.val)) +
  geom_ribbon(aes(ymin = low, ymax = high), fill = 'gray') +
  geom_line() +
  theme_bw(base_size = 14) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = 'Cone Density (cones/ha)', y = 'Occupancy Probability', title = '(a)') + 
  theme(text = element_text(family = 'LM Roman 10'))
# Figure 5b
load('data/spOccupancy-data.rda')
n.pred <- 500
occ.pred.vals <- seq(from = min(c(data.list$occ.covs$mean.live.dbh), na.rm = TRUE),
                     to = max(c(data.list$occ.covs$mean.live.dbh), na.rm = TRUE), length.out = n.pred)
X.pred <- matrix(0, n.pred, dim(out$X)[3])
n.post <- out$n.post * out$n.chains
psi.pred <- matrix(NA, n.post, n.pred)
# Intercept
X.pred[, 1] <- 1
# Mean Live DBH 
X.pred[, 5] <- (occ.pred.vals - mean(c(data.list$occ.covs$mean.live.dbh), na.rm = TRUE)) / 
  sd(c(data.list$occ.covs$mean.live.dbh), na.rm = TRUE)
# Setting all other values to 0 (the mean)
for (i in 1:n.post) {
  psi.pred[i, ] <- plogis(X.pred %*% as.matrix(out$beta.samples[i, ]))
}
# Make the plot
plot.df <- data.frame(mean.val = apply(psi.pred, 2, mean),
                      low = apply(psi.pred, 2, quantile, 0.025),
                      high = apply(psi.pred, 2, quantile, 0.975),
		      dbh = occ.pred.vals)
plot.fig.5.b <- ggplot(plot.df, aes(x = dbh, y = mean.val)) +
  geom_ribbon(aes(ymin = low, ymax = high), fill = 'gray') +
  geom_line() +
  theme_bw(base_size = 14) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = 'Mean Live DBH (cm)', y = 'Occupancy Probability', title = '(b)') + 
  theme(text = element_text(family = 'LM Roman 10'))
# Figure 5c
load('data/spAbundance-data.rda')
n.pred <- 500
mu.pred.vals <- seq(from = min(c(data.list$covs$cones), na.rm = TRUE),
                     to = max(c(data.list$covs$cones), na.rm = TRUE), length.out = n.pred)
X.pred <- matrix(0, n.pred, dim(out.glmm$X)[3])
n.post <- out.glmm$n.post * out.glmm$n.chains
mu.pred <- matrix(NA, n.post, n.pred)
# Intercept
X.pred[, 1] <- 1
# Cones
X.pred[, 3] <- (mu.pred.vals - mean(c(data.list$covs$cones), na.rm = TRUE)) / 
  sd(c(data.list$covs$cones), na.rm = TRUE)
# Setting all other values to 0 (the mean)
for (i in 1:n.post) {
  mu.pred[i, ] <- exp(X.pred %*% as.matrix(out.glmm$beta.samples[i, ]))
}
# Make the plot
plot.df <- data.frame(mean.val = apply(mu.pred, 2, mean),
                      low = apply(mu.pred, 2, quantile, 0.025),
                      high = apply(mu.pred, 2, quantile, 0.975),
		      cones = mu.pred.vals / 0.4)
plot.fig.5.c <- ggplot(plot.df, aes(x = cones, y = mean.val)) +
  geom_ribbon(aes(ymin = low, ymax = high), fill = 'gray') +
  geom_line() +
  theme_bw(base_size = 14) +
  labs(x = 'Cone Density (cones/ha)', y = 'Expected Number of Vocalizations', title = '(c)') + 
  theme(text = element_text(family = 'LM Roman 10'))
# Figure 5d
load('data/spAbundance-data.rda')
n.pred <- 500
mu.pred.vals <- seq(from = min(c(data.list$covs$mean.live.dbh), na.rm = TRUE),
                     to = max(c(data.list$covs$mean.live.dbh), na.rm = TRUE), length.out = n.pred)
X.pred <- matrix(0, n.pred, dim(out.glmm$X)[3])
n.post <- out.glmm$n.post * out.glmm$n.chains
mu.pred <- matrix(NA, n.post, n.pred)
# Intercept
X.pred[, 1] <- 1
# Mean Live DBH
X.pred[, 4] <- (mu.pred.vals - mean(c(data.list$covs$mean.live.dbh), na.rm = TRUE)) / 
  sd(c(data.list$covs$mean.live.dbh), na.rm = TRUE)
# Setting all other values to 0 (the mean)
for (i in 1:n.post) {
  mu.pred[i, ] <- exp(X.pred %*% as.matrix(out.glmm$beta.samples[i, ]))
}
# Make the plot
plot.df <- data.frame(mean.val = apply(mu.pred, 2, mean),
                      low = apply(mu.pred, 2, quantile, 0.025),
                      high = apply(mu.pred, 2, quantile, 0.975),
		      dbh = mu.pred.vals)
plot.fig.5.d <- ggplot(plot.df, aes(x = dbh, y = mean.val)) +
  geom_ribbon(aes(ymin = low, ymax = high), fill = 'gray') +
  geom_line() +
  theme_bw(base_size = 14) +
  labs(x = 'Mean Live DBH (cm)', y = 'Expected Number of Vocalizations', title = '(d)') + 
  theme(text = element_text(family = 'LM Roman 10'))
plot.fig.5 <- ggarrange(plot.fig.5.a, plot.fig.5.b, plot.fig.5.c, plot.fig.5.d, nrow = 2, ncol = 2)
ggsave(plot.fig.5, file = 'figures/Figure-5.png', device = 'png', units = 'in', 
       width = 8, height = 8)

# Generate Figure S5 ------------------------------------------------------
plot.df <- data.frame(true = y.site.true.means,
		      fit = y.site.rep.means)
plot.fig.s5 <- ggplot(plot.df, aes(x = true, y = fit)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, lty = 2) + 
  theme_bw(base_size = 14) + 
  labs(x = 'True Vocalization Site Averages', y = 'Estimated Vocalization Site Averages') + 
  theme(text = element_text(family = 'LM Roman 10'))
ggsave(plot.fig.s5, file = 'figures/Figure-S5.png', device = 'png', units = 'in', 
       width = 5, height = 5)

# Generate Figure S6 ------------------------------------------------------
# Figure S6a
load("data/spOccupancy-data.rda")
psi.means <- apply(out$psi.samples, c(2), mean)
coords.sf <- st_as_sf(as.data.frame(data.list$coords), 
		      coords = c('X', 'Y'), 
		      crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs")
coords.sf$psi <- psi.means

plot.fig.s6.a <- ggplot() + 
  geom_sf(data = coords.sf, aes(fill = psi), pch = 21) +
  theme_bw(base_size = 14) + 
  scale_fill_gradientn("Occupancy\nProbability", colors = ocean.tempo(1000), limits = c(0, 1),
                        guide = guide_colourbar(title.position="top", reverse = FALSE), 
  			na.value = NA) +
  labs(x = 'Longitude', y = 'Latitude', title = '(a)') +
  theme(text = element_text(family = 'LM Roman 10'), 
        legend.position = c(0.25, 0.35), 
        legend.background = element_rect(fill = NA))
# Figure S6b
load("data/spAbundance-data.rda")
tmp.df <- data.frame(mu = abund.means, 
		     site = site.names)
mu.means <- tmp.df %>%
  group_by(site) %>%
  summarize(mu = mean(mu)) %>%
  pull(mu)
# Only grab sites with acoustic data.
site.indx <- rep(NA, length(unique(site.names.glmm)))
for (i in 1:length(site.indx)) {
  site.indx[i] <- which(site.names.occ == unique(site.names.glmm)[i])
}
coords.sf.glmm <- coords.sf[site.indx, ]
coords.sf.glmm$mu <- mu.means

plot.fig.s6.b <- ggplot() + 
  geom_sf(data = coords.sf.glmm, aes(fill = mu), pch = 21) +
  theme_bw(base_size = 14) + 
  scale_fill_gradientn("Number of\nVocalizations", colors = ocean.tempo(1000), 
                        guide = guide_colourbar(title.position="top", reverse = FALSE), 
  			na.value = NA) +
  labs(x = 'Longitude', y = 'Latitude', title = '(b)') +
  theme(text = element_text(family = 'LM Roman 10'), 
        legend.position = c(0.25, 0.35), 
        legend.background = element_rect(fill = NA))
plot.fig.s6 <- ggarrange(plot.fig.s6.a, plot.fig.s6.b)
ggsave(plot.fig.s6, file = 'figures/Figure-S6.png', device = 'png', units = 'in', 
       width = 8, height = 5, bg = 'white')
