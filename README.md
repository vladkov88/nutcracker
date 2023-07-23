# Paired acoustic recordings and point count surveys reveal Clark’s Nutcracker and Whitebark Pine associations across a protected area

### In Review

### Vladimir Kovalenko, Jeffrey W. Doser, Lisa J. Bate, Diana L. Six

### Please contact the first author for questions about the code or data used in the manuscript: vladimir.kovalenko@umconnect.umt.edu

---------------------------------

## Abstract

Global declines of tree species have led to dramatic shifts in forest ecosystem composition, biodiversity, and functioning. These changes have consequences for both the forest plant and wildlife communities, particularly when declining species are parties to coevolved mutualisms. Whitebark pine (Pinus albicaulis) is a rapidly declining keystone species in high-elevation ecosystems and an obligate mutualist of Clark’s Nutcracker, an avian seed predator and disperser. Using passive acoustic monitoring paired with traditional point count surveys, we investigated how stand characteristics of whitebark pine in a protected area (Glacier National Park, Montana, USA) influence occurrence and activity patterns in Clark’s Nutcracker. Through the use of Bayesian spatial occupancy models and generalized linear mixed models, we found that habitat use of Clark’s Nutcracker is primarily influenced by cone production and increasing diameter of live whitebark pine. Additionally, we demonstrate the effectiveness of combining traditional point count surveys with passive acoustic monitoring to assess habitat use and its drivers for highly mobile species. These findings emphasize the importance of restoration strategies that work to protect existing large cone-bearing whitebark pines to preserve natural tree regeneration potential as well as foraging resources for Clark’s Nutcracker.   

## Repository Directory

+ `main-spatial-occ.R`: script to fit the spatial multi-season occupancy model to assess Clark's Nutcracker use of whitebark pine forests.
+ `main-glmm.R`: script to fit the generalized linear mixed model to the acoustic vocalization data to assess Clark's Nutcracker use of whitebark pine forests.
+ `acoustic-covariates-1.csv`: some of the covariates used for the acoustic data GLMM.
+ `acoustic-covariates-2.csv`: some of the covariates used for the acoustic data GLMM.
+ `acoustic-det-covs.csv`: covariates related to variation in detection probability of vocalizations for the acoustic GLMM.
+ `acoustic-vocalization-data.csv`: the vocalization data for use in the acoustic GLMM that are derived from BirdNet.
+ `coordinates.csv`: spatial coordinates of the sites.
+ `covariate-point-count-data.csv`: covariates used for the point count occupancy model.
+ `point-count-data.csv`: the point count data for use in the occupancy model.
+ `site-indx-acoustic.csv`: index to indicate which sites have acoustic data associated with them.
+ `year-acoustic.csv`: indicates the years of sampling for each site in the acoustic data set.
+ `year-point-count.csv`: indicates the years of sampling for all sites in the count data set.
