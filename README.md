# Paired acoustic recordings and point count surveys reveal Clark’s Nutcracker and Whitebark Pine associations across Glacier National Park

### In Review

### Vladimir Kovalenko, Jeffrey W. Doser, Lisa J. Bate, Diana L. Six

### Please contact the first author for questions about the code or data used in the manuscript: vladimir.kovalenko@umconnect.umt.edu

---------------------------------

## Abstract

Global declines in tree populations have led to dramatic shifts in forest ecosystem composition, biodiversity, and functioning. These changes have consequences for both forest plant and wildlife communities, particularly when declining species are involved in coevolved mutualisms. Whitebark pine (Pinus albicaulis) is a rapidly declining keystone species in high-elevation ecosystems and an obligate mutualist of Clark’s Nutcracker (Nucifraga columbiana), an avian seed predator and disperser. By augmenting traditional point count surveys with passive acoustic monitoring, we investigated how stand characteristics of whitebark pine in a protected area (Glacier National Park, Montana, USA) influenced occupancy and vocal activity patterns in Clark’s Nutcracker. Using Bayesian spatial occupancy models and generalized linear mixed models, we found that habitat use of Clark’s Nutcracker was primarily supported by greater cone production and increasing diameter of live whitebark pine. Additionally, we demonstrated the value of performing parallel analyses with traditional point count surveys and passive acoustic monitoring to provide multiple lines of evidence for relationships between Clark’s Nutcracker and whitebark pine forest characteristics. These findings emphasize the importance of management strategies that work to protect existing large cone-bearing whitebark pines to preserve natural tree regeneration potential as well as foraging resources for Clark’s Nutcracker.   

## Repository Directory

### [code/](./code)

+ `main-spatial-occ.R`: script to fit the spatial multi-season occupancy model to assess Clark's Nutcracker use of whitebark pine forests.
+ `main-glmm.R`: script to fit the generalized linear mixed model to the acoustic vocalization data to assess Clark's Nutcracker use of whitebark pine forests.
+ `summary.R`: script to summarize the results from both models and generate all figures used in the manuscript. 

### [data/](.data)

+ `acoustic-covariates-1.csv`: some of the covariates used for the acoustic data GLMM.
+ `acoustic-covariates-2.csv`: some of the covariates used for the acoustic data GLMM.
+ `acoustic-det-covs.csv`: covariates related to variation in detection probability of vocalizations for the acoustic GLMM.
+ `acoustic-vocalization-data.csv`: the vocalization data for use in the acoustic GLMM that are derived from BirdNet.
+ `coordinates.csv`: spatial coordinates of the sites.
+ `covariate-point-count-data.csv`: covariates used for the point count occupancy model.
+ `hr_cones.csv`: home range estimates discussed in the Supplemental Material and used to generate an informed prior in the spatial occupancy model.
+ `point-count-data.csv`: the point count data for use in the occupancy model.
+ `site-indx-acoustic.csv`: index to indicate which sites have acoustic data associated with them.
+ `spAbundance-data.rda`: the acoustic vocalization data formatted for fitting a GLMM using `spAbundance`.
+ `spOccupancy-data.rda`: the point count data formatted for fitting a spatial multi-season occupancy model using `spOccupancy`.
+ `year-acoustic.csv`: indicates the years of sampling for each site in the acoustic data set.
+ `year-point-count.csv`: indicates the years of sampling for all sites in the count data set.

## [figures/](.figures)

Contains all figures included in the manuscript.

## [results/](.results)

+ `glmm-vocalization-results.rda`: results from the GLMM fit in `spAbundance` using the acoustic data.
+ `spatial-occ-results.rda`: results from the spatial occupancy model fit using `spOccupancy`.
