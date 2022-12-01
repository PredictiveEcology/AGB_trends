##############################################################################################################################################
##############################################################################################################################################
##
## Date Created: Nov 30, 2022
## Auteur: Tyler Rudolph, biologist M.Sc., CFS/NRCAN, Trade, Economics & Industry Branch
##
## Name of script : "01_ABoVE_AGB_time_series_analysis.R"
##
## Description : R script serving to estimate cell-wise linear regression coefficients for ABoVE AGB 31-year time series (1984-2014)
##
##############################################################################################################################################
##############################################################################################################################################

Require::Require(c('terra'))

tile_folders <- list.files('inputs/ABoVE_AGB_30m', pattern='Bh', full.names = T)

slope <- function(x, y = 1:length(x)) {
  if(sum(!is.na(x)) >= 2) {
    return(sum((x - mean(unlist(x), na.rm=T)) * (y - mean(y)), na.rm=T) / sum((x - mean(unlist(x), na.rm=T)) ^ 2, na.rm=T))
  } else {
    return(NA)
  }
}

# slope <- function(x) {
#   if(sum(!is.na(x)) >= 2) {
#     x <- x[!is.na(x)]
#     y <- 1:length(x)
#     return(sum((x - mean(x)) * (y - mean(y))) / sum((x - mean(x)) ^ 2))
#   } else {
#     return(NA)
#   }
# }

i=1

library(terra)
terraOptions(tempdir='scratch', todisc=TRUE)

tile_folder <- tile_folders[i]
ragb <- rast(file.path(tile_folder, list.files(tile_folder, pattern='ragb')))
  
# core 8 = 153
# core 25 = 121
# core 32 = 117

system.time(app(ragb, fun=slope, cores=36,
    filename=file.path(tile_folder, 'slopetest.tif'),
    overwrite=TRUE))

################################################################################
## Estimate cell-wise linear regression coefficient for undisrupted time series
## aka "geographically weighted regression" or GWR



## Group pixel values into (layer-specific) 5-year age classes based on spatiotemporal disturbance dynamics,
## eliminating pixel-years occurring pre-disturbance
# rcl=c(cuts=seq.int(from=0, to=515, by=5)), right=FALSE, others=NA)