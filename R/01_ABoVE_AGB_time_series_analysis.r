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

Require::Require(c('terra','stringr'))
terraOptions(tempdir='scratch', todisc=TRUE)

tile_folders <- sort(list.files('inputs/ABoVE_AGB_30m', pattern='Bh', full.names = T))

slope <- function(x) {
  y <- which(!is.na(x))
  if(length(y) >= 2) {
    x <- unlist(x[!is.na(x)])
    return(sum((x - mean(x)) * (y - mean(y))) / sum((x - mean(x)) ^ 2))
  } else {
    return(NA)
  }
}

################################################################################
## 1) Estimate cell-wise linear regression coefficients for undisrupted time series
## aka "local" or "geographically weighted regression (GWR)"

system.time({
  
  sapply(1:length(tile_folders), function(i) {
    
    app(rast(file.path(tile_folders[i], list.files(tile_folders[i], pattern='ragb'))), 
        fun=slope, cores=32,
        filename=file.path(tile_folders[i], paste0('agb_slopes_', str_sub(tile_folders[i], start=-7L),'.tif')),
        overwrite=T)
    
  })
  
}) # 3.8 hrs


############################################################################################################
## 2) Group pixel values into (layer-specific) 5-year age classes based on spatiotemporal disturbance dynamics
##    and calculate slopes per 5-year period

Require::Require(c('terra','stringr'))
terraOptions(tempdir='scratch', todisc=TRUE)

tile_folders <- sort(list.files('inputs/ABoVE_AGB_30m', pattern='Bh', full.names = T))

slope <- function(x) {
  y <- which(!is.na(x))
  if(length(y) >= 2) {
    x <- unlist(x[!is.na(x)])
    return(sum((x - mean(x)) * (y - mean(y))) / sum((x - mean(x)) ^ 2))
  } else {
    return(NA)
  }
}

classify(rage, rcl=c(cuts=seq.int(from=0, to=515, by=5)), right=FALSE, others=NA)

