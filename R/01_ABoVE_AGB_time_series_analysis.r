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

Require::Require(c('terra','stringr','gdalUtilities'))
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
## 1 a) Estimate cell-wise linear regression coefficients for undisrupted time series
## aka "local" or "geographically weighted regression (GWR)"

system.time({
  
  sapply(1:length(tile_folders), function(i) {
    
    app(rast(file.path(tile_folders[i], list.files(tile_folders[i], pattern='ragb'))), 
        fun=slope, cores=32,
        filename=file.path(tile_folders[i], paste0('agb_slopes_', str_sub(tile_folders[i], start=-7L),'.tif')),
        overwrite=T)
    
  })
  
}) # 3.8 hrs

## 1 b) Combine tiled slope rasters into unified mosaic

## Build virtual raster
gdalbuildvrt(gdalfile = unname(sapply(tile_folders, function(dsn) file.path(dsn, list.files(dsn, pattern='agb_slopes_Bh')))),
             output.vrt = 'cache/AGB_slope_mosaic.vrt')

## Write to raster mosaic
gdalwarp(srcfile = 'cache/AGB_slope_mosaic.vrt', dstfile = 'inputs/ABoVE_AGB_30m/AGB_slope_mosaic.tif')


#################################################################################################################
## 2 a) Calculate cell-specific slopes per 5-year time interval (n=6)

Require::Require(c('terra','stringr'))
terraOptions(tempdir='scratch', todisc=TRUE)

tile_folders <- sort(list.files('inputs/ABoVE_AGB_30m', pattern='Bh', full.names = T))

## slope derivation function
slope <- function(x) {
  y <- which(!is.na(x))
  if(length(y) >= 2) {
    x <- unlist(x[!is.na(x)])
    return(sum((x - mean(x)) * (y - mean(y))) / sum((x - mean(x)) ^ 2))
  } else {
    return(NA)
  }
}

## define time intervals (year ranges between 1984-2014)
timeint <- list(t1=1:5, t2=6:10, t3=11:15, t4=16:20, t5=21:25, t6=26:31)

system.time({
  
  for(i in rev(1:length(tile_folders))) {
    
    cat(i, '/', length(tile_folders), '\n')
    
    tile_folder <- tile_folders[i]
    
    sapply(1:length(timeint), function(timestep) {
      
      ## calculate local slope coefficient for specified time interval
      app(rast(file.path(tile_folders[i], list.files(tile_folders[i], pattern='ragb')))[[timeint[[timestep]]]],
          fun=slope, cores=45,
          filename=file.path(tile_folders[i], paste0('agb_slopes_', names(timeint)[timestep], '_', str_sub(tile_folders[i], start=-7L),'.tif')),
          overwrite=T)
      
    })
    
  }
}) # 10 hours ?

########
## 2 b) Combine tiled slope rasters into numerous unified mosaics

Require::Require('parallel')
no_cores <- 6
cl <- makeCluster(no_cores)
clusterExport(cl, varlist=c('tile_folders'))

parLapply(cl, names(timeint), function(tp) {

  ## Build virtual raster
  Require::Require('gdalUtilities')
  gdalbuildvrt(gdalfile = unname(sapply(tile_folders, function(dsn) file.path(dsn, list.files(dsn, pattern=tp)))),
               output.vrt = paste0('cache/AGB_slope_mosaic_', tp, '.vrt'))
  
  ## Write to raster mosaic
  gdalwarp(srcfile = paste0('cache/AGB_slope_mosaic_', tp, '.vrt'), 
           dstfile = paste0('inputs/ABoVE_AGB_30m/AGB_slope_mosaic_', tp, '.tif'))

    
}) # 23 min

## Visual examination of results
par(mfrow=c(3,2))
for(tp in names(timeint)) {
  # plot(rast(paste0('inputs/ABoVE_AGB_30m/AGB_slope_mosaic_', tp, '.tif')), main=tp)
}

sapply(names(timeint), function(tp) digest(paste0('inputs/ABoVE_AGB_30m/AGB_slope_mosaic_', tp, '.tif'), algo='xxhash64'))
sapply(names(timeint), function(tp) file.size(paste0('inputs/ABoVE_AGB_30m/AGB_slope_mosaic_', tp, '.tif')))

system.time(test <- diff(rast(sapply(names(timeint), function(tp) rast(paste0('inputs/ABoVE_AGB_30m/AGB_slope_mosaic_', tp, '.tif')))),
            filename='cache/diff.tif'))

#########################################################################################
## 3) 


# classify(rast(file.path(tile_folders[i], list.files(tile_folders[i], pattern='rage'))),
#          rcl=c(cuts=seq.int(from=0, to=515, by=5)), right=FALSE, others=NA)



