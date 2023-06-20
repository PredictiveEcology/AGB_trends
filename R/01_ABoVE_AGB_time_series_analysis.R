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

Require::Require(c('terra','stringr','gdalUtilities','dplyr','sf','reproducible'), purge = T)

file.remove(file.path('/mnt/scratch/trudolph/AGB_trends/terra', list.files('/mnt/scratch/trudolph/AGB_trends/terra')))

oldTmpDir <- tempdir()
newTmpDir <- file.path("/mnt/scratch/trudolph/AGB_trends/tmp")
if (!dir.exists(newTmpDir)) dir.create(newTmpDir, recursive = TRUE)
newTmpDir <- tools::file_path_as_absolute(newTmpDir)
Sys.setenv(TMPDIR = newTmpDir)
unlink(oldTmpDir, recursive = TRUE)
tempdir(check = TRUE)

terraOptions(tempdir = '/mnt/scratch/trudolph/AGB_trends/terra',
             memmax = 25,
             memfrac = 0.8,
             progress = 1,
             verbose = TRUE)

tile_folders <- sort(list.files('inputs/clean/tiled', pattern='Bh', full.names = TRUE))

################################################################################
## 1.1) Estimate cell-wise linear regression coefficients for undisrupted time series
## aka "local" or "geographically weighted regression (GWR)"

sapply(1:length(tile_folders), function(i) {
  ## 1.1.1) Derive slope of numerical vector across a time series
  terraOptions(datatype = 'FLT4S')
  app(rast(file.path(tile_folders[i], list.files(tile_folders[i], pattern = 'ragb'))),
      fun = function(x, ff) ff(x), cores = 32, ff = slope,
      filename = file.path(tile_folders[i], paste0('agb_slopes_', str_sub(tile_folders[i], start = -7L), '.tif')),
      overwrite = TRUE)

  ## 1.1.2) stock number of non-NA values for subsequent weighted standard deviation
  terraOptions(datatype = 'INT1U')
  app(rast(file.path(tile_folders[i], list.files(tile_folders[i], pattern = 'ragb'))),
      fun = function(x, ff) ff(x), cores = 10, ff = nsamp,
      filename = file.path(tile_folders[i], paste0('agb_sample_size_', str_sub(tile_folders[i], start = -7L),'.tif')),
      overwrite = TRUE)

  return(invisible(NULL))
}) # ~ 7 hrs

################################################################################
## 1.2) Combine tiled slope rasters into unified mosaics

## 1.2.1) Build virtual rasters
gdalbuildvrt(gdalfile = unname(sapply(tile_folders, function(dsn) file.path(dsn, list.files(dsn, pattern='agb_slopes_Bh')))),
             output.vrt = 'cache/AGB_slope_mosaic.vrt')
gdalbuildvrt(gdalfile = unname(sapply(tile_folders, function(dsn) file.path(dsn, list.files(dsn, pattern='agb_sample_size_Bh')))),
             output.vrt = 'cache/AGB_sampleSize_mosaic.vrt')

## 1.2.2) Write to raster mosaics
gdalwarp(srcfile = 'cache/AGB_slope_mosaic.vrt', dstfile = 'outputs/AGB_slope_mosaic.tif', overwrite=T)
gdalwarp(srcfile = 'cache/AGB_sampleSize_mosaic.vrt', dstfile = 'outputs/AGB_sample_size_mosaic.tif', overwrite=T)

#################################################################################################################
## 2.1) Calculate cell-specific slopes per 5-year time interval (n=6)

## define time intervals (year ranges between 1984-2014)
timeint <- list(t1=1:5, t2=6:10, t3=11:15, t4=16:20, t5=21:25, t6=26:31)

for(i in 1:length(tile_folders)) {

  cat(i, '/', length(tile_folders), '\n')

  tile_folder <- tile_folders[i]

  sapply(1:length(timeint), function(timestep) {

    ## 2.1.1) calculate local slope coefficient for specified time interval
    terraOptions(datatype = 'FLT4S')
    app(rast(file.path(tile_folders[i], list.files(tile_folders[i], pattern='ragb')))[[timeint[[timestep]]]],
        function(x, ff) ff(x), cores=32, ff=slope,
        filename = file.path(tile_folders[i], paste0('agb_slopes_', names(timeint)[timestep], '_', str_sub(tile_folders[i], start=-7L),'.tif')),
        overwrite = T)

    ## 2.1.2) stock number of non-NA values for subsequent weighted standard deviation
    terraOptions(datatype = 'INT1U')
    app(rast(file.path(tile_folders[i], list.files(tile_folders[i], pattern='ragb')))[[timeint[[timestep]]]],
        fun = function(x, ff) ff(x), cores = 10, ff = nsamp,
        filename = file.path(tile_folders[i], paste0('agb_sample_size_', names(timeint)[timestep], '_', str_sub(tile_folders[i], start=-7L),'.tif')),
        overwrite = T)

    return(invisible(NULL))

  })

} # 16 hours ?

########
## 2.2) Combine tiled slope rasters into numerous unified mosaics
Require::Require('parallel')
no_cores <- 6 ## TODO: why 6? length(timeint)?
cl <- parallelly::makeClusterPSOCK(no_cores,
                                   default_packages = c("terra","gdalUtilities","stringr"),
                                   rscript_libs = .libPaths(),
                                   autoStop = TRUE)
clusterExport(cl, varlist=c('tile_folders'))
parallel::clusterEvalQ(cl, {
  terraOptions(tempdir = '/mnt/scratch/trudolph/AGB_trends/terra',
               memmax = 25,
               memfrac = 0.6,
               progress = 1,
               verbose = TRUE)
})

parLapply(cl, names(timeint), function(tp) {

  ## 2.2.1) Build virtual rasters
  flist <- unname(sapply(tile_folders, function(dsn) file.path(dsn, list.files(dsn, pattern=tp))))

  gdalbuildvrt(gdalfile = flist[str_detect(flist, 'slope')],
               output.vrt = paste0(terraOptions(print=F)$tempdir, '/AGB_slope_mosaic_', tp, '.vrt'))

  gdalbuildvrt(gdalfile = flist[str_detect(flist, 'sample_size')],
               output.vrt = paste0(terraOptions(print=F)$tempdir, '/AGB_sample_size_mosaic_', tp, '.vrt'))

  ## 2.2.2) Write to raster mosaics
  gdalwarp(srcfile = paste0(terraOptions(print=F)$tempdir, '/AGB_slope_mosaic_', tp, '.vrt'),
           dstfile = paste0('outputs/AGB_slope_mosaic_', tp, '.tif'))

  gdalwarp(srcfile = paste0(terraOptions(print=F)$tempdir, '/AGB_sample_size_mosaic_', tp, '.vrt'),
           dstfile = paste0('outputs/AGB_sample_size_mosaic_', tp, '.tif'))

  return(invisible(NULL))

}) # 23 min

stopCluster(cl)

# #########################################################################################
# ## Visual examination of results !! finesse and write to png !!
#
# par(mfrow=c(3,2))
# for(tp in names(timeint)) {
#   plot(rast(paste0('outputs/AGB_slope_mosaic_', tp, '.tif')), main=tp)
# }
#
# sapply(names(timeint), function(tp) digest::digest(paste0('outputs/AGB_slope_mosaic_', tp, '.tif'), algo='xxhash64'))
# sapply(names(timeint), function(tp) file.size(paste0('outputs/AGB_slope_mosaic_', tp, '.tif')))
#
# ## write to PNG !!
# par(mfrow=c(2,3))
# hist(rast(paste0('outputs/AGB_slope_mosaic_t1.tif')), main='t1')
# hist(rast(paste0('outputs/AGB_slope_mosaic_t2.tif')), main='t2')
# hist(rast(paste0('outputs/AGB_slope_mosaic_t3.tif')), main='t3')
# hist(rast(paste0('outputs/AGB_slope_mosaic_t4.tif')), main='t4')
# hist(rast(paste0('outputs/AGB_slope_mosaic_t5.tif')), main='t5')
# hist(rast(paste0('outputs/AGB_slope_mosaic_t6.tif')), main='t6')

############################################################################################################
## 3) Group slopes by age at time x (band argument determines reference layer/year),
##    effectively masking out pixels disturbed mid-time series
no_cores <- 6
cl <- parallelly::makeClusterPSOCK(no_cores,
                                   default_packages = c("terra","gdalUtilities"),
                                   rscript_libs = .libPaths(),
                                   autoStop = TRUE)

parallel::clusterExport(cl, varlist = c("tile_folders","no_cores"), envir = environment())
parallel::clusterEvalQ(cl, {
  terraOptions(tempdir = '/mnt/scratch/trudolph/AGB_trends/terra',
               memmax = 25,
               memfrac = 0.6 / no_cores,
               progress = 1,
               verbose = T)
})

parallel::parLapply(cl, 1:6, function(i) {

  ## 3.1) Build virtual raster of estimated stand age at time 0 (i.e. 1984)
  gdalbuildvrt(
    gdalfile = sapply(tile_folders, function(dsn) unname(file.path(dsn, list.files(dsn, pattern='rage')))),
    output.vrt = paste0('cache/AGB_age_mosaic_t', i, '.vrt'),
    ## !! MODIFY if doing alternative time steps is a DESIRED functionality (e.g. 10-year, etc.)
    b=c(1,6,11,16,21,26)[i], # time 1 = 1984 etc.
    overwrite=T)

  ## 3.2) Write to raster mosaic (stand age, kNN 2020)
  gdalwarp(srcfile = paste0('cache/AGB_age_mosaic_t', i, '.vrt'),
           dstfile = paste0('inputs/clean/AGB_age_mosaic_t', i, '.tif'),
           overwrite = T)

  ## 3.3) Group into 5 discrete age classes
  ageRast <- classify(rast(paste0('inputs/clean/AGB_age_mosaic_t', i, '.tif')),
                      rcl=cbind(from=c(0, 25, 50, 80, 125), to=c(25, 50, 80, 125, 500), becomes=1L:5L),
                      right=FALSE, others=NA_integer_)

  names(ageRast) <- 'ageClass'
  levels(ageRast) <- data.frame(value=1:5, ageClass=paste0(c(0, 25, 50, 80, 125), '-', to=c(25, 50, 80, 125, 500)))
  writeRaster(ageRast, paste0('inputs/clean/AGB_age_mosaic_classes_t', i, '.tif'), overwrite=T)

  return(invisible(NULL))

})

parallel::stopCluster(cl)

###################################################################################
## 4) rasterize study area by categorical zones of interest (basis of subsequent results comparison)
##    WBI and ecozones by default

# targetCRS <- paste0("PROJCRS[\"Canada_Albers_Equal_Area_Conic\",\n",
#        "    BASEGEOGCRS[\"NAD83\",\n",
#        "        DATUM[\"North American Datum 1983\",\n",
#        "            ELLIPSOID[\"GRS 1980\",6378137,298.257222101004,\n",
#        "                LENGTHUNIT[\"metre\",1]]],\n",
#        "        PRIMEM[\"Greenwich\",0,\n",
#        "            ANGLEUNIT[\"degree\",0.0174532925199433]],\n",
#        "        ID[\"EPSG\",4269]],\n    CONVERSION[\"unnamed\",\n",
#        "        METHOD[\"Albers Equal Area\",\n",
#        "            ID[\"EPSG\",9822]],\n",
#        "        PARAMETER[\"Latitude of false origin\",40,\n",
#        "            ANGLEUNIT[\"degree\",0.0174532925199433],\n",
#        "            ID[\"EPSG\",8821]],\n",
#        "        PARAMETER[\"Longitude of false origin\",-96,\n",
#        "            ANGLEUNIT[\"degree\",0.0174532925199433],\n",
#        "            ID[\"EPSG\",8822]],\n",
#        "        PARAMETER[\"Latitude of 1st standard parallel\",50,\n",
#        "            ANGLEUNIT[\"degree\",0.0174532925199433],\n",
#        "            ID[\"EPSG\",8823]],\n",
#        "        PARAMETER[\"Latitude of 2nd standard parallel\",70,\n",
#        "            ANGLEUNIT[\"degree\",0.0174532925199433],\n",
#        "            ID[\"EPSG\",8824]],\n",
#        "        PARAMETER[\"Easting at false origin\",0,\n",
#        "            LENGTHUNIT[\"metre\",1],\n",
#        "            ID[\"EPSG\",8826]],\n",
#        "        PARAMETER[\"Northing at false origin\",0,\n",
#        "            LENGTHUNIT[\"metre\",1],\n",
#        "            ID[\"EPSG\",8827]]],\n",
#        "    CS[Cartesian,2],\n",
#        "        AXIS[\"easting\",east,\n",
#        "            ORDER[1],\n",
#        "            LENGTHUNIT[\"metre\",1,\n",
#        "                ID[\"EPSG\",9001]]],\n",
#        "        AXIS[\"northing\",north,\n",
#        "            ORDER[2],\n",
#        "            LENGTHUNIT[\"metre\",1,\n",
#        "                ID[\"EPSG\",9001]]]]")
#
# source('modules/AGB_dataPrep/R/analysisZones.R')
# zoi <- createAnalysisZones(st_read('outputs/WBI_studyArea.gpkg'), targetCRS, "inputs")
# st_write(zoi, dsn='outputs/WBI_studyArea.gpkg', delete_layer=T)

no_cores <- 6
cl <- parallelly::makeClusterPSOCK(no_cores,
                                   default_packages = c("terra","gdalUtilities"),
                                   rscript_libs = .libPaths(),
                                   autoStop = TRUE)
parallel::clusterExport(cl, varlist = c('prepZones', 'no_cores'), envir = environment())
parallel::clusterEvalQ(cl, terraOptions(tempdir = '/mnt/scratch/trudolph/AGB_trends/terra', memfrac = 0.5 / no_cores))

parallel::parLapply(cl, 1:6, function(i) {

  prepZones(#zoi = zoi,
    field = 'ECOZONE', ## CAN ALTER FOR ECOREGION OR ECOPROVINCE, (which I did so results on file)
    ageClass = rast(paste0('inputs/clean/AGB_age_mosaic_classes_t', i, '.tif')),
    file.id = paste0('WBI_ecozone_t', i),
    ow = T)

  return(invisible(NULL))

}) # 17 min

parallel::stopCluster(cl)

####################################################################################################
## 5) Calculate comparative summary statistics by categorical zone of interest for all time periods

## Note in following that age at beginning of the 31 year time series (1984-2014) is identical to age at beginning of 't1' time interval (i.e. 1984-1988)
irast <- list(
  slope = list.files('outputs', pattern = "slope_mosaic", full.names = TRUE),
  w = list.files('outputs', pattern = "sample_size", full.names = TRUE),
  ## these last values '2' (below) refer to ageClass at 'time 0' (i.e. 1984) used
  ## for the complete time series slope raster mosaic stats assessment.
  ## this should be '1', but currently set to 6 years in b/c disturbed pixels are only known as of 1987
  ecozone = list.files('outputs', pattern = 'ZOIxageClass_WBI_ecozone', full.names = TRUE)[-c(2,4,6,8,10,12)][c(1:6,2)],
  ecoregion = list.files('outputs', pattern = 'ZOIxageClass_WBI_ecoregion', full.names = TRUE)[-c(2,4,6,8,10,12)][c(1:6,2)],
  ecoprovince = list.files('outputs', pattern = 'ZOIxageClass_WBI_ecoprovince', full.names = TRUE)[-c(2,4,6,8,10,12)][c(1:6,2)]
)

ncores = 7 ## TODO: length(timeint) + 1
cl <- parallelly::makeClusterPSOCK(ncores, default_packages = c('terra','gdalUtilities','dplyr'),
                                   rscript_libs = .libPaths(), autoStop = TRUE)
parallel::clusterExport(cl, varlist = c('zoneStats','irast','ncores'), envir = environment())
parallel::clusterEvalQ(cl, terraOptions(tempdir = '/mnt/scratch/trudolph/AGB_trends/terra', memfrac = 0.5 / ncores))

system.time({
  parallel::parLapply(cl, 1:7, function(i, svar = 'ecozone', maskRaster = NULL) {
    ## TODO: use maskRaster file name to qualify file.id writeRaster tag
    if (i == 7) {
      file.id <- paste0('WBI_', svar)
      # file.id <- paste0('WBI_distMask_', svar)
    } else {
      file.id <- paste0('WBI_', svar, '_t', i)
      # file.id <- paste0('WBI_distMask_', svar, '_t', i)
    }

    zoneStats(slopeRaster = rast(irast$slope[i]),
              weightRaster = rast(irast$w[i]),
              zoneRaster = rast(irast[[svar]][i]),
              ## maskRaster arg can be either e.g. was it disturbed? or e.g. is it forested? or both....
              maskRaster = maskRaster, # e.g. use rast('inputs/raw/ABoVE_ForestDisturbance_Agents/binary_disturbed_mosaic.tif') for pixels disturbed over course of time series (according to ABoVE)
              file.id = file.id)

    return(invisible(NULL))
  })
}) # ~ 2 hrs

parallel::stopCluster(cl)


###################################################################################
## 6) Plot differences

#####
## 6 a) without disturbance mask

## i=1 corresponds to 31-year time series, i=2 corresponds to time interval t1 (1984-1988), and so on and so forth
plotZoneStats(file2plot = paste0('outputs/zoneStats_summary_WBI_ecozone.rds'))

## x Ecozone x ageClass
pdf(file = paste0('outputs/AGB_temporal_trends_x_ECOZONE_x_ageClass_', Sys.Date(), '.pdf'))
lapply(plotZoneStatsIntervals(files2plot = file.path('outputs', list.files('outputs/', pattern = 'zoneStats_summary_WBI_ecozone_')),
                              weighted = T, xVar = 'tp', groupVar = 'ageClass', ptype = 2), plot)
dev.off()

## x ageClass x Ecozone (THIS ONE CRASHES ATM, not sure why yet)
pdf(file = paste0('outputs/AGB_temporal_trends_x_ageClass_x_ECOZONE_', Sys.Date(), '.pdf'))
lapply(plotZoneStatsIntervals(files2plot = file.path('outputs', list.files('outputs/', pattern = 'zoneStats_summary_WBI_ecozone_')),
                              weighted = T, xVar = 'tp', catVar = 'ageClass', groupVar='ECOZONE', ptype = 2), plot)
dev.off()

#########
## 6 b) with disturbance mask
plotZoneStats(file2plot = paste0('outputs/zoneStats_summary_WBI_distMask_ecozone.rds'))

## x ageClass x Ecozone
pdf(file = paste0('outputs/AGB_temporal_trends_x_ECOZONE_distMask_', Sys.Date(), '.pdf'))
lapply(plotZoneStatsIntervals(files2plot = file.path('outputs', list.files('outputs/', pattern = 'zoneStats_summary_WBI_distMask_ecozone_')),
                              weighted = T, xVar = 'tp', catVar = 'ageClass', groupVar='ECOZONE', ptype = 2), plot)
dev.off()

gp <- plotZoneStatsIntervals(files2plot = file.path('outputs', list.files('outputs/', pattern='WBI_distMask_ecozone')))

############################################################
## 7) Test for significant differences between groups




######################
## Functions

## Derive slope of numerical vector across a time series
slope <- function(x) {
  y <- which(!is.na(x))
  if (length(y) >= 2) {
    x <- unlist(x[!is.na(x)])
    if (length(unique(x)) == 1) {
      return(0)
    } else {
      return(sum((x - mean(x)) * (y - mean(y))) / sum((x - mean(x)) ^ 2))
    }
  } else {
    return(NA)
  }
}

## stock number of non-NA values for subsequent weighted standard deviation
nsamp <- function(x) {
  x <- sum(!is.na(x))
  if (x > 1) {
    return(x)
  } else {
    return(NA)
  }
}

############################################################################
## Create unique categorical zones for comparative analysis
## (e.g. study area ZOI x age class for desired time period)

prepZones <- function(zoi = st_read('outputs/WBI_studyArea.gpkg', quiet=T),
                      field = 'ECOZONE',
                      ageClass = rast('inputs/clean/AGB_age_mosaic_classes_t1.tif'),
                      file.id = 'WBI_ecozone_t1',
                      cropRaster = NULL, #rast('inputs/clean/tiled/Bh11v06/agb_slopes_Bh11v06.tif'),
                      ow = T) {

  Require::Require(c('sf','terra','dplyr','stringr'))

  if(!is.null(cropRaster)) {
    ageClass <- crop(ageClass, cropRaster)
  }

  ## QC of field to rasterize
  # if(!inherits(class(pull(zoi, field)), 'numeric')) zoi <- mutate(zoi, !!field := as.integer(factor(!!as.name(field))))

  ## define comparative zones of interest (ZOI x age class for specified time period) within study area
  rasterize(x = zoi %>% vect(.),
            y = ageClass,
            field = field,
            filename = paste0('outputs/ZOI_', file.id, '.tif'),
            overwrite = ow)

  # safeguard
  maxchar <- max(nchar(levels(rast(paste0('outputs/ZOI_', file.id, '.tif')))[[1]]$value))
  rcoef <- as.numeric(str_c(c(1, rep(0, maxchar)), collapse=""))

  ## combine study area zoi and age class mosaic into unique combined raster categories (when age class from 100 - 500 [n=5])
  zoneCat <- rast(paste0('outputs/ZOI_', file.id, '.tif'))
  system.time(jointRast <- as.int(rcoef * ageClass + zoneCat))

  ## re-assign factor levels
  ftab <- data.frame(value = unique(jointRast)) %>%
    rename(value = ageClass) %>%
    mutate(ageClass = levels(ageClass)[[1]]$ageClass[match(as.integer(str_sub(value, start=1L, end=1L)), levels(ageClass)[[1]]$value)],
           !!field := levels(zoneCat)[[1]][,field][match(as.integer(str_sub(value, start=-maxchar)), levels(zoneCat)[[1]]$value)])

  ftab <- bind_cols(ftab, zoneCat = paste0(ftab[,3], ' (', ftab[,2], ' yrs)')) %>%
    relocate(zoneCat, .after = value)

  levels(jointRast) <- ftab

  writeRaster(jointRast, filename = paste0('outputs/ZOIxageClass_', file.id, '.tif'), overwrite = ow)

  return(invisible(NULL))

}

#####################################################################
## Calculate summary statistics by categorical zone of interest
## (weighted mean & sd of slopes by age class, time period & study area zones of interest)
## - WBI and ecozones by default for every time period

zoneStats <- function(slopeRaster = rast('outputs/AGB_slope_mosaic.tif'),
                      weightRaster = rast('outputs/AGB_sample_size_mosaic.tif'),
                      zoneRaster = rast('outputs/ZOIxageClass_WBI_ecozone_t1.tif'),
                      cropRaster = NULL, #rast('inputs/clean/tiled/Bh11v06/agb_slopes_Bh11v06.tif'),
                      maskRaster = NULL,
                      file.id = 'WBI_ecozones') {

  names(weightRaster) <- 'w'
  names(slopeRaster) <- 'slope'

  ## crop to smaller (e.g. test) area, if one provided
  if(!is.null(cropRaster)) {
    slopeRaster <- crop(slopeRaster, cropRaster)
    weightRaster <- crop(weightRaster, cropRaster)
    zoneRaster <- crop(zoneRaster, cropRaster)
    if(!is.null(maskRaster)) maskRaster <- crop(maskRaster, cropRaster)
  }

  ## apply mask, if provided
  if(!is.null(maskRaster)) {
    slopeRaster <- mask(slopeRaster, maskRaster)
    weightRaster <- mask(weightRaster, maskRaster)
    zoneRaster <- mask(zoneRaster, maskRaster)
  } # 13 min

  a <- cats(zoneRaster)[[1]]
  levels(zoneRaster) <- NULL
  names(zoneRaster) <- 'value'

  ## 1) count of slope observations by zone
  a <- a %>%
    left_join(., zonal(slopeRaster, zoneRaster, fun='notNA'), by='value') %>%
    rename(count = slope)

  if(any(!is.na(a$count))) {

    ## 2) geometric (unweighted) mean by zone
    meanRast <- zonal(slopeRaster, zoneRaster, fun='mean', as.raster=T, na.rm=T) # 12 min

    a <- a %>%
      left_join(., zonal(meanRast, zoneRaster, fun='min', na.rm=T), by='value') %>%
      rename(mean.slope = slope)

    ## 3) sd by zone
    a <- a %>%
      left_join(., zonal(x = (slopeRaster - meanRast)^2,
                         z = zoneRaster, fun='sum', na.rm=T), by='value') %>%
      mutate(sd = sqrt(slope / count)) %>%
      select(!slope)

    ## 4.1) sum of wi*xi by zone (numerator of weighted mean)
    b <- zonal(slopeRaster * weightRaster, zoneRaster, fun='sum', na.rm=T) %>%
      rename(b = slope) # ~ 11 min

    ## 4.2) sum of weights by zone (denominator of weighted mean)
    w <- zonal(weightRaster, zoneRaster, fun = 'sum', na.rm = T) # ~ 3.7 min

    ## 4.3) derive weighted mean
    a <- a %>%
      left_join(., b, by = 'value') %>%
      left_join(., w, by = 'value') %>%
      mutate(wtd.mean.slope = b / w, .before = w)

    rm(w)

    ## 5) derive weighted sd
    a <- a %>%
      left_join(.,
                zonal(x = weightRaster * (slopeRaster - classify(zoneRaster, rcl=a[,c('value','wtd.mean.slope')])) ^ 2,
                      z = zoneRaster, fun = 'sum', na.rm = T) %>% rename(x = w), by='value') %>%
      mutate(wtd.sd = sqrt(x / w)) %>%
      filter(!is.na(count)) %>%
      select(!c(b, w, x))

    ##  6 b) write to file
    saveRDS(a, file = paste0('outputs/zoneStats_summary_', file.id, '.rds'))

  }

  return(invisible(NULL))

}

###################################################################
## Plot summary statistics

plotZoneStats <- function(file2plot = 'outputs/zoneStats_summary_WBI_ecozone.rds',
                          weighted = T, out = 1) {

  Require::Require(c('ggplot2','stringr','dplyr'))

  sumtab <- readRDS(file2plot) %>%
    mutate(ageClass = factor(ageClass, levels=unique(ageClass), ordered=T))
  catvar <- names(sumtab)[4]

  # y-axis label corresponding to examined time period
  tref <- cbind.data.frame(tp = paste0('t', 1:6),
                           yr1 = c(1984, 1989, 1994, 1999, 2004, 2009),
                           yr2 = c(1988, 1993, 1998, 2003, 2008, 2014))

  tp <- ifelse(str_detect(file2plot, '_t'),
               paste0('(', str_c(filter(tref, tp == str_sub(file2plot, start=-6L, end=-5L)) %>% select(yr1, yr2) %>% unlist(), collapse = '-'), ')'),
               '(1984-2014)')

  if(weighted) meanVar <- 'wtd.mean.slope' else meanVar <- 'mean.slope'
  if(weighted) sdVar <- 'wtd.sd' else sdVar <- 'sd'

  gp <- ggplot(sumtab, aes(x=ageClass, y=!!as.name(meanVar), group=!!as.name(catvar), color=!!as.name(catvar))) +
    xlab('stand age') +
    ylab(paste0('n-weighted mean AGB trend ', tp)) +
    scale_x_discrete(breaks = sumtab$ageClass, labels = sumtab$ageClass) +
    geom_hline(yintercept=0, linetype = 2) +
    geom_line(linetype = 3) +
    geom_errorbar(aes(ymin=!!as.name(meanVar) - (!!as.name(sdVar) / sqrt(count)), ymax=!!as.name(meanVar) + (!!as.name(sdVar) / sqrt(count))),
                  width=.2)

  if(out==1) plot(gp) else return(gp)

}

plotZoneStatsIntervals <- function(files2plot = file.path('outputs', list.files('outputs/', pattern='zoneStats_summary_WBI_ecozone_')),
                          weighted = T, xVar = 'ageClass', groupVar = 'tp', catVar = NULL, ptype = 1, plotResult = TRUE) {

  Require::Require(c('ggplot2','stringr','dplyr'))

  ## compile summary tables
  sumtab <- do.call(rbind, lapply(files2plot, function(x) {
    y <- readRDS(x) %>%
      mutate(ageClass = factor(ageClass, levels=unique(ageClass), ordered=T))
    y$tp = factor(rep(str_sub(x, start=-6L, end=-5L), nrow(y)),
                  levels = paste0('t', 1:6), ordered=T)
    return(y %>% relocate(tp, .before=value))
  }))

  ## validate arguments
  if(weighted) meanVar <- 'wtd.mean.slope' else meanVar <- 'mean.slope'
  if(weighted) sdVar <- 'wtd.sd' else sdVar <- 'sd'

  if(is.null(catVar)) catVar <- names(sumtab)[5]
  if(!is.null(groupVar)) groupVar <- match.arg(groupVar, names(sumtab))
  xVar <- match.arg(xVar, names(sumtab))
  if(xVar == 'ageClass') xlabel <- 'stand age'
  if(xVar == 'tp') xlabel <- '5-year time period (1984-2014)'

  if(ptype == 1) {

    gp <- ggplot(sumtab, aes(x=!!as.name(xVar), y=!!as.name(meanVar), group=!!as.name(groupVar), color=!!as.name(groupVar))) +
      labs(x = xlabel,
           y = 'n-weighted mean AGB trend') +
      geom_hline(yintercept=0, linetype='dotted') +
      geom_line() +
      facet_grid(~ .data[[catVar]])

  } else {

    if(ptype == 2) {

      gp <- lapply(unique(sumtab[,catVar]), function(cat) {

        x <- filter(sumtab, !!as.name(catVar) == cat)

        return(ggplot(x, aes(x=!!as.name(xVar), y=!!as.name(meanVar), group=!!as.name(groupVar), color=!!as.name(groupVar))) +
                 geom_line() +
                 xlab(xlabel) +
                 ylab('n-weighted mean AGB trend') +
                 geom_hline(yintercept = 0, linetype='dotted') +
                 geom_line() +
                 ggtitle(cat))

      })
      names(gp) <- unique(sumtab[,catVar])

      # pdf("all.pdf")
      # invisible(lapply(gp, print))
      # dev.off()

    }

    if(!plotResult) return(gp)
    if(ptype == 1) print(gp) else lapply(gp, print)

  }

}








## Test code using extracted cell values. Do results agree? Tested and yes :)
# rstack <- as.data.frame(c(slopeRaster, weightRaster, zoneRaster))
#
# group_by(rstack, zone) %>%
#   filter(!is.na(slope) & !is.na(zone)) %>%
#   dplyr::summarize(mean = format(mean(slope), scientific=T),
#                    sd = format(sd(slope), scientific=T),
#                    wmean = format(sum(w*slope) / sum(w), scientific=T),
#                    meanXw = sum(w * slope),
#                    wmean2 = format(weighted.mean(slope, w, na.rm=T), scientific=T),
#                    wSD1 = format(weighted.sd(slope, w, na.rm=T), scientific=T),
#                    wSD2 = format(sdwt(slope, w), scientific=T))
