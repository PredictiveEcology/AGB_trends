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

Require::Require(c('terra','stringr','gdalUtilities','dplyr','sf'))
terraOptions(tempdir='/mnt/scratch/trudolph/AGB_trends/cache', todisc=TRUE, progress=1)

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

################################################################################
## 1 a) Estimate cell-wise linear regression coefficients for undisrupted time series
## aka "local" or "geographically weighted regression (GWR)"

system.time({
  
  sapply(1:length(tile_folders), function(i) {
    
    ## calculate local slope coefficient
    # app(rast(file.path(tile_folders[i], list.files(tile_folders[i], pattern='ragb'))), 
    #     fun=slope, cores=32,
    #     filename=file.path(tile_folders[i], paste0('agb_slopes_', str_sub(tile_folders[i], start=-7L),'.tif')),
    #     overwrite=T)
    
    ## stock number of non-NA values for subsequent weighted standard deviation
    nsamp <- function(x) sum(!is.na(x))
    app(rast(file.path(tile_folders[i], list.files(tile_folders[i], pattern='ragb'))), 
        fun=nsamp, cores=32, 
        filename=file.path(tile_folders[i], paste0('agb_sample_size_', str_sub(tile_folders[i], start=-7L),'.tif')),
        overwrite=T)
    
  })
  
}) # 3.8 hrs + 3 hrs for nsamp

## 1 b) Combine tiled slope rasters into unified mosaic

## Build virtual rasters
gdalbuildvrt(gdalfile = unname(sapply(tile_folders, function(dsn) file.path(dsn, list.files(dsn, pattern='agb_slopes_Bh')))),
             output.vrt = 'cache/AGB_slope_mosaic.vrt')

## Write to raster mosaics
gdalwarp(srcfile = 'cache/AGB_slope_mosaic.vrt', dstfile = 'inputs/ABoVE_AGB_30m/AGB_slope_mosaic.tif')


#################################################################################################################
## 2 a) Calculate cell-specific slopes per 5-year time interval (n=6)

## define time intervals (year ranges between 1984-2014)
timeint <- list(t1=1:5, t2=6:10, t3=11:15, t4=16:20, t5=21:25, t6=26:31)

system.time({
  
  for(i in 1:length(tile_folders)) {
    
    cat(i, '/', length(tile_folders), '\n')
    
    tile_folder <- tile_folders[i]
    
    sapply(1:length(timeint), function(timestep) {
      
      # ## calculate local slope coefficient for specified time interval
      # app(rast(file.path(tile_folders[i], list.files(tile_folders[i], pattern='ragb')))[[timeint[[timestep]]]],
      #     fun=slope, cores=45,
      #     filename=file.path(tile_folders[i], paste0('agb_slopes_', names(timeint)[timestep], '_', str_sub(tile_folders[i], start=-7L),'.tif')),
      #     overwrite=T)
      
      ## stock number of non-NA values for subsequent weighted standard deviation
      nsamp <- function(x) sum(!is.na(x))
      app(rast(file.path(tile_folders[i], list.files(tile_folders[i], pattern='ragb'))), 
          fun=nsamp, cores=32, 
          filename=file.path(tile_folders[i], paste0('agb_sample_size_', names(timeint)[timestep], '_', str_sub(tile_folders[i], start=-7L),'.tif')),
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

#########################################################################################
## Visual examination of results !! write to png !!

par(mfrow=c(3,2))
for(tp in names(timeint)) {
  plot(rast(paste0('inputs/ABoVE_AGB_30m/AGB_slope_mosaic_', tp, '.tif')), main=tp)
}

sapply(names(timeint), function(tp) digest::digest(paste0('inputs/ABoVE_AGB_30m/AGB_slope_mosaic_', tp, '.tif'), algo='xxhash64'))
sapply(names(timeint), function(tp) file.size(paste0('inputs/ABoVE_AGB_30m/AGB_slope_mosaic_', tp, '.tif')))

# system.time(test <- diff(rast(sapply(names(timeint)[c(1,6)], function(tp) rast(paste0('inputs/ABoVE_AGB_30m/AGB_slope_mosaic_', tp, '.tif')))),
#             filename='cache/diff.tif', overwrite=T))

## write to PNG !!
par(mfrow=c(2,3))
hist(rast(paste0('inputs/ABoVE_AGB_30m/AGB_slope_mosaic_t1.tif')), main='t1')
hist(rast(paste0('inputs/ABoVE_AGB_30m/AGB_slope_mosaic_t2.tif')), main='t2')
hist(rast(paste0('inputs/ABoVE_AGB_30m/AGB_slope_mosaic_t3.tif')), main='t3')
hist(rast(paste0('inputs/ABoVE_AGB_30m/AGB_slope_mosaic_t4.tif')), main='t4')
hist(rast(paste0('inputs/ABoVE_AGB_30m/AGB_slope_mosaic_t5.tif')), main='t5')
hist(rast(paste0('inputs/ABoVE_AGB_30m/AGB_slope_mosaic_t6.tif')), main='t6')

############################################################################################################
## 3) Group slopes by age at time 0 (band=1), masking out pixels disturbed mid-time series

## 3.1) Build virtual raster of estimated stand age at time 0 (i.e. 1984)
Require::Require('gdalUtilities')
gdalbuildvrt(
  gdalfile = sapply(tile_folders, function(dsn) unname(file.path(dsn, list.files(dsn, pattern='rage')))),
  output.vrt = 'cache/AGB_age_mosaic.vrt',
  b=1 # time 0 (1984) = band or layer 1
)

## 3.2) Write to raster mosaic (stand age, kNN 2020)
gdalwarp(srcfile = 'cache/AGB_age_mosaic.vrt', 
         dstfile = 'inputs/ABoVE_AGB_30m/AGB_age_mosaic.tif')

## 3.3) Group into 20-year age classes
classify(rast('inputs/ABoVE_AGB_30m/AGB_age_mosaic.tif'),
         rcl=cbind(from=c(0, 25, 50, 80, 125), to=c(25, 50, 80, 125, 500), becomes=1:5), 
         right=FALSE, others=NA, overwrite=TRUE,
         filename='inputs/ABoVE_AGB_30m/AGB_age_mosaic_classes.tif')

####################################################
## 4) Calculate mean of slopes by ecozone

age.class <- data.frame(from=c(0, 25, 50, 80, 125), to=c(25, 50, 80, 125, 500), becomes=1:5)
wbi.class <- st_read('inputs/WBI/WBI_studyArea.gpkg', quiet=T) %>% st_drop_geometry()

## 4.1) rasterize WBI ecozones
rasterize(x = st_read('inputs/WBI/WBI_studyArea.gpkg', quiet=T) %>% vect(.), 
          y = rast('inputs/ABoVE_AGB_30m/AGB_age_mosaic_classes.tif'), 
          field = 'BCR',
          filename = 'inputs/WBI/WBI_studyArea.tif', 
          overwrite = T)

## 4.2) combine ecozone & age class into unique combined raster categories 
writeRaster(rast('inputs/ABoVE_AGB_30m/AGB_age_mosaic_classes.tif') + 
              rast('inputs/WBI/WBI_studyArea.tif') * 10,
            filename='inputs/ABoVE_AGB_30m/AGB_age_ecozone_combined.tif',
            overwrite=T) # ? min

## 4.3) compute mean slope per pixel by unique ecozone + age category
zonal(rast('inputs/ABoVE_AGB_30m/AGB_slope_mosaic.tif'),
      rast('inputs/ABoVE_AGB_30m/AGB_age_ecozone_combined.tif'),
      fun='mean', na.rm=T, as.raster=T, 
      filename='inputs/ABoVE_AGB_30m/mean_slope_x_ecozone_x_age.tif') # 12 min

## 4.4) compute tabular mean, sample size, and sd (via formula numerator) per unique group
system.time({
  
  zonal(rast('inputs/ABoVE_AGB_30m/AGB_slope_mosaic.tif'),
        rast('inputs/ABoVE_AGB_30m/AGB_age_ecozone_combined.tif'),
        fun='mean', na.rm=T) %>%
    rename(code = AGB_age_mosaic,
           mean.trend = AGB_slope_mosaic) %>%
    
    ## sample size
    left_join(freq(rast('inputs/ABoVE_AGB_30m/AGB_age_ecozone_combined.tif'))[,-1], by=c('code'='value')) %>%
    
    ## sd formula numerator
    left_join(zonal(x = (rast('inputs/ABoVE_AGB_30m/AGB_slope_mosaic.tif') - 
                           rast('inputs/ABoVE_AGB_30m/mean_slope_x_ecozone_x_age.tif')) ^ 2, 
                    z = rast('inputs/ABoVE_AGB_30m/AGB_age_ecozone_combined.tif'),
                    fun = 'sum', na.rm=T) %>%
                rename(sd.form.num = AGB_slope_mosaic), by=c('code'='AGB_age_mosaic')) %>%
    
    ## calculate SD
    mutate(sd = sqrt(sd.form.num / count)) %>%
    
    mutate(ecozone.id = str_sub(code, end=1L),
           age.id = str_sub(code, start=2L),
           age.min = age.class$from[match(age.id, age.class$becomes)],
           age.max = age.class$to[match(age.id, age.class$becomes)],
           ecozone = wbi.class$Label[match(ecozone.id, wbi.class$BCR)], .before=mean.trend) %>%
    
    ## save to file
    saveRDS(., file='output/slopeXageXecozone.rds')
  
})

##############################################################
## Visualize
Require::Require('ggplot2')
readRDS('output/slopeXageXecozone.rds') %>%
  ggplot(aes(x=age.id, y=AGB, group=ecozone, color=ecozone)) +
  xlab('stand age') + 
  ylab('mean AGB trend (1984-2014)') +
  scale_x_discrete(breaks=slopesum$age.id,
                   labels=paste0(slopesum$age.min, '-', slopesum$age.max)) +
  geom_hline(yintercept=0, linetype='dotted') +
  geom_line()





app(rast('inputs/ABoVE_AGB_30m/AGB_slope_mosaic.tif'), fun=function(x) {
  
  x - 
  
}, cores=, filename=, )


for(i in slopesum$code)

  i <- slopesum$code[1]

  imean <- slopesum$AGB[slopesum$code==i]
  # test <- mask(rast('inputs/ABoVE_AGB_30m/AGB_slope_mosaic.tif'), rast('inputs/ABoVE_AGB_30m/AGB_age_ecozone_combined.tif')==i)







