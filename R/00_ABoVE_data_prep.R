##############################################################################################################################################
##############################################################################################################################################
##
## Date Created: Nov 17, 2022
## Auteur: Tyler Rudolph, biologist M.Sc., CFS/NRCAN, Trade, Economics & Industry Branch
##
## Name of script : "00_ABoVE_data_prep.R"
##
## Description : R script serving to estimate cell-specific stand age for raster tiles from the ABoVE AGB time series
##               using a combination of ABoVE Disturbance Agents and the CaNFIR stand age mosaic (kNN 2020)
##
##############################################################################################################################################
##############################################################################################################################################

Require::Require(c('sf','terra','stringr'))

######################################################################
## 1.1) Import ABoVE product tiles
dist_tiles <- st_read('outputs/ABoVE_DistAgents_study_area.gpkg', 'tileset')
agb_tiles <- st_read('outputs/ABoVE_AGB_study_area.gpkg', 'tileset')

## 1.2) Index 30m ABG input rasters
agbdsn <- 'inputs/raw/ABoVE_AGB_30m/data'
agbfiles <- list.files(agbdsn, pattern='AGB_B')

## 1.3) Index ABoVE disturbance history raster masks
distdsn <- 'inputs/raw/ABoVE_ForestDisturbance_Agents/data/'
distfiles <- list.files(distdsn, pattern='.tif')

## 1.4) Import WBI study area
wbi <- st_read('outputs/WBI_studyArea.gpkg')

## 1.5) Reduce input file selection to tiles contained within WBI study area & align filename indices
agbfiles <- agbfiles[apply(relate(x=vect(agb_tiles), y=vect(wbi), relation='intersects'), 1, any)]
tilenames <- str_sub(str_remove_all(agbfiles, 'ABoVE_AGB_'), end=-5L)
distfiles <- distfiles[match(tilenames, str_sub(distfiles, end=7L))]


###########################################################################
## 2) Estimate cell-specific stand age using a combination of 
##    ABoVE Disturbance Agents and the CaNFIR stand age mosaic (kNN 2020)

library(parallel)
no_cores <- 10
cl <- makeCluster(no_cores)
clusterExport(cl, varlist=c('agbdsn', 'agbfiles', 'distdsn', 'distfiles', 'tilenames'))

ptime <- system.time({
  
  ## 1 b) import raster tiles
  parLapply(cl, 1:length(tilenames), function(i) {
    
    library(terra)
    terraOptions(tempdir='scratch', todisc=TRUE)
    
    tdir <- file.path('inputs/clean/tiled', tilenames[i])
    if(!file.exists(tdir)) dir.create(tdir)
    
    ## Import AGB raster tile (0 = NA)
    ragb <- rast(file.path(agbdsn, agbfiles[i]))
    names(ragb) <- 1984:2014
    
    ## Restrict selection to available disturbance history years
    NAflag(ragb) <- 0
    
    ## 1 c) step i) project stand age mosaic (kNN 2020) to ABoVE CRS & crop to tile
    ## 1 c) step ii) determine stand age as of year associated with raster layer (negative values become NA)
    ## 1 c) step iii) split into 5 year age categories
    rage <- classify(as.integer(names(ragb)) -
                       (2020 - crop(project(rast('inputs/raw/CaNFIR/mosaic_age.tif'),
                                            ragb[[1]], method='near'), ragb[[1]])),
                     rcl=cbind(-100, 0, NA), include.lowest=T)
    
    names(rage) <- 1984:2014 # 3.63 min
    
    ###############################################################################
    ## Import disturbance history (ABoVE) & mask pixels in PRE-disturbance years 
    rdist <- rast(file.path(distdsn, distfiles[i]))
    names(rdist) <- 1987:2012
    
    ## identify year of earliest disturbance for each cell in tile (according to ABoVE product)
    chid <- app(
      rast(lapply(names(rdist), function(k) {
        classify(rdist[[k]], rcl = cbind(1, 3, as.integer(k)), include.lowest=TRUE)
      })), 
      fun=min, na.rm=T) # 106 sec
    
    ## Render NA all cell values pre-dating disturbances (where stand age cannot be known)
    ## & modify age for cells post-dating disturbances (according to ABoVE disturbance product)
    ## * only possible for years 1987-2012 *
    rage <- writeRaster(
      rast(lapply(1:nlyr(rage), function(j) {
      return(ifel(is.na(chid), rage[[j]], 
                  ifel(chid < as.integer(names(rage)[j]), as.integer(names(rage)[j]) - chid, NA)))
    })), filename = file.path(tdir, paste0('rage_', tilenames[i], '.tif')), overwrite=T) # 371 sec

    mask(ragb, rage,
         filename = paste0(tdir, '/ragb_', tilenames[i], '.tif'),
         overwrite=T) # ~ 2 min
    
    rm(list=c('ragb','rage','rdist','chid'))
    gc()
    
  })
  
  stopCluster(cl)
  
}) # 2.05 hours
  

###################
## END OF SCRIPT ##
###################
  
