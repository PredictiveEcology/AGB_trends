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

Require::Require(c('sf','terra','stringr','dplyr'))

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

Require::Require('parallel')
no_cores <- 8
cl <- parallelly::makeClusterPSOCK(no_cores,
                                   default_packages = c("terra","gdalUtilities","stringr"),
                                   rscript_libs = .libPaths(),
                                   autoStop = TRUE)
clusterExport(cl, varlist=c('agbdsn', 'agbfiles', 'distdsn', 'distfiles', 'tilenames'))
parallel::clusterEvalQ(cl, {
  terraOptions(tempdir = '/mnt/scratch/trudolph/AGB_trends/terra',
               memmax = 25,
               memfrac = 0.6,
               progress = 1,
               verbose = T)
})

ptime <- system.time({

  ## 1 b) import raster tiles
  parLapply(cl, 1:length(tilenames), function(i) {

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

    ## derive pixel- and year-specific age since disturbance (ABoVE disturbed pixels exclusively; *only possible for years 1987-2012*),
    ageSinceDist <- rast(lapply(1984:2014, function(k) k - chid))

    ## create mask serving to identify relevant cells
    ageMask <- ifel(!is.na(ageSinceDist), 9999, NA)

    ## replace rage values with 9999 where ABoVE disturbances overlap
    rageMask <- mask(rage, ageMask, maskvalue=9999, updatevalue=9999)
    rage <- cover(rageMask, ageSinceDist, values=9999)

    ## verify that ABoVE Disturbance product is correctly integrated
    if(FALSE) {
      k = sample(1:31, size=1)
      x <- mask(rage[[k]], ageMask[[k]])
      if(global(x, fun='sum', na.rm=T)$sum != global(ageSinceDist[[k]], fun='sum', na.rm=T)$sum) stop('PROBLEM ...!')
    }

    ## render NA all cell values pre-dating disturbances (where stand age cannot be known)
    rage <- ifel(rage > 0, rage, NA,
                 filename = file.path(tdir, paste0('rage_', tilenames[i], '.tif')), overwrite=T)

    ## remove AGB values where (pre-disturbance) year unknown
    mask(ragb, rage,
         filename = paste0(tdir, '/ragb_', tilenames[i], '.tif'),
         overwrite=T)

    rm(list=c('ragb','rage','rageMask','ageMask','ageSinceDist','rdist','chid'))

    gc()

  })

  stopCluster(cl)
  file.remove(file.path('/mnt/scratch/trudolph/AGB_trends/terra', list.files('/mnt/scratch/trudolph/AGB_trends/terra')))

}) # 1.6 hours


################################################################################
## Create mask limiting input values to cells disturbed during time series
## (according to ABoVE)

distdsn <- 'inputs/raw/ABoVE_ForestDisturbance_Agents/data/'
distfiles <- list.files(distdsn, pattern='.tif')

no_cores <- 25
cl <- parallelly::makeClusterPSOCK(no_cores,
                                   default_packages = c('terra','stringr'),
                                   rscript_libs = .libPaths(),
                                   autoStop = TRUE)

parallel::clusterExport(cl, varlist = 'distdsn', envir = environment())
parallel::clusterEvalQ(cl, {
  terraOptions(tempdir = '/mnt/scratch/trudolph/AGB_trends/terra',
               memmax = 25,
               memfrac = 0.8,
               progress = 1,
               verbose = T)
})

parallel::parLapply(cl, distfiles, function(distfile) {
  writeRaster(any(rast(file.path(distdsn, distfile))),
              filename = paste0('inputs/raw/ABoVE_ForestDisturbance_Agents/binary_disturbed/',
                                str_sub(distfile, end=7L), '_disturbed.tif'), overwrite = T)
})

parallel::stopCluster(cl)

## build raster mosaic
gdalUtilities::gdalbuildvrt(file.path('inputs/raw/ABoVE_ForestDisturbance_Agents/binary_disturbed',
                                      dir('inputs/raw/ABoVE_ForestDisturbance_Agents/binary_disturbed')),
                            '/mnt/scratch/trudolph/AGB_trends/terra/binary_dist_mosaic.vrt')
gdalUtilities::gdalwarp('/mnt/scratch/trudolph/AGB_trends/terra/binary_dist_mosaic.vrt',
                        'inputs/raw/ABoVE_ForestDisturbance_Agents/binary_disturbed_mosaic.tif')




################################################################################
## Import and pre-process ABoVE Landcover time series product (1984 - 2014)

## Set up access to EarthData webserver (in bash)
#Set up your ~/.netrc file as listed here: https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+cURL+And+Wget

# cd /home/trudolph
# touch .netrc
# echo "machine urs.earthdata.nasa.gov login trudolph password 4RrdVh.KqN?Rz7G" > .netrc
# chmod 0600 .netrc
# touch .urs_cookies

# Require::Require("httr")
# library(httr)
# netrc_path <- "/home/trudolph/.netrc"
# cookie_path <- "/home/trudolph/.urs_cookies"
# downloaded_file_path <- "/home/trudolph/test"
# set_config(config(followlocation=1,netrc=1,netrc_file=netrc_path,cookie=cookie_path,cookiefile=cookie_path,cookiejar=cookie_path))

## Get metadata
httr::GET(url = "https://daac.ornl.gov/daacdata/above/Annual_Landcover_ABoVE/comp/Annual_Landcover_ABoVE.pdf",
          write_disk("/home/trudolph/GitHub/AGB_trends/inputs/raw/ABoVE_Landcover/Annual_Landcover_ABoVE.pdf", overwrite = TRUE))

httr::GET(url = 'https://daac.ornl.gov/daacdata/above/Annual_Landcover_ABoVE/data/accuracy_assess_1984-2014.csv',
          write_disk("/home/trudolph/GitHub/AGB_trends/inputs/raw/ABoVE_Landcover/accuracy_assess_1984-2014.csv", overwrite = TRUE))

httr::GET(url = 'https://daac.ornl.gov/daacdata/above/Annual_Landcover_ABoVE/data/accuracy_summary_1984-2014.csv',
          write_disk("/home/trudolph/GitHub/AGB_trends/inputs/raw/ABoVE_Landcover/accuracy_summary_1984-2014.csv"))

## Download individually tiled rasters
sfiles <- paste0('https://daac.ornl.gov/daacdata/above/Annual_Landcover_ABoVE/data/ABoVE_LandCover_', stringr::str_sub(dir("inputs/raw/ABoVE_AGB_30m/data"), start=-11L, end=-5L), '.tif')

lapply(sfiles, function(sfile) {
  httr::GET(url = sfile,
            write_disk(paste0("/home/trudolph/GitHub/AGB_trends/inputs/raw/ABoVE_Landcover/ABoVE_LandCover_", str_sub(sfile, start=-11L, end=-5L), ".tif"), overwrite = TRUE))
  httr::GET(url = sfile,
            write_disk(paste0("/home/trudolph/GitHub/AGB_trends/inputs/raw/ABoVE_Landcover/ABoVE_LandCover_Simplified_", str_sub(sfile, start=-11L, end=-5L), ".tif"), overwrite = TRUE))
})


##############################################################################################
## Evaluate frequency of forested land cover types by ecozone * FOR IMPLEMENTATION AS MASK *

## Create land cover mosaic
gdalbuildvrt(gdalfile = file.path('inputs/raw/ABoVE_Landcover', paste0('ABoVE_LandCover_Simplified_', dir('inputs/clean/tiled'), '.tif')),
             output.vrt = 'cache/AGB_landCover_mosaic.vrt')

gdalwarp(srcfile = 'cache/AGB_landCover_mosaic.vrt', dstfile = 'outputs/AGB_landCover_mosaic.tif', overwrite=T)


## Rasterize ecozones to fit
ecoRast <- terra::rasterize(x = st_read('outputs/WBI_studyArea.gpkg', quiet=T),
                            y = rast('outputs/AGB_landCover_mosaic.tif'),
                            field = 'ECOZONE',
                            filename = 'cache/ecoRast.tif')

Require::Require('parallel')
no_cores <- 6
cl <- parallelly::makeClusterPSOCK(no_cores,
                                   default_packages = c("terra"),
                                   rscript_libs = .libPaths(),
                                   autoStop = TRUE)

parallel::clusterEvalQ(cl, {
  terraOptions(tempdir = '/mnt/scratch/trudolph/AGB_trends/terra',
               memmax = 25,
               memfrac = 0.5,
               progress = 1,
               verbose = T)
})

system.time({

  ftab <- parLapply(cl, cats(rast('cache/ecoRast.tif'))[[1]]$value, function(ezone) {
    return(freq(mask(rast('outputs/AGB_landCover_mosaic.tif'),
                     rast('cache/ecoRast.tif'),
                     maskvalues = ezone,
                     inverse = T)))
  })

})

stopCluster(cl)

saveRDS(ftab, 'outputs/ABoVE_LandCover_freq_tables.rds')

ftab <- readRDS('outputs/ABoVE_LandCover_freq_tables.rds')

  ## example looking at mean and sd count by ecozone
  group_by(ftab$`Boreal Shield`, value) %>%  summarize(meanVal = mean(count), sdVal = sd(count))

reftab <- data.frame(value = 1:15,
                     cat = c('Evergreen Forest', 'Deciduous Forest', 'Mixed Forest', 'Woodland',
                             'Low Shrub', 'Tall Shrub', 'Open Shrubs', 'Herbaceous', 'Tussock Tundra',
                             'Sparsely Vegetated', 'Fen', 'Bog', 'Shallows/littoral', 'Barren', 'Water'))

## Evaluate frequency distributions
ftab <- do.call(rbind, lapply(1:length(ftab), function(i) {
  return(bind_cols(ecozone = cats(rast('cache/ecoRast.tif'))[[1]]$ECOZONE[i],
                   group_by(ftab[[i]], value) %>%
                     summarize(count = mean(count, na.rm=T)) %>%
                     left_join(., reftab, by='value') %>%
                     relocate(cat, .before = count)))
}))

mutate(ftab, cat2 = ifelse(value %in% c(1:4), 'Forested', 'Non-Forested')) %>%
  group_by(ecozone, cat2) %>%
  summarize(catCount = sum(count)) %>%
  group_by(ecozone) %>%
  summarize(pcntForestedPixels = catCount[1] / sum(catCount))

##############################################################################################################


###################
## END OF SCRIPT ##
###################
