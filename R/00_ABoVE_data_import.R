##############################################################################################################################################
##############################################################################################################################################
##
## Date Created: Nov 17, 2022
## Auteur: Tyler Rudolph, biologist M.Sc., CFS/NRCAN, Trade, Economics & Industry Branch
##
## Name of script : "00_ABoVE_data_import.R"
## Description : Import 1) above-ground biomass (AGB) and 2) disturbance (year/type) time series datasets
##               created through the Arctic-Boreal Vulnerability Experiment (ABoVE) research project
##               URL = https://daac.ornl.gov/cgi-bin/dataset_lister.pl?p=34;
##               3) Import CaNFIR kNN stand age estimation (2020)
##               4) Create spatial reference polygons corresponding to individual ABoVE tiles
##
##############################################################################################################################################
##############################################################################################################################################

## global parameters
Require::Require(c('terra','sf','stringr','dplyr','googledrive','reproducible'))
terraOptions(tempdir='scratch', todisc=TRUE)
googledrive::drive_auth(path="~/forprod-52bae8224d11.json")

#########################################
## 1) download tiled ABoVE AGB rasters
agbtiles <- drive_ls('https://drive.google.com/drive/folders/1CdJ4t_Ja0qcQgk5DIkFi8QSYQVcP9EV4')

## 1 a) setup parallel processing (one thread per tile)
library(parallel)
no_cores <- 32
cl <- makeCluster(no_cores)
clusterExport(cl, varlist=c('agbtiles'))

## 1 b) import raster tiles
parLapply(cl, 1:nrow(agbtiles), function(m) {
  library(googledrive)
  googledrive::drive_auth(path="forprod-52bae8224d11.json")
  retry(quote({
    httr::with_config(config = httr::config(http_version = 2), {
      googledrive::drive_download(file=as_id(agbtiles[m,]),
                                  path=file.path('inputs/ABoVE_AGB_30m/data', agbtiles$name[m]),
                                  overwrite=TRUE) # use check files instead + check sums
    })
  }))
})

stopCluster(cl)

#############################################################
## 2) download tiled ABoVE Forest Disturbance Agents
googledrive::drive_auth(path="forprod-52bae8224d11.json")
distiles <- drive_ls('https://drive.google.com/drive/folders/1CNalAGmw9fO0-TuMMqLt3HrMNTj6cnER')

## 2 a) setup parallel processing (one thread per tile)
library(parallel)
no_cores <- 32
cl <- makeCluster(no_cores)
clusterExport(cl, varlist=c('distiles'))

## 2 b) import raster tiles
parLapply(cl, 1:nrow(distiles), function(m) {
  library(googledrive)
  googledrive::drive_auth(path="forprod-52bae8224d11.json")
  retry(quote({
    httr::with_config(config = httr::config(http_version = 2), {
      drive_download(file=as_id(distiles[m,]),
                     path=file.path('inputs/ABoVE_ForestDisturbance_Agents/data', distiles$name[m]),
                     overwrite=TRUE)
    })
  }))
})

stopCluster(cl)

####################################################################################
## 3) download CaNFIR stand age product + supporting file types
terraOptions(tempdir='scratch', todisc=TRUE)

## create destination folder if non-existent
if(!dir.exists('ForProd/inputs/CaNFIR')) dir.create('ForProd/inputs/CaNFIR')

## authenticate w/Google Cloud
googledrive::drive_auth(path="forprod-52bae8224d11.json")
preProcess(url = 'https://drive.google.com/file/d/1Thvxc9I8d1DzE7_hMovKESaCQGuHEYKy',
           targetFile = "mosaic_age.tif",
           fun = "raster::raster",
           alsoExtract='similar',
           # .tempPath='scratch', # not necessary since postProcessing (i.e. crop & re-project) will be wrapped in analysis (much better performance)
           destinationPath = "inputs/CaNFIR",
           userTags = c("bigRasters", "CaNFIR"), ## keywords to use in the cache db to make it easy to find later
           overwrite=TRUE)

##########################################################################
## 4 a) Create vector polygons corresponding to ABoVE AGB raster tiles
dsn <- 'inputs/ABoVE_AGB_30m/data'
rfiles <- list.files(dsn, pattern='AGB_B')

## individual tile bounding boxes
boxes <- do.call(rbind, lapply(rfiles, function(x) as.vector(ext(rast(file.path(dsn, x))))))

## wide area bounding box
Abox <- st_as_sf(st_as_sfc(st_bbox(c(apply(boxes, 2, min)[c(1,3)], apply(boxes, 2, max)[c(2,4)])[c(1,3,2,4)])))
st_crs(Abox) <- crs(rast(file.path(dsn, rfiles[1])))

st_write(merge(Abox, data.frame(description='ABoVE_AGB_study_area')),
         dsn='inputs/ABoVE_AGB_30m/ABoVE_AGB_study_area.gpkg', layer='study_area',
         driver='GPKG', delete_layer=T)

## append individual tiles
vtiles <- do.call(rbind, lapply(1:nrow(boxes), function(m) st_as_sf(st_as_sfc(st_bbox(c(boxes[m,1], boxes[m,3], boxes[m,2], boxes[m,4]))), data.frame(tile_name=rfiles[m]))))
st_crs(vtiles) <- crs(rast(file.path(dsn, rfiles[1])))
st_write(vtiles, dsn='inputs/ABoVE_AGB_30m/ABoVE_AGB_study_area.gpkg', layer='tileset', driver='GPKG', delete_layer=T)


###########################################################################
## 4 b) Create vector polygons corresponding to ABoVE disturbance history raster tiles
dsn <- 'inputs/ABoVE_ForestDisturbance_Agents/data'
rfiles <- list.files(dsn, pattern='.tif')

## individual tile bounding boxes
boxes <- do.call(rbind, lapply(rfiles, function(x) as.vector(ext(rast(file.path(dsn, x))))))

## wide area bounding box
Abox <- st_as_sf(st_as_sfc(st_bbox(c(apply(boxes, 2, min)[c(1,3)], apply(boxes, 2, max)[c(2,4)])[c(1,3,2,4)])))
st_crs(Abox) <- st_read(dsn='inputs/ABoVE_AGB_30m/ABoVE_AGB_study_area.gpkg', layer='study_area') %>% st_crs() # identical CRS to AGB product, but mis-specified

st_write(merge(Abox, data.frame(description='ABoVE_ForestDisturbance_Agents_study_area')),
         dsn='inputs/ABoVE_ForestDisturbance_Agents/ABoVE_DistAgents_study_area.gpkg', layer='study_area',
         driver='GPKG', delete_layer=T)

## append individual tiles
vtiles <- do.call(rbind, lapply(1:nrow(boxes), function(m) st_as_sf(st_as_sfc(st_bbox(c(boxes[m,1], boxes[m,3], boxes[m,2], boxes[m,4]))), data.frame(tile_name=rfiles[m]))))
st_crs(vtiles) <- st_read(dsn='inputs/ABoVE_AGB_30m/ABoVE_AGB_study_area.gpkg', layer='study_area') %>% st_crs() # identical CRS to AGB product, but mis-specified
st_write(vtiles, dsn='inputs/ABoVE_ForestDisturbance_Agents/ABoVE_DistAgents_study_area.gpkg', layer='tileset', driver='GPKG', delete_layer=T)

###############################################################################
## 5) Import WBI study area (TO DO)
Require::Require(c('magrittr','sf','reproducible'))

cPath <- 'cache'
if(!dir.exists('inputs/WBI')) dir.create('inputs/WBI')

## ABoVE default CRS (Canada_Albers_Equal_Area_Conic)
targetCRS <- st_read(dsn='inputs/ABoVE_AGB_30m/ABoVE_AGB_study_area.gpkg', layer='study_area') %>% st_crs()

bcrzip <- "https://www.birdscanada.org/download/gislab/bcr_terrestrial_shape.zip"

bcrWB <- Cache(prepInputs,
                url = bcrzip,
                cacheRepo = cPath,
                destinationPath = 'inputs/WBI',
                targetCRS = targetCRS,
                fun = "sf::st_read") %>%
  filter(BCR %in% c(4, 6:8))

provsWB <- Cache(prepInputs,
                  url = 'https://www12.statcan.gc.ca/census-recensement/2011/geo/bound-limit/files-fichiers/2016/lpr_000b16a_e.zip',
                  destinationPath = 'cache',
                  targetFile='lpr_000b16a_e.shp',
                  alsoExtract='similar',
                  targetCRS = targetCRS,
                  fun = "sf::st_read") %>%
  filter(PREABBR %in% c('B.C.','Alta.','Sask.','Man.','Y.T.','N.W.T.','Nvt.')) %>%
  st_cast(., 'POLYGON')

## crop WBI to western Cdn provinces and save to file
studyArea <- bcrWB %>%
  st_crop(., provsWB[provsWB$PREABBR != 'Nvt.',]) %>%
  st_intersection(., provsWB) %>%
  group_by(BCR, Label) %>%
  summarize() %>%
  st_write(., dsn='inputs/WBI/WBI_studyArea.gpkg', driver='GPKG', delete_layer=T)

#################################################################################################################
## OPTIONAL: Visually compare available tiles between ABoVE products and WBI study area

## 1) AGB tiles
agb_sa <- st_read('inputs/ABoVE_AGB_30m/ABoVE_AGB_study_area.gpkg', 'study_area')
agb_tiles <- st_read('inputs/ABoVE_AGB_30m/ABoVE_AGB_study_area.gpkg', 'tileset')

## 2) Disturbance history tiles
dist_sa <- st_read('inputs/ABoVE_ForestDisturbance_Agents/ABoVE_DistAgents_study_area.gpkg', 'study_area')
dist_tiles <- st_read('inputs/ABoVE_ForestDisturbance_Agents/ABoVE_DistAgents_study_area.gpkg', 'tileset')

## 3) WBI study area
wbi <- st_read('inputs/WBI/WBI_studyArea.gpkg')

## 4) Comparison plot

# dir.create('figures')
# png ...

plot(wbi %>% st_geometry)
plot(filter(dist_tiles, apply(relate(x=vect(dist_tiles), y=vect(wbi), relation='intersects'), 1, any)) %>% st_geometry, border='red', add=T)
plot(filter(agb_tiles, apply(relate(x=vect(agb_tiles), y=vect(wbi), relation='intersects'), 1, any)) %>% st_geometry,
     border='lightblue', add=T)
plot(wbi %>% st_geometry, add=T)


###################
## END OF SCRIPT ##
###################
