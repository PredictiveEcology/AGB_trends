defineModule(sim, list(
  name = "AGB_dataPrep",
  description = paste("Import 1) above-ground biomass (AGB) and 2) disturbance (year/type) time",
                      "series datasets created through the Arctic-Boreal Vulnerability Experiment",
                      "(ABoVE) research project (https://daac.ornl.gov/cgi-bin/dataset_lister.pl?p=34);",
                      "import 3) CaNFIR kNN stand age estimation (2020) data and create",
                      "spatial reference polygons corresponding to individual ABoVE tiles."),
  keywords = "", ## TODO
  authors = c(
    person(c("Tyler", "D"), "Rudolph", email = "tyler.rudolph@nrcan-rncan.gc.ca", role = c("aut", "cre")),
    person("Alex M", "Chubaty", email = "achubaty@for-cast.ca", role = c("ctb"))
  ),
  childModules = character(0),
  version = list(AGB_dataPrep = "0.0.0.9000"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.md", "AGB_dataPrep.Rmd"), ## same file
  reqdPkgs = list("dplyr", "geodata", "ggplot2", "ggspatial", "googledrive", "parallel", "purrr",
                  "sf", "stringr", "terra",
                  "PredictiveEcology/reproducible@development",
                  "PredictiveEcology/SpaDES.core@development (>= 1.1.0.9017)"),
  parameters = bindrows(
    defineParameter("analysisZonesType", "character", "ecozone", NA, NA,
                    paste("spatial scale at which to apply trend analyses;",
                          "one of 'ecozone', 'ecoprovince', 'ecoregion', or 'ecodistrict'.")),
    defineParameter("urlAGBTiles", "character", "https://drive.google.com/drive/folders/1CdJ4t_Ja0qcQgk5DIkFi8QSYQVcP9EV4", NA, NA,
                    "Google Drive URL to ABoVE Above Ground Biomass tiles."),
    defineParameter("urlForDistTiles", "character", "https://drive.google.com/drive/folders/1CNalAGmw9fO0-TuMMqLt3HrMNTj6cnER", NA, NA,
                    "Google Drive URL to ABoVE Forest Disturbance tiles."),
    defineParameter("urlCaNFIR", "character", "https://drive.google.com/file/d/1Thvxc9I8d1DzE7_hMovKESaCQGuHEYKy", NA, NA,
                    "Google Drive URL to CaNFIR stand age raster."),
    defineParameter(".plots", "character", c("screen", "png", "raw"), NA, NA,
                    "Used by Plots function, which can be optionally used here"),
    defineParameter(".plotInitialTime", "numeric", start(sim), NA, NA,
                    "Describes the simulation time at which the first plot event should occur."),
    defineParameter(".plotInterval", "numeric", NA, NA, NA,
                    "Describes the simulation time interval between plot events."),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA,
                    "Describes the simulation time at which the first save event should occur."),
    defineParameter(".saveInterval", "numeric", NA, NA, NA,
                    "This describes the simulation time interval between save events."),
    defineParameter(".studyAreaName", "character", NA, NA, NA,
                    "Human-readable name for the study area used - e.g., a hash of the study",
                          "area obtained using `reproducible::studyAreaName()`"),
    ## .seed is optional: `list('init' = 123)` will `set.seed(123)` for the `init` event only.
    defineParameter(".seed", "list", list(), NA, NA,
                    "Named list of seeds to use for each event (names)."),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    "Should caching of events or module be used?")
  ),
  inputObjects = bindrows(
    expectsInput("studyArea", "sf",
                 desc = paste("Polygon to use as the study area. Must be provided by the user"),
                 sourceURL = NA),
  ),
  outputObjects = bindrows(
    createsOutput("analysisZones", "sf",
                  desc = "spatial units for which trend analyses will be applied and reported.")
  )
))

## event types
#   - type `init` is required for initialization

doEvent.AGB_dataPrep = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      # do stuff for this event
      sim <- Init(sim)

      # schedule future event(s)
      sim <- scheduleEvent(sim, start(sim), "AGB_dataPrep", "download", .first())
      sim <- scheduleEvent(sim, start(sim), "AGB_dataPrep", "createABoVEPolys")
      sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "AGB_dataPrep", "plot")
      #sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "AGB_dataPrep", "save")
    },
    download = {
      sim <- downloadFromGoogleDrive(sim)
    },
    createABoVEPolys = {
      sim <- createAGBPolys(sim)
      sim <- createForDistPolys(sim)
    },
    plot = {
      # ! ----- EDIT BELOW ----- ! #

      plotFun(sim)

      # ! ----- STOP EDITING ----- ! #
    },
    save = {
      # ! ----- EDIT BELOW ----- ! #
      # do stuff for this event

      # e.g., call your custom functions/methods here
      # you can define your own methods below this `doEvent` function

      # schedule future event(s)

      # e.g.,
      # sim <- scheduleEvent(sim, time(sim) + P(sim)$.saveInterval, "AGB_dataPrep", "save")

      # ! ----- STOP EDITING ----- ! #
    },
    warning(paste("Undefined event type: \'", current(sim)[1, "eventType", with = FALSE],
                  "\' in module \'", current(sim)[1, "moduleName", with = FALSE], "\'", sep = ""))
  )
  return(invisible(sim))
}

## event functions
#   - keep event functions short and clean, modularize by calling subroutines from section below.

### template initialization
Init <- function(sim) {
  mod$dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)

  ## "Canada_Albers_Equal_Area_Conic" - no recognized EPSG code, using wkt:
  mod$targetCRS <- paste0("PROJCRS[\"Canada_Albers_Equal_Area_Conic\",\n",
                          "    BASEGEOGCRS[\"NAD83\",\n",
                          "        DATUM[\"North American Datum 1983\",\n",
                          "            ELLIPSOID[\"GRS 1980\",6378137,298.257222101004,\n",
                          "                LENGTHUNIT[\"metre\",1]]],\n",
                          "        PRIMEM[\"Greenwich\",0,\n",
                          "            ANGLEUNIT[\"degree\",0.0174532925199433]],\n",
                          "        ID[\"EPSG\",4269]],\n    CONVERSION[\"unnamed\",\n",
                          "        METHOD[\"Albers Equal Area\",\n",
                          "            ID[\"EPSG\",9822]],\n",
                          "        PARAMETER[\"Latitude of false origin\",40,\n",
                          "            ANGLEUNIT[\"degree\",0.0174532925199433],\n",
                          "            ID[\"EPSG\",8821]],\n",
                          "        PARAMETER[\"Longitude of false origin\",-96,\n",
                          "            ANGLEUNIT[\"degree\",0.0174532925199433],\n",
                          "            ID[\"EPSG\",8822]],\n",
                          "        PARAMETER[\"Latitude of 1st standard parallel\",50,\n",
                          "            ANGLEUNIT[\"degree\",0.0174532925199433],\n",
                          "            ID[\"EPSG\",8823]],\n",
                          "        PARAMETER[\"Latitude of 2nd standard parallel\",70,\n",
                          "            ANGLEUNIT[\"degree\",0.0174532925199433],\n",
                          "            ID[\"EPSG\",8824]],\n",
                          "        PARAMETER[\"Easting at false origin\",0,\n",
                          "            LENGTHUNIT[\"metre\",1],\n",
                          "            ID[\"EPSG\",8826]],\n",
                          "        PARAMETER[\"Northing at false origin\",0,\n",
                          "            LENGTHUNIT[\"metre\",1],\n",
                          "            ID[\"EPSG\",8827]]],\n",
                          "    CS[Cartesian,2],\n",
                          "        AXIS[\"easting\",east,\n",
                          "            ORDER[1],\n",
                          "            LENGTHUNIT[\"metre\",1,\n",
                          "                ID[\"EPSG\",9001]]],\n",
                          "        AXIS[\"northing\",north,\n",
                          "            ORDER[2],\n",
                          "            LENGTHUNIT[\"metre\",1,\n",
                          "                ID[\"EPSG\",9001]]]]")

  # # ! ----- EDIT BELOW ----- ! #

  sim$studyArea <- st_transform(sim$studyArea, mod$targetCRS)

  mod$AGBTilesPath <- file.path(mod$dPath, "ABoVE_AGB_30m", "data")
  mod$ForDistTilesPath <- file.path(mod$dPath, "ABoVE_ForestDisturbance_Agents", "data")

  mod$AGBTiles <- drive_ls(P(sim)$urlAGBTiles)
  mod$ForDistTiles <- drive_ls(P(sim)$urlForDistTiles)

  urlList <- list(
    ecodistrict = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/district/ecodistrict_shp.zip",
    ecoregion = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/region/ecoregion_shp.zip",
    ecoprovince = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/province/ecoprovince_shp.zip",
    ecozone = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/zone/ecozone_shp.zip"
  )
  url <- urlList[[tolower(P(sim)$analysisZonesType)]]

  sim$analysisZones <- prepInputs(url = url,
                                  destinationPath = mod$dPath,
                                  studyArea = sim$studyArea,
                                  fun = "sf::st_read",
                                  overwrite = TRUE) %>%
    st_transform(mod$targetCRS)

  # ! ----- STOP EDITING ----- ! #

  return(invisible(sim))
}

### template for save events
Save <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #

  sim <- saveFiles(sim)

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

downloadFromGoogleDrive <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #

  ##  NOTE: we _could_ do parallel downloads, but bandwidth will always bottleneck;
  ##        it's likely faster to download each file at max speed, rather than multiple slowly

  ## ABoVE above ground biomass tiles --------------------------------------------------------------
  fnames1 <- file.path(mod$AGBTilesPath, mod$AGBTiles$name)

  needDownload1 <- which(!file.exists(fnames1)) ## TODO: also checksum the files

  if (length(needDownload1) > 0) {
    retry(quote({
      httr::with_config(config = httr::config(http_version = 2), {
        map2_dfr(mod$AGBTiles$id[needDownload1], fnames1[needDownload1], drive_download)
      })
    }))
  }

  ## ABoVE forest disturbance tiles ----------------------------------------------------------------
  fnames2 <- file.path(mod$ForDistPath, mod$ForDistTiles$name)
  needDownload2 <- which(!file.exists(fnames2)) ## TODO: also checksum the files

  if (length(needDownload2) > 0) {
    retry(quote({
      httr::with_config(config = httr::config(http_version = 2), {
        map2_dfr(mod$ForDistTiles$id[needDownload2], fnames2[needDownload2], drive_download)
      })
    }))
  }

  ## CaNFIR:
  ## NOTE: postProcessing (i.e. crop & re-project) will be wrapped in analysis (much better performance)
  out <- Cache(preProcess,
               url = P(sim)$urlCaNFIR,
               targetFile = "mosaic_age.tif",
               fun = "raster::raster",
               alsoExtract = "similar",
               destinationPath = file.path(mod$dPath, "CaNFIR"),
               #overwrite = TRUE,
               userTags = c("CaNFIR")) ## keywords to use in the cache db to make it easy to find later

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

## Create vector polygons corresponding to ABoVE AGB raster tiles
createAGBPolys <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #

  gpkgFile <- file.path(dirname(mod$AGBTilesPath), "ABoVE_AGB_study_area.gpkg")

  rfiles <- list.files(mod$AGBTilesPath, pattern = "AGB_B")

  ## individual tile bounding boxes
  boxes <- do.call(rbind, lapply(rfiles, function(x) as.vector(ext(rast(file.path(mod$AGBTilesPath, x))))))

  ## wide area bounding box
  Abox <- st_as_sf(st_as_sfc(st_bbox(c(apply(boxes, 2, min)[c(1, 3)], apply(boxes, 2, max)[c(2, 4)])[c(1, 3, 2, 4)])))
  st_crs(Abox) <- crs(rast(file.path(mod$AGBTilesPath, rfiles[1])))

  if (!file.exists(gpkgFile)) {
    st_write(merge(Abox, data.frame(description = "ABoVE_AGB_study_area")),
             dsn = gpkgFile, layer = "study_area",  driver = "GPKG", delete_layer = TRUE)
  }

  ## append individual tiles
  vtiles <- do.call(rbind, lapply(1:nrow(boxes), function(m) {
    st_as_sf(st_as_sfc(st_bbox(c(boxes[m, 1], boxes[m, 3], boxes[m, 2], boxes[m, 4]))),
             data.frame(tile_name = rfiles[m]))
  }))
  st_crs(vtiles) <- crs(rast(file.path(mod$AGBTilesPath, rfiles[1])))

  if (!file.exists(gpkgFile)) {
    st_write(vtiles, dsn = gpkgFile, layer = "tileset", driver = "GPKG", delete_layer = TRUE)
  }

  mod$agbTiles <- vtiles

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

## Create vector polygons corresponding to ABoVE disturbance history raster tiles
createForDistPolys <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #

  gpkgFile <- file.path(dirname(mod$ForDistTilesPath), "ABoVE_DistAgents_study_area.gpkg")


  rfiles <- list.files(mod$ForDistTilesPath, pattern = ".tif")

  ## individual tile bounding boxes
  boxes <- do.call(rbind, lapply(rfiles, function(x) {
    as.vector(ext(rast(file.path(mod$ForDistTilesPath, x))))
  }))

  ## wide area bounding box
  Abox <- st_as_sf(st_as_sfc(st_bbox(c(apply(boxes, 2, min)[c(1, 3)], apply(boxes, 2, max)[c(2, 4)])[c(1, 3, 2, 4)])))
  st_crs(Abox) <- st_read(dsn = gpkgFile, layer = "study_area") %>%
    st_crs() # identical CRS to AGB product, but mis-specified

  if (!file.exists(gpkgFile)) {
    st_write(merge(Abox, data.frame(description = "ABoVE_ForestDisturbance_Agents_study_area")),
             dsn = gpkgFile, layer = "study_area", driver = "GPKG", delete_layer = TRUE)
  }

  ## append individual tiles
  vtiles <- do.call(rbind, lapply(1:nrow(boxes), function(m) {
    st_as_sf(st_as_sfc(st_bbox(c(boxes[m, 1], boxes[m, 3], boxes[m, 2], boxes[m, 4]))),
             data.frame(tile_name = rfiles[m]))
  }))
  st_crs(vtiles) <- st_read(dsn = gpkgFile, layer = "study_area") %>% st_crs() # identical CRS to AGB product, but mis-specified

  if (!file.exists(gpkgFile)) {
    st_write(vtiles, dsn = gpkgFile, layer = "tileset", driver = "GPKG", delete_layer = TRUE)
  }

  mod$dstTiles <- vtiles

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

plotFun <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #

  agbTilesToPlot <- filter(mod$agbTiles, apply(relate(x = vect(mod$agbTiles), y = vect(sim$studyArea),
                                                      relation = "intersects"), 1, any))
  dstTilesToPlot <- filter(mod$dstTiles, apply(relate(x = vect(mod$dstTiles), y = vect(sim$studyArea),
                                                      relation = "intersects"), 1, any))

  Plots(sim$analysisZones, agbTilesToPlot, dstTilesToPlot, fn = plotstudyAreaCoverage)

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

plotstudyAreaCoverage <- function(studyArea, agbTiles, dstTiles, ...) {
  alpha <- 0.3
  ggplot(studyArea) +
    geom_sf(aes(fill = ZONE_NAME), colour = "black", alpha = alpha) +
    theme_bw() +
    annotation_north_arrow(location = "bl", which_north = "true",
                           pad_x = unit(0.25, "in"), pad_y = unit(0.25, "in"),
                           style = north_arrow_fancy_orienteering) +
    xlab("Longitude") + ylab("Latitude") +
    ggtitle("Analysis zones within study area") +
    geom_sf(data = st_as_sf(dstTiles), colour = "red", alpha = alpha) +
    geom_sf(data = st_as_sf(agbTiles), colour = "lightblue", alpha = alpha)
}

.inputObjects <- function(sim) {
  cacheTags <- c(currentModule(sim), "function:.inputObjects")
  mod$dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": using dataPath '", mod$dPath, "'.")

  # ! ----- EDIT BELOW ----- ! #

  if (!suppliedElsewhere("studyArea", sim)) {
    stop("studyArea must be supplied by user")
  }

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}
