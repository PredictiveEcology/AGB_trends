# packages ------------------------------------------------------------------------------------

library(data.table)
library(dplyr)
library(reproducible)
library(sf)
library(stringr)
library(terra)

library(AGBtrends)

# global parameters for project setup ---------------------------------------------------------

projName <- workflowtools::findProjectName()
studyAreaName <- "AGB_WBI"
user <- Sys.info()[["user"]]

## specify which WBI runs to use for AGB analyses
if (!exists("climateGCM")) {
  climateGCM <- "CanESM5" ## "CNRM-ESM2-1"
}
if (!exists("climateSSP")) {
  climateSSP <- 370 ## 585
}

climateScenario <- paste0(climateGCM, "_SSP", climateSSP)

allReps <- sprintf("%02d", 1:5)

lapply(allReps, function(thisRep) {
  message(paste("Processing", climateScenario, "rep", thisRep))
  allStudyAreas <- c("AB", "BC", "NT", "SK", "YT") ## MB doesn't overlap with ABoVE (see #5 below)

  allOutputDirs <- paste0(allStudyAreas, "_", climateScenario) |>
    rep(length(thisRep)) |>
    sort() |>
    paste0("_run", thisRep)

  paths <- list(
    project = workflowtools::findProjectPath(),
    cache = "cache_wbi",
    inputs = "inputs_wbi",
    outputs = file.path("outputs_wbi", studyAreaName, paste0(climateScenario, "_run", thisRep)),
    scratch = ifelse(dir.exists("/mnt/scratch"), file.path("/mnt/scratch", user, projName), "scratch")
  )
  paths$mosaics <- file.path(paths$outputs, "mosaics")
  paths$terra <- checkPath(file.path(paths$scratch, "terra", climateScenario, "00"), create = TRUE)

  ## set the max number of cores to use for parallel computations
  # options(parallelly.availableCores.fallback = 4L) ## set to limit the number of cores
  no_cores <- AGBtrends::getNumCores() ## use up to half the number cores or fallback

  terraOptions(tempdir = paths$terra, todisk = TRUE)

  ## WBI default CRS
  targetCRS <- paste("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95",
                     "+x_0=0 +y_0=0 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0") |> crs()

  ## define time intervals (year ranges between 2011-2100)
  timeint <- list(
    t1 = 1:5, t2 = 6:10, t3 = 11:15, t4 = 16:20, t5 = 21:25, t6 = 26:30,
    t7 = 31:35, t8 = 36:40, t9 = 41:45, t10 = 46:50, t11 = 51:55, t12 = 56:60,
    t13 = 61:65, t14 = 66:70, t15 = 71:75, t16 = 76:80, t17 = 81:85, t18 = 86:90
  )
  timeint_all <- timeint |> unlist() |> unname() |> list(all = _)

  years <- (timeint_all |> unlist() |> unname()) + 2010
  n_int <- length(timeint)

  # 1) data import ------------------------------------------------------------------------------

  ## previously-created ABoVE and WBI study area polygons
  agb_gpkg <- file.path("outputs", "studyArea_WBI", "ABoVE_AGB_study_area.gpkg")
  dstagnt_gpkg <- file.path("outputs", "studyArea_WBI", "ABoVE_DistAgents_study_area.gpkg")
  studyArea_gpkg <- file.path("outputs", "studyArea_WBI", "WBI_studyArea.gpkg")

  ## 1.1) Import ABoVE product tiles ------------------------------------------------------------
  # agb_tiles <- st_read(agb_gpkg, "tileset", quiet = TRUE) |>
  #   st_transform(targetCRS)
  # dist_tiles <- st_read(dstagnt_gpkg, "tileset", quiet = TRUE) |>
  #   st_transform(targetCRS)

  ## 1.2) Index WBI AGB input rasters -----------------------------------------------------------
  agbdsn <- file.path(paths$inputs, allOutputDirs, "postprocess")
  agbfiles <- list.files(agbdsn, pattern = "simulatedBiomassMap_redux_.*[.]tif$", full.names = TRUE)

  agb_out <- file.path(paths$outputs, "agb") |>
    checkPath(create = TRUE)

  ## 1.3) Index WBI disturbance history raster masks --------------------------------------------
  distdsn <- file.path(paths$inputs, allOutputDirs)
  distfiles <- list.files(distdsn, pattern = "burnMap_.*[.]tif$", full.names = TRUE)

  # all(distfiles %in% distrasters)

  ## 1.4) Import WBI study area -----------------------------------------------------------------
  wbi <- st_read(studyArea_gpkg, quiet = TRUE) |>
    st_transform(targetCRS)

  # 2) simulated cell-specific stand age --------------------------------------------------------

  agedsn <- file.path(paths$inputs, allOutputDirs, "postprocess")
  agefiles <- list.files(agedsn, pattern = "timeSinceFire_.*[.]tif$", full.names = TRUE)
  sadirs <- file.path(paths$inputs, allOutputDirs)

  age_out <- file.path(paths$outputs, "age") |>
    checkPath(create = TRUE)

  ptime <- system.time({
    ## this prepares both agb and age rasters, creating "stacks" of timeseries for each WBI prov;
    ## this works on the rasters "as-is", which includes the buffered area used for simulation.
    source("R/WBI_standAge.R")
    outfiles <- WBI_standAge(agbfiles, agefiles, sadirs, years)
    terra::tmpFiles(remove = TRUE)
  })

  ## 2.1) mask study area raster with the study area polygon to remove buffered portion

  source("R/WBI_prepRasters.R")

  agb_src <- WBI_prepRasters(
    type = "agb",
    studyAreas = allStudyAreas,
    years = years,
    paths = paths,
    dst = agb_out
  )

  age_src <- WBI_prepRasters(
    type = "age",
    studyAreas = allStudyAreas,
    years = years,
    paths = paths,
    dst = age_out
  )

  ## Create mask limiting input values to cells disturbed during time series --------------------

  dist_out <- file.path(paths$outputs, "binary_disturbed") |>
    checkPath(create = TRUE)

  no_cores <- AGBtrends::getNumCores(0.95) ## TODO: adjust as needed
  mempercore <- as.integer(terra::free_RAM() / 1024^2 / no_cores) ## TODO: propagate this; find better ram detection fun
  cl <- parallelly::makeClusterPSOCK(
    no_cores,
    default_packages = c("stringr", "terra"),
    rscript_libs = .libPaths(),
    autoStop = TRUE
  )

  parallel::clusterExport(cl, varlist = c("dist_out", "mempercore", "no_cores", "paths"),
                          envir = environment())
  parallel::clusterEvalQ(cl, {
    terraOptions(
      tempdir = paths$terra,
      memmax = min(25, mempercore),
      memfrac = 0.8 / no_cores,
      progress = 1,
      verbose = TRUE
    )

    invisible(NULL)
  })

  parallel::parLapply(cl, distfiles, function(f) {
    sa <- str_sub(basename(dirname(f)), start = 1, end = 2L)
    yr <- str_sub(basename(f), start = 9, end = 12L)
    r <- any(rast(f))
    names(r) <- yr
    writeRaster(
      r,
      filename = file.path(dist_out, paste0(sa, "_disturbed_", yr, ".tif")),
      overwrite = TRUE
    )

    invisible(NULL)
  })

  parallel::stopCluster(cl)
  terra::tmpFiles(remove = TRUE)

  ## mask study area raster with the study area polygon to remove buffered portion

  source("R/WBI_prepRasters.R")

  dist_src <- WBI_prepRasters(
    type = "disturbed",
    studyAreas = allStudyAreas,
    years = years,
    paths = paths,
    dst = dist_out
  )

  ## build raster mosaic

  dist_mosaic <- file.path(paths$mosaics) |>
    checkPath(create = TRUE)

  source("R/WBI_buildMosaics.R")

  f0a <- WBI_buildMosaics("binary_disturbed", intervals = NULL, src = dist_src, dst = dist_mosaic)

  # 3) Evaluate frequency of forested land cover types by ecozone -------------------------------

  ## TODO: calculate LCC from species layers maps for WBI
  if (FALSE) {
    tileIDs <- dir(file.path(paths$outputs, "tiles"))

    ## Create land cover mosaic
    lc_src <- list.files(file.path(paths$inputs, "ABoVE_LandCover"), full.names = TRUE)
    lc_dst <- file.path(paths$mosaics)

    ## NOTE: we want LandCover_Simplified (10 classes instead of 15)
    f0b <- AGBtrends::buildMosaics("LandCover_Simplified", intervals = timeint_all, src = lc_src, dst = lc_dst)

    ## Rasterize ecozones to fit
    ecoRast_tif <- file.path(paths$terra, "ecoRast.tif")

    ecoRast <- terra::rasterize(
      x = wbi,
      y = rast(file.path(paths$mosaics, "landcover_simplified_mosaic.tif")),
      field = "ECOZONE",
      filename = ecoRast_tif
    )

    ecoRast_lvls <- if (exists("ecoRast", .GlobalEnv)) {
      cats(ecoRast)[[1]]
    } else {
      cats(rast(ecoRast_tif))[[1]]
    }

    no_cores <- AGBtrends::getNumCores(nrow(ecoRast_lvls))
    cl <- parallelly::makeClusterPSOCK(no_cores,
                                       default_packages = c("terra"),
                                       rscript_libs = .libPaths(),
                                       autoStop = TRUE
    )

    parallel::clusterEvalQ(cl, {
      terraOptions(
        tempdir = paths$terra,
        memmax = 25,
        memfrac = 0.5 / no_cores,
        progress = 1,
        verbose = TRUE
      )

      invisible(NULL)
    })

    system.time({
      ftab <- parLapply(cl, ecoRast_lvls$value, function(ezone) {
        return(freq(mask(
          rast(file.path(paths$outputs, "landcover_simplified_mosaic.tif")),
          rast(ecoRast_tif),
          maskvalues = ezone,
          inverse = TRUE
        )))
      })
    })

    stopCluster(cl)

    names(ftab) <- ecoRast_lvls$ECOZONE
    ltab <- data.frame(layer = 1:31, year = years)
    reftab <- data.frame(
      value = 1:15,
      cat = c(
        "Evergreen Forest", "Deciduous Forest", "Mixed Forest", "Woodland",
        "Low Shrub", "Tall Shrub", "Open Shrubs", "Herbaceous", "Tussock Tundra",
        "Sparsely Vegetated", "Fen", "Bog", "Shallows/littoral", "Barren", "Water"
      )
    )

    ftab <- lapply(ftab, function(x) {
      x |>
        mutate(year = ltab$year[match(x$layer, ltab$layer)], .after = layer) |>
        mutate(ecozone = reftab$cat[match(x$value, reftab$value)], .after = value)
    })

    saveRDS(ftab, file.path(paths$outputs, "ABoVE_LandCover_freq_tables.rds"))
  }
})

# cleanup -------------------------------------------------------------------------------------
terra::tmpFiles(remove = TRUE)
