defineModule(sim, list(
  name = "AGB_analyses",
  description = paste(
    "Estimate cell-wise linear regression coefficients for ABoVE AGB 31-year time series (1984-2014)"
  ),
  keywords = "", ## TODO
  authors = c(
    person("Tyler D", "Rudolph", email = "tyler.rudolph@nrcan-rncan.gc.ca", role = c("aut")),
    person("Céline", "Boisvenue", email = "celine.boisvenue@nrcan-rncan.gc.ca", role = c("aut")),
    person("Alex M", "Chubaty", email = "achubaty@for-cast.ca", role = c("aut", "cre"))
  ),
  childModules = character(0),
  version = list(AGB_analyses = "0.0.1"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("NEWS.md", "README.md", "AGB_analyses.Rmd"),
  reqdPkgs = list("cowplot", "dplyr", "ggplot2", "ggspatial", "gridGraphics",
                  "parallel", "parallelly (>= 1.33.0)", "purrr", "sf", "stringr", "terra",
                  "PredictiveEcology/AGBtrends (>= 0.0.4)",
                  "PredictiveEcology/reproducible@development",
                  "PredictiveEcology/SpaDES.core@development (>= 1.1.0.9017)"),
  parameters = bindrows(
    defineParameter("analysisZonesType", "character", "ecozone", NA, NA,
                    paste("spatial scale at which to apply trend analyses;",
                          "one of 'ecozone', 'ecoprovince', 'ecoregion', or 'ecodistrict'.")),
    defineParameter("dataYears", "integer", 1984L:2014L, NA, NA,
                    "ABoVE AGB time series years"),
    defineParameter("nCores", "integer", min(32L, parallelly::availableCores(constraints = "connections")), NA, NA,
                    "Number of cpu threads be used where possible."),
    defineParameter("summaryIntervals", "list",
                    list(t1 = 1:5, t2 = 6:10, t3 = 11:15, t4 = 16:20, t5 = 21:25, t6 = 26:31), NA, NA,
                    paste("named list of summary intervals, with indicies specified using",
                          "`1:length(dataYears)` rather than specific years.")),
    defineParameter(".plots", "character", "screen", NA, NA,
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
    defineParameter(".seed", "list", list(), NA, NA,
                    "Named list of seeds to use for each event (names)."),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    "Should caching of events or module be used?")
  ),
  inputObjects = bindrows(
    expectsInput("analysisZones", "sf",
                 desc = "spatial units for which trend analyses will be applied and reported.",
                 sourceURL = NA),
    expectsInput("studyArea", "sf",
                 desc = paste("Polygon to use as the study area. Must be provided by the user"),
                 sourceURL = NA),
  ),
  outputObjects = bindrows(
    #createsOutput("objectName", "objectClass", "output object description", ...),
    createsOutput(objectName = NA, objectClass = NA, desc = NA)
  )
))

## event types
#   - type `init` is required for initialization

doEvent.AGB_analyses = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      # do stuff for this event
      sim <- Init(sim)

      # schedule future event(s)
      sim <- scheduleEvent(sim, start(sim), "AGB_analyses", "geoWghtReg")
      sim <- scheduleEvent(sim, start(sim), "AGB_analyses", "buildVRTs")
      sim <- scheduleEvent(sim, start(sim), "AGB_analyses", "slopesPerTime")
      sim <- scheduleEvent(sim, start(sim), "AGB_analyses", "createMosaicRasts")
      sim <- scheduleEvent(sim, start(sim), "AGB_analyses", "slopesByAgeTime")
      sim <- scheduleEvent(sim, start(sim), "AGB_analyses", "SAbyZOI")
      sim <- scheduleEvent(sim, start(sim), "AGB_analyses", "statsByZOI")
      sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "AGB_analyses", "plotting")
    },
    geoWghtReg = {
      ## 1.1.1) Derive slope of numerical vector across a time series
      f1 <- AGBtrends::gwr(mod$tile_folders, type = "slopes", cores = P(sim)$nCores)
      # sim <- registerOutputs(sim, f1) ## TODO: enable once implement in SpaDES.core

      ## 1.1.2) stock number of non-NA values for subsequent weighted standard deviation
      f2 <- AGBtrends::gwr(mod$tile_folders, type = "sample_size", cores = P(sim)$nCores)
      # sim <- registerOutputs(sim, f2) ## TODO: enable once implement in SpaDES.core
    },
    buildVRTs = {
      outputDir <- outputPath(sim)
      scratchDir <- scratchPath(sim)
      tileDirs <- mod$tile_folders

      tifs1 <- AGBtrends::buildMosaics(type = "slopes", intervals = c(all = 1:31), paths = paths)
      tifs2 <- AGBtrends::buildMosaics(type = "sample_size", intervals = c(all = 1:31), paths = paths)

      # sim <- registerOutputs(sim, c(tifs1, tifs2)) ## TODO: enable once implement in SpaDES.core
    },
    slopesPerTime = {
      ## 2.1.1) calculate local slope coefficient for specified time interval
      f1 <- AGBtrends::gwrt(mod$tile_folders, type = "slopes", cores = P(sim)$nCores, intervals = P(sim)$summaryIntervals)
      # sim <- registerOutputs(sim, f1) ## TODO: enable once implement in SpaDES.core

      ## 2.1.2) stock number of non-NA values for subsequent weighted standard deviation
      f2 <- AGBtrends::gwrt(mod$tile_folders, type = "sample_size", cores = P(sim)$nCores, intervals = P(sim)$summaryIntervals)
      # sim <- registerOutputs(sim, f2) ## TODO: enable once implement in SpaDES.core
    },
    createMosaicRasts = {
      intervals <- P(sim)$summaryIntervals
      ncores <- length(P(sim)$summaryIntervals)
      paths <- list(
        outputs = outputPath(sim),
        scratch = scratchPath(sim),
        tiles = mod$tile_folders
      )

      f3a <- buildMosaics(type = "slopes", intervals = intervals, paths = paths)
      f3b <- buildMosaics(type = "sample_size", intervals = intervals, paths = paths)

      # sim <- registerOutputs(sim, c(f3a, f3b)) ## TODO: enable once implement in SpaDES.core
    },
    slopesByAgeTime = {
      sim <- groupSlopes(sim)
    },
    SAbyZOI = {
      sim <- studyAreaByZOI(sim)
    },
    statsByZOI = {
      sim <- summarizeZOI(sim)
    },
    plotting = {
      plotAll(sim)
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

  # ! ----- EDIT BELOW ----- ! #

  mod$cleanTilePath <- checkPath(file.path(mod$dPath, "clean", "tiled"), create = TRUE)
  mod$tile_folders <- sort(list.files(mod$cleanTilePath, pattern = "Bh", full.names = TRUE))

  mod$targetCRS <- AGBtrends::Canada_Albers_Equal_Area_Conic

  # ! ----- STOP EDITING ----- ! #

  return(invisible(sim))
}

groupSlopes <- function(sim) {
  cacheDir <- cachePath(sim)
  cleanTilePath <- mod$cleanTilePath
  ncores <- length(P(sim)$summaryIntervals)
  outputDir <- outputPath(sim)
  scratchDir <- scratchPath(sim)
  terraDir <- terraPath(sim)
  tileDirs <- mod$tile_folders
  timeint <- P(sim)$summaryIntervals
browser() ## TODO: replace this chunk with AGBtrends::buildMosaics()
  cl <- parallelly::makeClusterPSOCK(ncores,
                                     default_packages = c("sf", "terra"),
                                     rscript_libs = .libPaths(),
                                     autoStop = TRUE)
  on.exit(stopCluster(cl), add = TRUE)

  parallel::clusterExport(cl, varlist = c("cacheDir", "ncores", "terraDir", "tileDirs"), envir = environment())
  parallel::clusterEvalQ(cl, {
    terraOptions(tempdir = terraDir,
                 memmax = 25,
                 memfrac = 0.6 / ncores,
                 progress = 1,
                 verbose = TRUE)
  })

  fout <- parallel::parLapply(cl, seq(length(timeint)), function(i) {
    ## 3.1) Build virtual raster of estimated stand age at time 0 (i.e. 1984)
    sf::gdal_utils(
      util = "buildvrt",
      source = sapply(paths$tiles, function(dsn) unname(file.path(dsn, list.files(dsn, pattern = "rage")))),
      destination = file.path(paths$cache, paste0("AGB_age_mosaic_t", i, ".vrt")),
      ## !! MODIFY if doing alternative time steps is a DESIRED functionality (e.g. 10-year, etc.)
      options = c("-b", c(1, 6, 11, 16, 21, 26)[i]) # time 1 = 1984 etc.
    )

    ## 3.2) Write to raster mosaic (stand age, kNN 2020)
    f_ageMosaic <- file.path(paths$outputs, "mosaics", paste0("AGB_age_mosaic_t", i, ".tif"))
    sf::gdal_utils(
      util = "warp",
      source = file.path(paths$cache, paste0("AGB_age_mosaic_t", i, ".vrt")),
      destination = f_ageMosaic
    )

    ## 3.3) Group into 5 discrete age classes
    ages_from <- c(25, 50, 80, 125, 500)
    ages_to <- c(24, 49, 79, 124, 500)

    ageRast <- classify(
      rast(f_ageMosaic),
      rcl = cbind(from = ages_from, to = ages_to, becomes = 1L:5L),
      right = FALSE, others = NA_integer_
    )

    names(ageRast) <- "ageClass"
    levels(ageRast) <- data.frame(value = 1:5, ageClass = paste0(ages_from, "-", to = ages_to))

    f_ageMosaicClass <- file.path(paths$outputs, "mosaics", paste0("AGB_age_mosaic_classes_t", i, ".tif"))
    writeRaster(ageRast, f_ageMosaicClass, overwrite = TRUE)

    return(f_ageMosaicClass)
  }) |>
    unlist()

  # sim <- registerOutputs(sim, fout) ## TODO: enable once implement in SpaDES.core

  return(invisble(sim))
}

studyAreaByZOI <- function(sim) {
  analysisZones <- sim$analysisZones
  analysisZonesType <- P(sim)$analysisZonesType
  cleanTilePath <- mod$cleanTilePath
  ncores <- length(P(sim)$summaryIntervals)
  terraDir <- terraPath(sim)
browser()
  cl <- parallelly::makeClusterPSOCK(ncores,
                                     default_packages = c("AGBtrends", "sf", "terra"),
                                     rscript_libs = .libPaths(),
                                     autoStop = TRUE)
  on.exit(stopCluster(cl), add = TRUE)

  parallel::clusterExport(cl, varlist = c("ncores"), envir = environment())
  parallel::clusterEvalQ(cl, {
    terraOptions(tempdir = terraDir, memfrac = 0.5 / ncores)
  })

  parallel::parLapply(cl, seq_len(length(P(sim)$summaryIntervals)), function(i) {
    prepZones(
      zoi = analysisZones,
      field = toupper(analysisZonesType),
      ageClass = rast(file.path(cleanTilePath, paste0("AGB_age_mosaic_classes_t", i, ".tif"))),
      fileID = paste0("WBI_ecozone_t", i),
      ow = TRUE
    )

    return(invisible(NULL))
  })

  return(invisible(sim))
}

summarizeZOI <- function(sim) {
  ncores <- length(P(sim)$summaryIntervals) + 1 ## includes full time series, not just summary intervals
  terraDir <- terraPath(sim)
browser()
  ## NOTE: age at beginning of the 31 year time series (1984-2014) is identical
  ## to age at beginning of 't1' time interval (i.e. 1984-1988)

  files <- list(
    list.files(outputPath(sim), pattern = "slope_mosaic", full.names = TRUE),
    list.files(outputPath(sim), pattern = "sample_size", full.names = TRUE),
    list.files(outputPath(sim), pattern = "ZOIxageClass_WBI_ecozone", full.names = TRUE),
    list.files(outputPath(sim), pattern = "ZOIxageClass_WBI_ecoregion", full.names = TRUE),
    list.files(outputPath(sim), pattern = "ZOIxageClass_WBI_ecoprovince", full.names = TRUE)
  )
  omit <- c(2, 4, 6, 8, 10, 12) ## TODO: omit aux files by grepping tifs only above
  irast <- list(
    slope = files[[1]],
    w = files[[2]],
    ## these last values '2' (below) refer to ageClass at 'time 0' (i.e. 1984) used
    ## for the complete time series slope raster mosaic stats assessment.
    ## this should be '1', but currently set to 6 years in b/c disturbed pixels are only known as of 1987
    ecozone = files[[3]][-omit][c(1:6, 2)],
    ecoregion = files[[4]][-omit][c(1:6, 2)],
    ecoprovince = files[[5]][-omit][c(1:6, 2)]
  )

  cl <- parallelly::makeClusterPSOCK(ncores,
                                     default_packages = c("AGBtrends", "dplyr", "terra"),
                                     rscript_libs = .libPaths(),
                                     autoStop = TRUE)
  on.exit(stopCluster(cl), add = TRUE)

  parallel::clusterExport(cl, varlist = c("irast", "ncores"), envir = environment())
  parallel::clusterEvalQ(cl, {
    terraOptions(tempdir = terraDir, memfrac = 0.5 / ncores)
  })

  parallel::parLapply(cl, 1:7, function(i, svar = "ecozone", maskRaster = NULL) {
    ## TODO: use maskRaster file name to qualify fileID writeRaster tag
    if (i == 7) {
      fileID <- paste0("WBI_", svar)
      # fileID <- paste0("WBI_distMask_, svar)
    } else {
      fileID <- paste0("WBI_", svar, "_t", i)
      # fileID <- paste0("WBI_distMask_, svar, "_t, i)
    }

    zoneStats(slopeRaster = rast(irast$slope[i]),
              weightRaster = rast(irast$w[i]),
              zoneRaster = rast(irast[[svar]][i]),
              ## maskRaster arg can be either e.g. was it disturbed, is it forested or both?
              maskRaster = maskRaster, # e.g. use rast("inputs/raw/ABoVE_ForestDisturbance_Agents/binary_disturbed_mosaic.tif") for pixels disturbed over course of time series (according to ABoVE)
              fileID = fileID)

    return(invisible(NULL))
  })

  return(invisible(sim))
}

plotAll <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #

  ## 2) visual examination of results ------------------------------------------------------------
  ## TODO

  # par(mfrow=c(3,2))
  # for(tp in names(timeint)) {
  #   plot(rast(paste0("outputs/AGB_slope_mosaic_", tp, ".tif")), main=tp)
  # }
  #
  # sapply(names(timeint), function(tp) digest::digest(paste0("outputs/AGB_slope_mosaic_", tp, ".tif"), algo="xxhash64"))
  # sapply(names(timeint), function(tp) file.size(paste0("outputs/AGB_slope_mosaic_", tp, ".tif")))
  #
  # ## write to PNG !!
  # par(mfrow=c(2,3))
  # hist(rast(paste0("outputs/AGB_slope_mosaic_t1.tif")), main="t1")
  # hist(rast(paste0("outputs/AGB_slope_mosaic_t2.tif")), main="t2")
  # hist(rast(paste0("outputs/AGB_slope_mosaic_t3.tif")), main="t3")
  # hist(rast(paste0("outputs/AGB_slope_mosaic_t4.tif")), main="t4")
  # hist(rast(paste0("outputs/AGB_slope_mosaic_t5.tif")), main="t5")
  # hist(rast(paste0("outputs/AGB_slope_mosaic_t6.tif")), main="t6")

  #Plots(sampleData, fn = ggplotFn)


  ## 6) diagnostic plots -------------------------------------------------------------------------

  ## 6.1) range, mode and mean of AGB values by ageClass, year 2000


  ## 7) plot differences -------------------------------------------------------------------------

  fout <- file.path(outputPath(sim), "figures", c(
    paste0("AGB_temporal_trends_x_ECOZONE_x_ageClass_", Sys.Date(), ".pdf"), ## TODO: use .png
    paste0("AGB_temporal_trends_x_ageClass_x_ECOZONE_", Sys.Date(), ".pdf"), ## TODO: use .png
    paste0("AGB_temporal_trends_x_ECOZONE_distMask_", Sys.Date(), ".pdf") ## TODO: use .png
  ))
  f2p <- list.files(outputPath(sim), pattern = "zoneStats_summary_WBI_ecozone_", full.names = TRUE) ## TODO: use fs pkg
  f2pd <- list.files(outputPath(sim), pattern = "zoneStats_summary_WBI_distMask_ecozone_", full.names = TRUE) ## TODO: use fs pkg
  f2pd2 <- list.files(outputPath(sim), pattern = "WBI_distMask_ecozone", full.names = TRUE) ## TODO: use fs pkg

  ## 7a) without disturbance mask

  ## TODO: update based on the updated script

  ## i=1 corresponds to 31-year time series, i=2 corresponds to time interval t1 (1984-1988), etc.
  png() ## TODO
  plotZoneStats(file2plot = file.path(outputPath(sim), paste0("zoneStats_summary_WBI_ecozone.rds")))
  dev.off()

  ## x Ecozone x ageClass
  pdf(file = fout[1]) ## TODO: use png() with ptype = 1
  lapply(plotZoneStatsIntervals(files2plot = f2p, weighted = TRUE, xVar = "tp", groupVar = "ageClass", ptype = 2), plot)
  dev.off()

  ## x ageClass x Ecozone
  ## TODO: this one crashes; not sure why yet
  pdf(file = fout[2]) ## TODO: use png() with ptype = 1
  lapply(plotZoneStatsIntervals(files2plot = f2p, weighted = TRUE, xVar = "tp", catVar = "ageClass", groupVar = "ECOZONE", ptype = 2), plot)
  dev.off()

  ## 7b) with disturbance mask
  ## TODO: make this a box plot per ecozone
  png() ## TODO
  plotZoneStats(file2plot = file.path(outputPath(sim), paste0("zoneStats_summary_WBI_distMask_ecozone.rds")))
  dev.off()

  ## x ageClass x Ecozone
  pdf(file = fout[3]) ## TODO: use png() with ptype = 1
  lapply(plotZoneStatsIntervals(files2plot = f2pd, weighted = TRUE, xVar = "tp", catVar = "ageClass", groupVar = "ECOZONE", ptype = 2), plot)
  dev.off()

  gp <- plotZoneStatsIntervals(files2plot = f2pd2) ## TODO: save as png

  #Plots(sampleData, fn = ggplotFn)

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  cacheTags <- c(currentModule(sim), "function:.inputObjects")
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")

  # ! ----- EDIT BELOW ----- ! #

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}
