defineModule(sim, list(
  name = "AGB_analyses",
  description = paste(
    "Estimate cell-wise linear regression coefficients for ABoVE AGB 31-year time series (1984-2014)"
  ),
  keywords = "", ## TODO
  authors = c(
    person("Tyler D", "Rudolph", email = "tyler.rudolph@nrcan-rncan.gc.ca", role = c("aut", "cre")),
    person("Alex M", "Chubaty", email = "achubaty@for-cast.ca", role = c("ctb"))
  ),
  childModules = character(0),
  version = list(AGB_analyses = "0.0.1"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("NEWS.md", "README.md", "AGB_analyses.Rmd"),
  reqdPkgs = list("dplyr", "geodata", "gdalUtilities", "ggplot2", "ggspatial",
                  "parallel", "parallelly (>= 1.33.0)", "purrr", "sf", "stringr", "terra",
                  "PredictiveEcology/reproducible@development",
                  "PredictiveEcology/SpaDES.core@development (>= 1.1.0.9017)"),
  parameters = bindrows(
    defineParameter("analysisZonesType", "character", "ecozone", NA, NA,
                    paste("spatial scale at which to apply trend analyses;",
                          "one of 'ecozone', 'ecoprovince', 'ecoregion', or 'ecodistrict'.")),
    defineParameter("dataYears", "integer", 1984L:2014L, NA, NA,
                    "ABoVE AGB time series years"),
    defineParameter("nCores", "integer", max(32L, parallelly::availableCores(constraints = "connections")), NA, NA,
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
      sim <- GWR(sim)
    },
    buildVRTs = {
      sim <- VRT(sim)
    },
    slopesPerTime = {
      sim <- calcSlopes(sim)
    },
    createMosaicRasts = {
      sim <- buildMosaics(sim)
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

  # ! ----- STOP EDITING ----- ! #

  return(invisible(sim))
}

GWR <- function(sim) {
  odt <- terraOptions(print = FALSE)[["datatype"]]
  on.exit(terraOptions(datatype = odt), add = TRUE)

  ncores <- P(sim)$nCores

  tileDirs <- mod$tile_folders
  sapply(seq_along(tileDirs), function(i) {
    ## 1.1.1) Derive slope of numerical vector across a time series
    terraOptions(datatype = "FLT4S")
    app(rast(file.path(tileDirs[i], list.files(tileDirs[i], pattern = "ragb"))),
        fun = function(x, ff) ff(x),
        cores = ncores,
        ff = ts_slope, ## defined in R/helpers.R
        filename = file.path(tileDirs[i], paste0("agb_slopes_", str_sub(tileDirs[i], start = -7L), ".tif")),
        overwrite = TRUE)

    ## 1.1.2) stock number of non-NA values for subsequent weighted standard deviation
    terraOptions(datatype = "INT1U")
    app(rast(file.path(tileDirs[i], list.files(tileDirs[i], pattern = "ragb"))),
        fun = function(x, ff) ff(x),
        cores = ncores,
        ff = ts_nsamp, ## defined in R/helpers.R
        filename = file.path(tileDirs[i], paste0("agb_sample_size_", str_sub(tileDirs[i], start = -7L), ".tif")),
        overwrite = TRUE)

    return(invisible(NULL))
  })

  return(invisible(sim))
}

VRT <- function(sim) {
  tileDirs <- mod$tile_folders
  vrts <- file.path(cachePath(sim), c("AGB_slope_mosaic.vrt", "AGB_sample_size_mosaic.vrt")) ## TODO: why cache, not tmp dir??
  tifs <- file.path(outputPath(sim), c("AGB_slope_mosaic.tif", "AGB_sample_size_mosaic.tif"))
browser()
  ## 1.2.1) Build virtual rasters
  gdalbuildvrt(gdalfile = unname(sapply(tileDirs, function(dsn) file.path(dsn, list.files(dsn, pattern = "agb_slopes_Bh")))),
               output.vrt = vrts[1])
  gdalbuildvrt(gdalfile = unname(sapply(tileDirs, function(dsn) file.path(dsn, list.files(dsn, pattern = "agb_sample_size_Bh")))),
               output.vrt = vrts[2])

  ## 1.2.2) Write to raster mosaics
  gdalwarp(srcfile = vrts[1], dstfile = tifs[1], overwrite = TRUE)
  gdalwarp(srcfile = vrts[2], dstfile = tifs[2], overwrite = TRUE)

  return(invisible(sim))
}

calcSlopes <- function(sim) {
  ncores <- P(sim)$nCores
  tileDirs <- mod$tile_folders
  timeint <- P(sim)$summaryIntervals

  odt <- terraOptions(print = FALSE)[["datatype"]]
  on.exit(terraOptions(datatype = odt), add = TRUE)

  for (i in seq_along(tileDirs)) {
    message(paste0(i, "/", length(tileDirs), "\n")) ## simple progress updates

    tile_folder <- tileDirs[i]

    sapply(1:length(timeint), function(timestep) {
      ## 2.1.1) calculate local slope coefficient for specified time interval
      terraOptions(datatype = "FLT4S")
      app(rast(file.path(tileDirs[i], list.files(tileDirs[i], pattern = "ragb")))[[timeint[[timestep]]]],
          function(x, ff) ff(x),
          cores = ncores,
          ff = ts_slope, ## defined in R/helpers.R
          filename = file.path(tileDirs[i], paste0("agb_slopes_", names(timeint)[timestep], "_", str_sub(tileDirs[i], start = -7L), ".tif")),
          overwrite = TRUE)

      ## 2.1.2) stock number of non-NA values for subsequent weighted standard deviation
      terraOptions(datatype = "INT1U")
      app(rast(file.path(tileDirs[i], list.files(tileDirs[i], pattern = "ragb")))[[timeint[[timestep]]]],
          fun = function(x, ff) ff(x),
          cores = ncores,
          ff = ts_nsamp, ## defined in R/helpers.R
          filename = file.path(tileDirs[i], paste0("agb_sample_size_", names(timeint)[timestep], "_", str_sub(tileDirs[i], start = -7L), ".tif")),
          overwrite = TRUE)

      return(invisible(NULL))
    })
  }

  return(invisible(sim))
}

buildMosaics <- function(sim) {
  ncores <- length(P(sim)$summaryIntervals)
  outputDir <- outputPath(sim)
  scratchDir <- scratchPath(sim)
  tileDirs <- mod$tile_folders
browser()
  cl <- parallelly::makeClusterPSOCK(ncores,
                                     default_packages = c("terra", "gdalUtilities", "stringr"),
                                     rscript_libs = .libPaths(),
                                     autoStop = TRUE)
  on.exit(stopCluster(cl), add = TRUE)

  clusterExport(cl, varlist = c("outputDir", "scratchDir", "tileDirs"))

  parallel::clusterEvalQ(cl, {
    terraOptions(tempdir = scratchDir,
                 memmax = 25,
                 memfrac = 0.6,
                 progress = 1,
                 verbose = TRUE) ## TODO: P(sim)$verbose
  })

  parLapply(cl, names(timeint), function(tp) {
    td <- terraOptions(print = FALSE)[["tempdir"]]
    od <- outputDir

    ## 2.2.1) Build virtual rasters
    flist <- unname(sapply(tileDirs, function(dsn) file.path(dsn, list.files(dsn, pattern = tp))))

    gdalbuildvrt(gdalfile = flist[str_detect(flist, "slope")],
                 output.vrt = paste0(td, "AGB_slope_mosaic_", tp, ".vrt"))

    gdalbuildvrt(gdalfile = flist[str_detect(flist, 'sample_size')],
                 output.vrt = paste0(td, "AGB_sample_size_mosaic_", tp, ".vrt"))

    ## 2.2.2) Write to raster mosaics
    gdalwarp(srcfile = file.path(td, paste0("AGB_slope_mosaic_", tp, ".vrt")),
             dstfile = file.path(od, paste0("AGB_slope_mosaic_", tp, ".tif")))

    gdalwarp(srcfile = file.path(td, paste0("AGB_sample_size_mosaic_", tp, ".vrt")),
             dstfile = file.path(od, paste0("AGB_sample_size_mosaic_", tp, ".tif")))

    return(invisible(NULL))
  })

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
browser()
  cl <- parallelly::makeClusterPSOCK(ncores,
                                     default_packages = c("terra", "gdalUtilities"),
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

  parallel::parLapply(cl, seq_len(length(P(sim)$summaryIntervals)), function(i) {
    ## 3.1) Build virtual raster of estimated stand age at time 0 (i.e. 1984)
    gdalbuildvrt(
      gdalfile = sapply(tile_folders, function(dsn) unname(file.path(dsn, list.files(dsn, pattern = "rage")))),
      output.vrt = file.path(cacheDir, paste0("AGB_age_mosaic_t", i, ".vrt")),
      b = c(1, 6, 11, 16, 21, 26)[i], # time 1 = 1984 etc. ## TODO: use sapply(timeint, `[`, 1)
      overwrite = TRUE
    )

    ## 3.2) Write to raster mosaic (stand age, kNN 2020)
    gdalwarp(srcfile = file.path(cacheDir, paste0("AGB_age_mosaic_t", i, ".vrt")),
             dstfile = file.path(cleanTilePath, paste0("AGB_age_mosaic_t", i, ".tif")),
             overwrite = TRUE)

    ## 3.3) Group into 5 discrete age classes
    RCL <- cbind(from = c(0, 25, 50, 80, 125), to = c(25, 50, 80, 125, 500), becomes = 1L:5L)
    ageRast <- classify(rast(file.path(cleanTilePath, paste0("AGB_age_mosaic_t", i, ".tif"))),
                        rcl = RCL, right = FALSE, others = NA_integer_)
    names(ageRast) <- "ageClass"
    levels(ageRast) <- data.frame(value = RCL$becomes, ageClass = paste0(RCL$from, '-', RCL$to))
    writeRaster(ageRast, file.path(cleanTilePath, paste0("AGB_age_mosaic_classes_t", i, ".tif")), overwrite = TRUE)

    return(invisible(NULL))
  })

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
                                     default_packages = c("terra", "gdalUtilities"),
                                     rscript_libs = .libPaths(),
                                     autoStop = TRUE)
  on.exit(stopCluster(cl), add = TRUE)

  parallel::clusterExport(cl, varlist = c("ncores", "prepZones"), envir = environment())
  parallel::clusterEvalQ(cl, {
    terraOptions(tempdir = terraDir, memfrac = 0.5 / ncores)
  })

  parallel::parLapply(cl, seq_len(length(P(sim)$summaryIntervals)), function(i) {
    prepZones(
      zoi = analysisZones,
      field = toupper(analysisZonesType),
      ageClass = rast(file.path(cleanTilePath, paste0("AGB_age_mosaic_classes_t", i, ".tif"))),
      file.id = paste0("WBI_ecozone_t", i),
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
                                     default_packages = c("terra", "gdalUtilities", "dplyr"),
                                     rscript_libs = .libPaths(),
                                     autoStop = TRUE)
  on.exit(stopCluster(cl), add = TRUE)

  parallel::clusterExport(cl, varlist = c("irast", "ncores", "zoneStats"), envir = environment())
  parallel::clusterEvalQ(cl, {
    terraOptions(tempdir = terraDir, memfrac = 0.5 / ncores)
  })

  parallel::parLapply(cl, 1:7, function(i, svar = 'ecozone', maskRaster = NULL) {
    ## TODO: use maskRaster file name to qualify file.id writeRaster tag
    if (i == 7) {
      file.id <- paste0("WBI_", svar)
      # file.id <- paste0("WBI_distMask_, svar)
    } else {
      file.id <- paste0("WBI_", svar, "_t", i)
      # file.id <- paste0("WBI_distMask_, svar, "_t, i)
    }

    zoneStats(slopeRaster = rast(irast$slope[i]),
              weightRaster = rast(irast$w[i]),
              zoneRaster = rast(irast[[svar]][i]),
              ## maskRaster arg can be either e.g. was it disturbed, is it forested or both?
              maskRaster = maskRaster, # e.g. use rast('inputs/raw/ABoVE_ForestDisturbance_Agents/binary_disturbed_mosaic.tif') for pixels disturbed over course of time series (according to ABoVE)
              file.id = file.id)

    return(invisible(NULL))
  })

  return(invisible(sim))
}

plotAll <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  ## TODO

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

  #Plots(sampleData, fn = ggplotFn)

  fout <- file.path(outputPath(sim), "figures", c(
    paste0("AGB_temporal_trends_x_ECOZONE_x_ageClass_", Sys.Date(), ".pdf"), ## TODO: use .png
    paste0("AGB_temporal_trends_x_ageClass_x_ECOZONE_", Sys.Date(), ".pdf"), ## TODO: use .png
    paste0("AGB_temporal_trends_x_ECOZONE_distMask_", Sys.Date(), ".pdf") ## TODO: use .png
  ))
  f2p <- list.files(outputPath(sim), pattern = "zoneStats_summary_WBI_ecozone_", full.names = TRUE) ## TODO: use fs pkg
  f2pd <- list.files(outputPath(sim), pattern = "zoneStats_summary_WBI_distMask_ecozone_", full.names = TRUE) ## TODO: use fs pkg
  f2pd2 <- list.files(outputPath(sim), pattern = "WBI_distMask_ecozone", full.names = TRUE) ## TODO: use fs pkg

  ## 6a) without disturbance mask

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

  ## 6b) with disturbance mask
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

ggplotFn <- function(data, ...) {
  ggplot(data, aes(TheSample)) +
    geom_histogram(...)
}

.inputObjects <- function(sim) {
  cacheTags <- c(currentModule(sim), "function:.inputObjects")
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")

  # ! ----- EDIT BELOW ----- ! #

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}
