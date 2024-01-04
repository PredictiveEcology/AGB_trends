## Date Created: Nov 30, 2022
## Auteurs: Tyler Rudolph, biologist M.Sc., CFS/NRCAN, Trade, Economics & Industry Branch
##          Alex M. Chubaty, PhD, FOR-CAST Research & Analytics
##
## Name of script : "01_ABoVE_AGB_time_series_analysis.R"
##
## Description :
##   R script serving to estimate cell-wise linear regression coefficients for
##   ABoVE AGB 31-year time series (1984-2014)

# package installation and loading ------------------------------------------------------------
Require::Require(c("dplyr", "ggplot2", "parallelly", "reproducible", "sf", "stringr", "terra",
                   "PredictiveEcology/AGBtrends"), upgrade = FALSE)

# global parameters for project setup ---------------------------------------------------------
projName <- "AGB_trends"
studyAreaName <- "studyArea_WBI"
user <- Sys.info()[["user"]]

paths <- list(
  project = getwd(),
  cache = "cache",
  inputs = "inputs",
  outputs = file.path("outputs", studyAreaName),
  scratch = ifelse(dir.exists("/mnt/scratch"), file.path("/mnt/scratch", user, projName), "scratch")
)
paths$terra <- checkPath(file.path(paths$scratch, "terra"), create = TRUE)

no_cores <- min(parallel::detectCores() / 2, 32L)

terraOptions(tempdir = paths$terra, todisk = TRUE)

file.remove(file.path(paths$terra, list.files(paths$terra)))

oldTmpDir <- tempdir()
newTmpDir <- file.path(paths$scratch, "tmp")
if (!dir.exists(newTmpDir)) dir.create(newTmpDir, recursive = TRUE)
newTmpDir <- tools::file_path_as_absolute(newTmpDir)
Sys.setenv(TMPDIR = newTmpDir)
unlink(oldTmpDir, recursive = TRUE)
tempdir(check = TRUE)

terraOptions(
  tempdir = paths$terra,
  memmax = 25,
  memfrac = 0.8,
  progress = 1,
  verbose = TRUE
)

tile_folders <- sort(fs::dir_ls("inputs/clean/tiled", regexp = "Bh", type = "directory"))

## 1.1) Estimate cell-wise linear regression coefficients for undisrupted time series ---------
## aka "local" or "geographically weighted regression (GWR)"

f1 <- AGBtrends::gwr(tile_folders, type = "slope", cores = length(tile_folders))
f2 <- AGBtrends::gwr(tile_folders, type = "nsamp", cores = length(tile_folders))

## 1.2) Combine tiled slope rasters into unified mosaics --------------------------------------

## TODO: why cache, not tmp dir??
vrts <- file.path(cachePath(sim), c("AGB_slope_mosaic.vrt", "AGB_sample_size_mosaic.vrt"))
tifs <- file.path(outputPath(sim), c("AGB_slope_mosaic.tif", "AGB_sample_size_mosaic.tif"))

agb_slopes_Bh <- vapply(tile_folders, function(dsn) {
  fs::dir_ls(dsn, regexp = "agb_slopes_Bh")
}, character(1)) |>
  unname()

agb_sample_size_Bh <- vapply(tile_folders, function(dsn) {
  fs::dir_ls(dsn, regexp = "agb_sample_size_Bh")
}, character(1)) |>
  unname()

## 1.2.1) Build virtual rasters
sf::gdal_utils(util = "buildvrt", source = agb_slopes_Bh, destination = vrts[1])
sf::gdal_utils(util = "buildvrt", source = agb_sample_size_Bh, destination = vrts[2])

## 1.2.2) Write to raster mosaics
sf::gdal_utils(utils = "warp", source = vrts[1], destination = tifs[1])
sf::gdal_utils(utils = "warp", source = vrts[2], destination = tifs[2])

f3 <- tifs

## 2.1) Calculate cell-specific slopes per 5-year time interval (n=6) -------------------------

## define time intervals (year ranges between 1984-2014)
timeint <- list(t1 = 1:5, t2 = 6:10, t3 = 11:15, t4 = 16:20, t5 = 21:25, t6 = 26:31)

### 2.1.1) calculate local slope coefficient for specified time interval ----------------------
f1 <- AGBtrends::gwrt(tile_folders, type = "slope", cores = no_cores, intervals = intervals)

### 2.1.2) stock number of non-NA values for subsequent weighted standard deviation -----------
f2 <- AGBtrends::gwrt(tile_folders, type = "nsamp", cores = no_cores, intervals = intervals)

## 2.2) Combine tiled slope rasters into numerous unified mosaics -----------------------------

no_cores <- length(timeint)
cl <- parallelly::makeClusterPSOCK(no_cores,
  default_packages = c("sf", "stringr", "terra"),
  rscript_libs = .libPaths(),
  autoStop = TRUE
)
clusterExport(cl, varlist = c("tile_folders"))
parallel::clusterEvalQ(cl, {
  terraOptions(
    tempdir = paths$terra,
    memmax = 25,
    memfrac = 0.6,
    progress = 1,
    verbose = TRUE
  )
})

parLapply(cl, names(timeint), function(tp) {
  ## 2.2.1) Build virtual rasters
  flist <- unname(sapply(tile_folders, function(dsn) file.path(dsn, list.files(dsn, pattern = tp))))

  gdalbuildvrt( ## TODO: use sf gdal utils
    gdalfile = flist[str_detect(flist, "slope")],
    output.vrt = file.path(terraOptions(print = FALSE)$tempdir, paste0("AGB_slope_mosaic_", tp, ".vrt"))
  )

  gdalbuildvrt( ## TODO: use sf gdal utils
    gdalfile = flist[str_detect(flist, "sample_size")],
    output.vrt = file.path(terraOptions(print = FALSE)$tempdir, paste0("AGB_sample_size_mosaic_", tp, ".vrt"))
  )

  ## 2.2.2) Write to raster mosaics
  gdalwarp( ## TODO: use sf gdal utils
    srcfile = file.path(terraOptions(print = FALSE)$tempdir, paste0("AGB_slope_mosaic_", tp, ".vrt")),
    dstfile = file.path(paths$outputs, paste0("AGB_slope_mosaic_", tp, ".tif"))
  )

  gdalwarp( ## TODO: use sf gdal utils
    srcfile = file.path(terraOptions(print = FALSE)$tempdir, paste0("AGB_sample_size_mosaic_", tp, ".vrt")),
    dstfile = file.path(paths$outputs, paste0("AGB_sample_size_mosaic_", tp, ".tif"))
  )

  return(invisible(NULL))
}) # 23 min

stopCluster(cl)

# Visual examination of results ---------------------------------------------------------------
## TODO: finesse and write to png
#
# par(mfrow = c(3, 2))
# for(tp in names(timeint)) {
#   plot(rast(file.path(paths$outputs, paste0("AGB_slope_mosaic_", tp, ".tif"))), main = tp)
# }
#
# sapply(names(timeint), function(tp) digest::digest(file.path(paths$outputs, paste0("AGB_slope_mosaic_", tp, ".tif")), algo = "xxhash64"))
# sapply(names(timeint), function(tp) file.size(file.path(paths$outputs, paste0("AGB_slope_mosaic_", tp, ".tif"))))
#
# ## write to PNG !!
# par(mfrow = c(2, 3))
# hist(rast(file.path(paths$outputs, paste0("AGB_slope_mosaic_t1.tif"))), main = "t1")
# hist(rast(file.path(paths$outputs, paste0("AGB_slope_mosaic_t2.tif"))), main = "t2")
# hist(rast(file.path(paths$outputs, paste0("AGB_slope_mosaic_t3.tif"))), main = "t3")
# hist(rast(file.path(paths$outputs, paste0("AGB_slope_mosaic_t4.tif"))), main = "t4")
# hist(rast(file.path(paths$outputs, paste0("AGB_slope_mosaic_t5.tif"))), main = "t5")
# hist(rast(file.path(paths$outputs, paste0("AGB_slope_mosaic_t6.tif"))), main = "t6")

# 3) Group slopes by age at time x ------------------------------------------------------------
##    (band argument determines reference layer/year),
##    effectively masking out pixels disturbed mid-time series
no_cores <- min(parallel::detectCores() / 2, 6L)
cl <- parallelly::makeClusterPSOCK(
  no_cores,
  default_packages = c("sf", "terra"),
  rscript_libs = .libPaths(),
  autoStop = TRUE
)

parallel::clusterExport(cl, varlist = c("tile_folders", "no_cores"), envir = environment())
parallel::clusterEvalQ(cl, {
  terraOptions(
    tempdir = paths$terra,
    memmax = 25,
    memfrac = 0.6 / no_cores,
    progress = 1,
    verbose = TRUE
  )
})

parallel::parLapply(cl, 1:6, function(i) {
  ## 3.1) Build virtual raster of estimated stand age at time 0 (i.e. 1984)
  gdalbuildvrt( ## TODO: use sf gdal utils
    gdalfile = sapply(tile_folders, function(dsn) unname(file.path(dsn, list.files(dsn, pattern = "rage")))),
    output.vrt = file.path(paths$cache, paste0("AGB_age_mosaic_t", i, ".vrt")),
    ## !! MODIFY if doing alternative time steps is a DESIRED functionality (e.g. 10-year, etc.)
    b = c(1, 6, 11, 16, 21, 26)[i], # time 1 = 1984 etc.
    overwrite = TRUE
  )

  ## 3.2) Write to raster mosaic (stand age, kNN 2020)
  gdalwarp( ## TODO: use sf gdal utils
    srcfile = file.path(paths$cache, paste0("AGB_age_mosaic_t", i, ".vrt")),
    dstfile = file.path(paths$outputs, "mosaics", paste0("AGB_age_mosaic_t", i, ".tif")),
    overwrite = TRUE
  )

  ## 3.3) Group into 5 discrete age classes
  ages_from <- c(25, 50, 80, 125, 500)
  ages_to <- c(24, 49, 79, 124, 500)

  ageRast <- classify(rast(file.path(paths$outputs, "mosaics", paste0("AGB_age_mosaic_t", i, ".tif"))),
    rcl = cbind(from = ages_from, to = ages_to, becomes = 1L:5L),
    right = FALSE, others = NA_integer_
  )

  names(ageRast) <- "ageClass"
  levels(ageRast) <- data.frame(value = 1:5, ageClass = paste0(ages_from, "-", to = ages_to))
  writeRaster(ageRast, file.path(paths$outputs, "mosaics", paste0("AGB_age_mosaic_classes_t", i, ".tif")), overwrite = TRUE)

  return(invisible(NULL))
})

parallel::stopCluster(cl)

# 4) rasterize study area by categorical zones of interest ------------------------------------
##   (basis of subsequent results comparison) WBI and ecozones by default

# targetCRS <- AGBtrends::Canada_Albers_Equal_Area_Conic
#
# source("modules/AGB_dataPrep/R/analysisZones.R")
# zoi <- createAnalysisZones(st_read("outputs/WBI_studyArea.gpkg"), targetCRS, "inputs")
# st_write(zoi, dsn="outputs/WBI_studyArea.gpkg", delete_layer = TRUE)

no_cores <- min(parallel::detectCores() / 2, 6L)
cl <- parallelly::makeClusterPSOCK(
  no_cores,
  default_packages = c("sf", "terra"),
  rscript_libs = .libPaths(),
  autoStop = TRUE
)
parallel::clusterExport(cl, varlist = c("prepZones", "no_cores"), envir = environment())
parallel::clusterEvalQ(cl, terraOptions(tempdir = paths$terra, memfrac = 0.5 / no_cores))

parallel::parLapply(cl, 1:6, function(i) {
  prepZones( # zoi = zoi,
    field = "ECOZONE", ## CAN ALTER FOR ECOREGION OR ECOPROVINCE, (which I did so results on file)
    ageClass = rast(file.path(paths$outputs, "mosaics", paste0("AGB_age_mosaic_classes_t", i, ".tif"))),
    file.id = paste0("WBI_ecozone_t", i),
    ow = TRUE
  )

  return(invisible(NULL))
}) # 17 min

parallel::stopCluster(cl)

# 5) Calculate comparative summary statistics by categorical zone of interest for all time periods -----

## Note in following that age at beginning of the 31 year time series (1984-2014) is identical to age at beginning of 't1' time interval (i.e. 1984-1988)
irast <- list(
  slope = file.path(paths$outputs, list.files(paths$outputs, pattern = "slope_mosaic")),
  w = file.path(paths$outputs, list.files(paths$outputs, pattern = "sample_size")),
  ## these last values '2' (below) refer to ageClass at 'time 0' (i.e. 1984) used for the complete time series slope raster mosaic stats assessment.
  ## this should be '1', but currently set to 6 years in b/c disturbed pixels are only known as of 1987
  ecozone = file.path(paths$outputs, list.files(paths$outputs, pattern = "ZOIxageClass_WBI_ecozone"))[-c(2, 4, 6, 8, 10, 12)][c(1:6, 2)],
  ecoregion = file.path(paths$outputs, list.files(paths$outputs, pattern = "ZOIxageClass_WBI_ecoregion"))[-c(2, 4, 6, 8, 10, 12)][c(1:6, 2)],
  ecoprovince = file.path(paths$outputs, list.files(paths$outputs, pattern = "ZOIxageClass_WBI_ecoprovince"))[-c(2, 4, 6, 8, 10, 12)][c(1:6, 2)]
)

no_cores <- length(timeint) + 1
cl <- parallelly::makeClusterPSOCK(no_cores,
  default_packages = c("dplyr", "sf", "terra"),
  rscript_libs = .libPaths(), autoStop = TRUE
)
parallel::clusterExport(cl, varlist = c("irast", "no_cores", "zoneStats"), envir = environment())
parallel::clusterEvalQ(cl, terraOptions(tempdir = paths$terra, memfrac = 0.5 / no_cores))

system.time({
  parallel::parLapply(cl, 1:7, function(i, svar = "ecozone", maskRaster = NULL) {
    ## TO DO: use maskRaster file name to qualify file.id writeRaster tag
    if (i == 7) {
      file.id <- paste0("WBI_", svar)
      # file.id <- paste0("WBI_distMask_", svar)
    } else {
      file.id <- paste0("WBI_", svar, "_t", i)
      # file.id <- paste0("WBI_distMask_", svar, "_t", i)
    }

    zoneStats(
      slopeRaster = rast(irast$slope[i]),
      weightRaster = rast(irast$w[i]),
      zoneRaster = rast(irast[[svar]][i]),
      ## maskRaster arg can be either e.g. was it disturbed? or e.g. is it forested? or both:
      ## e.g. use `rast(file.path(paths$outputs, "mosaics", "binary_disturbed_mosaic.tif"))`
      ##      for pixels disturbed over course of time series (according to ABoVE)
      maskRaster = maskRaster,
      file.id = file.id
    )

    return(invisible(NULL))
  })
}) # ~ 2 hrs

parallel::stopCluster(cl)

# 6) Diagnostic plots -------------------------------------------------------------------------

## Request 1: range, mode and mean of AGB values by ageClass, year 2000

## 1 a) range

## i) make mosaic for year 2000
tilePath <- list.files(file.path(paths$outputs, "tiles"), full.names = TRUE)
rPath <- unname(sapply(tilePath, function(x) list.files(x, pattern = "ragb", full.names = TRUE)))

agb_mosaic <- file.path(paths$outputs, "mosaics", "agb_mosaic_2000.tif")
agb_mosaic_classes <- file.path(paths$outputs, "mosaics", "agb_mosaic_2000_classes.tif")

## TODO: use sf gdal utils
sf::gdal_utils(
  util = "buildvrt",
  source = rPath,
  destination = file.path(paths$terra, "agb_2000.vrt"),
  options = c("-b", "17") ## band 17
)
sf::gdal_utils(
  util = "warp",
  source = file.path(paths$terra, "agb_2000.vrt"),
  destination = agb_mosaic
)

## i) rescale by 0.01 and classify into bins similar to Wang et al.
classify(rast(agb_mosaic) * 0.01,
  rcl = c(0, 50, 100, 150, 250),
  include.lowest = TRUE, brackets = TRUE, right = FALSE,
  filename = agb_mosaic_classes,
  overwrite = TRUE
)

## ii) compute sum of AGB (in Tg) by age class (t4 = 2000)
agbSum <- zonal(rast(agb_mosaic) * 0.09,
  rast(file.path(paths$outputs, "mosaics", "AGB_age_mosaic_classes_t4.tif")),
  fun = "sum", na.rm = TRUE
)

## iv) compute sum of AGB (in Mg) by AGB class as in Wang et al.
agbClass <- zonal(rast() * 0.09,
  rast(agb_mosaic_classes),
  fun = "sum", na.rm = TRUE
)

## iii) visualize AGB (in Tg * 0.01) by age class (Mg/ha * 0.01)
gg_agb_age_class <- ggplot(data = agbSum, aes(x = ageClass, y = I(agb_mosaic_2000 * 1e-6 * 0.01))) +
  scale_x_discrete(name = "Stand Age Class", labels = c("0-24", "25-49", "50-79", "80-124", ">= 125")) +
  scale_y_continuous(name = "AGB (Tg * 0.01)") +
  geom_bar(stat = "identity")

ggsave(file.path(paths$outputs, "figures", "AGB_distribution_x_ageClass.png"), gg_agb_age_class,
       width = 7.5, height = 4)

## iv) visualize AGB (in Tg) by AGB class as per Wang et al. (2021)
gg_agb_agb_class <- ggplot(
  data = agbClass |>
    rename(agbClass = agb_mosaic_2000, agb = agb_mosaic_2000.1) |>
    mutate(agbClass = factor(agbClass, levels = agbClass)),
  aes(x = agbClass, y = I(agb * 1e-6 * 0.01))
) +
  scale_x_discrete(
    name = "AGB Class (Mg ha-1 * 0.01)",
    labels = c("0-50", "50-100", "100-150", ">150")
  ) +
  scale_y_continuous(name = "AGB Stock (Tg * 0.01)") +
  geom_bar(stat = "identity")

ggsave(file.path(paths$outputs, "figures", "AGB_distribution_x_AGBClass.png"), gg_agb_agb_class,
       width = 7.5, height = 4)

## Request 2: cumulative delta AGB by ecozone
tilePath <- list.files(file.path(paths$outputs, "tiles"), full.names = TRUE)
agbPath <- unname(sapply(tilePath, function(x) list.files(x, pattern = "ragb", full.names = TRUE)))

no_cores <- min(parallelly::availableCores(constraints = c("connections")) / 2, 20L)
cl <- parallelly::makeClusterPSOCK(no_cores,
  default_packages = c("dplyr", "sf", "terra"),
  rscript_libs = .libPaths(),
  autoStop = TRUE
)
parallel::clusterExport(cl, varlist = c("no_cores"))
parallel::clusterEvalQ(cl, {
  terraOptions(
    tempdir = paths$terra,
    memmax = 25,
    memfrac = 0.6 / no_cores,
    progress = 1,
    verbose = TRUE
  )
})

do.call(rbind, parallel::parLapply(cl, agbPath, function(r) {
  ez <- file.path(paths$outputs, "WBI_studyArea.gpkg") |>
    st_read(quiet = TRUE) |>
    select("ECOZONE") |>
    rasterize(rast(r, lyr = "1984"), field = "ECOZONE")

  yrs <- names(rast(r))

  x <- do.call(rbind, lapply(yrs, function(yr) {
    return(zonal(rast(r, lyr = yr) * 0.09, ez, fun = "sum", na.rm = TRUE) |>
      mutate(Year = yr, .before = all_of(yr)) |>
      dplyr::rename(AGB = yr))
  }))

  return(x)
})) |>
  saveRDS(file = file.path(paths$outputs, "sum_AGB_x_Ecozone.rds"))

parallel::stopCluster(cl)

ptab <- readRDS(file.path(paths$outputs, "sum_AGB_x_Ecozone.rds")) |>
  arrange(ECOZONE, Year) |>
  group_by(ECOZONE, Year) |>
  summarize(AGB = sum(AGB, na.rm = TRUE) * 1e-6) |>
  group_by(ECOZONE) |>
  mutate(
    dAGB = c(0, diff(AGB)),
    Year = as.integer(Year)
  )

## plot

ecozones <- unique(ptab$ECOZONE)

lapply(1:length(ecozones), function(i) {
  df <- filter(ptab, ECOZONE == ecozones[i]) |> ungroup()

  # https://finchstudio.io/blog/ggplot-dual-y-axes/
  # scale and shift variables calculated based on desired mins and maxes
  max_first <- max(df$dAGB) # Specify max of first y axis
  max_second <- max(df$AGB) # Specify max of second y axis
  min_first <- min(df$dAGB) # Specify min of first y axis
  min_second <- min(df$AGB) # Specify min of second y axis

  # scale and shift variables calculated based on desired mins and maxes
  scale <- (max_second - min_second) / (max_first - min_first)
  shift <- min_first - min_second

  # Function to scale secondary axis
  scale_function <- function(x, scale, shift) {
    return((x) * scale - shift)
  }

  # Function to scale secondary variable values
  inv_scale_function <- function(x, scale, shift) {
    return((x + shift) / scale)
  }

  gg_ptab_ez <- ggplot(
    data = df |>
      bind_rows(data.frame(ECOZONE = ecozones[i], AGB = NA, Year = 1984, dAGB = 0)),
    aes(x = Year, y = dAGB, weight = c(rep(1, nrow(df)), 100), color = "dAGB (Tg)")
  ) +
    ggtitle(ecozones[i]) +
    geom_smooth(method = "loess") +
    geom_point(aes(y = inv_scale_function(c(df$AGB, NA), scale, shift), color = "Total AGB (Tg)"), pch = 20, cex = 0.75) +
    scale_x_continuous(breaks = seq.int(from = 1984, to = 2014, by = 5), labels = identity) +
    scale_y_continuous("Cumulative AGB change (Tg)",
      limits = c(min_first, max_first * 1.1),
      sec.axis = sec_axis(~ scale_function(scale, shift), name = "Total AGB (Tg)")
    ) +
    geom_hline(yintercept = 0, lty = "dashed") +
    labs(color = "Units")

  ggsave(file.path(paths$outputs, "figures", paste0("AGB_distribution_x_Year_", ecozones[i], ".png")),
         gg_ptab_ez, width = 7.5, height = 4)
})

# 7) Plot differences -------------------------------------------------------------------------

## 7 a) without disturbance mask --------------------------------------------------------------

## i=1 corresponds to 31-year time series, i=2 corresponds to time interval t1 (1984-1988), and so on and so forth
png(file = "outputs/AGB_global_trends_WBI_ecozone_x_ageClass.png", width = 7.5, height = 4, units = "in", res = 300)
plotZoneStats(file2plot = file.path(paths$outputs, "zoneStats_summary_WBI_ecozone.rds"))
dev.off()

## x Ecozone x ageClass
png(file = file.path(paths$outputs, paste0("AGB_temporal_trends_x_ECOZONE_x_ageClass_", Sys.Date(), ".png")))
lapply(plotZoneStatsIntervals(
  files2plot = file.path(paths$outputs, list.files(paths$outputs, pattern = "zoneStats_summary_WBI_ecozone_")),
  weighted = TRUE, xVar = "tp", groupVar = "ageClass", ptype = 1
), plot)
dev.off()

## x ageClass x Ecozone
## TODO: this one crashes :(
png(file = file.path(paths$outputs, paste0("AGB_temporal_trends_x_ageClass_x_ECOZONE_", Sys.Date(), ".png")))
lapply(plotZoneStatsIntervals(
  files2plot = file.path(paths$outputs, list.files(paths$outputs, pattern = "zoneStats_summary_WBI_ecozone_")),
  weighted = TRUE, xVar = "tp", catVar = "ageClass", groupVar = "ECOZONE", ptype = 2
), plot)
dev.off()

## 7 b) with disturbance mask -----------------------------------------------------------------
plotZoneStats(file2plot = file.path(paths$outputs, paste0("zoneStats_summary_WBI_distMask_ecozone.rds")))

## x ageClass x Ecozone
png(file = file.path(paths$outputs, paste0("AGB_temporal_trends_x_ECOZONE_distMask_", Sys.Date(), ".png")))
lapply(plotZoneStatsIntervals(
  files2plot = file.path(paths$outputs, list.files(paths$outputs, pattern = "zoneStats_summary_WBI_distMask_ecozone_")),
  weighted = TRUE, xVar = "tp", catVar = "ageClass", groupVar = "ECOZONE", ptype = 2
), plot)
dev.off()

gp <- plotZoneStatsIntervals(files2plot = file.path(paths$outputs, list.files(paths$outputs, pattern = "WBI_distMask_ecozone")))

# 8) Test for significant differences between groups ------------------------------------------



# Functions -----------------------------------------------------------------------------------

## TODO: move these to AGBtrends package

## Derive slope of numerical vector across a time series
slope <- function(x) {
  y <- which(!is.na(x))
  if (length(y) >= 2) {
    x <- unlist(x[!is.na(x)])
    if (length(unique(x)) == 1) {
      return(0)
    } else {
      return(sum((x - mean(x)) * (y - mean(y))) / sum((x - mean(x))^2))
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

#' Create unique categorical zones for comparative analysis
#' (e.g. study area ZOI x age class for desired time period)
prepZones <- function(zoi = st_read(file.path(paths$outputs, "WBI_studyArea.gpkg"), quiet = TRUE),
                      field = "ECOZONE",
                      ageClass = rast(file.path(paths$outputs, "mosaics", "AGB_age_mosaic_classes_t1.tif")),
                      file.id = "WBI_ecozone_t1",
                      cropRaster = NULL, # rast(file.path(paths$outputs, "tiles", "Bh11v06", "agb_slopes_Bh11v06.tif")),
                      ow = TRUE) {
  if (!is.null(cropRaster)) {
    ageClass <- crop(ageClass, cropRaster)
  }

  ## QC of field to rasterize
  # if(!inherits(class(pull(zoi, field)), "numeric")) zoi <- mutate(zoi, !!field := as.integer(factor(!!as.name(field))))

  ## define comparative zones of interest (ZOI x age class for specified time period) within study area
  rasterize(
    x = vect(zoi),
    y = ageClass,
    field = field,
    filename = file.path(paths$outputs, paste0("ZOI_", file.id, ".tif")),
    overwrite = ow
  )

  # safeguard
  maxchar <- max(nchar(levels(rast(file.path(paths$outputs, paste0("ZOI_", file.id, ".tif"))))[[1]]$value))
  rcoef <- as.numeric(str_c(c(1, rep(0, maxchar)), collapse = ""))

  ## combine study area zoi and age class mosaic into unique combined raster categories (when age class from 100 - 500 [n=5])
  zoneCat <- rast(file.path(paths$outputs, paste0("ZOI_", file.id, ".tif")))
  system.time({
    jointRast <- as.int(rcoef * ageClass + zoneCat)
  })

  ## re-assign factor levels
  ftab <- data.frame(value = unique(jointRast)) |>
    rename(value = ageClass) |>
    mutate(
      ageClass = levels(ageClass)[[1]]$ageClass[match(as.integer(str_sub(value, start = 1L, end = 1L)), levels(ageClass)[[1]]$value)],
      !!field := levels(zoneCat)[[1]][, field][match(as.integer(str_sub(value, start = -maxchar)), levels(zoneCat)[[1]]$value)]
    )

  ftab <- bind_cols(ftab, zoneCat = paste0(ftab[, 3], " (", ftab[, 2], " yrs)")) |>
    relocate(zoneCat, .after = value)

  levels(jointRast) <- ftab

  writeRaster(jointRast, filename = file.path(paths$outputs, paste0("ZOIxageClass_", file.id, ".tif")), overwrite = ow)

  return(invisible(NULL))
}

#' Calculate summary statistics by categorical zone of interest
#' (weighted mean & sd of slopes by age class, time period & study area zones of interest)
#' - WBI and ecozones by default for every time period
zoneStats <- function(slopeRaster = rast(file.path(paths$outputs, "AGB_slope_mosaic.tif")),
                      weightRaster = rast(file.path(paths$outputs, "AGB_sample_size_mosaic.tif")),
                      zoneRaster = rast(file.path(paths$outputs, "OIxageClass_WBI_ecozone_t1.tif")),
                      cropRaster = NULL, # rast(file.path(paths$outputs, "mosaics", "Bh11v06", "agb_slopes_Bh11v06.tif")),
                      maskRaster = NULL,
                      file.id = "WBI_ecozones") {
  names(weightRaster) <- "w"
  names(slopeRaster) <- "slope"

  ## crop to smaller (e.g. test) area, if one provided
  if (!is.null(cropRaster)) {
    slopeRaster <- crop(slopeRaster, cropRaster)
    weightRaster <- crop(weightRaster, cropRaster)
    zoneRaster <- crop(zoneRaster, cropRaster)
    if (!is.null(maskRaster)) maskRaster <- crop(maskRaster, cropRaster)
  }

  ## apply mask, if provided
  if (!is.null(maskRaster)) {
    slopeRaster <- mask(slopeRaster, maskRaster)
    weightRaster <- mask(weightRaster, maskRaster)
    zoneRaster <- mask(zoneRaster, maskRaster)
  } # 13 min

  a <- cats(zoneRaster)[[1]]
  levels(zoneRaster) <- NULL
  names(zoneRaster) <- "value"

  ## 1) count of slope observations by zone
  a <- a |>
    left_join(a, zonal(slopeRaster, zoneRaster, fun = "notNA"), by = "value") |>
    rename(count = slope)

  if (any(!is.na(a$count))) {
    ## 2) geometric (unweighted) mean by zone
    meanRast <- zonal(slopeRaster, zoneRaster, fun = "mean", as.raster = TRUE, na.rm = TRUE) # 12 min

    a <- a |>
      left_join(zonal(meanRast, zoneRaster, fun = "min", na.rm = TRUE), by = "value") |>
      rename(mean.slope = slope)

    ## 3) sd by zone
    a <- a |>
      left_join(
        zonal(
          x = (slopeRaster - meanRast)^2,
          z = zoneRaster, fun = "sum", na.rm = TRUE
        ), by = "value") |>
      mutate(sd = sqrt(slope / count)) |>
      select(!slope)

    ## 4.1) sum of wi*xi by zone (numerator of weighted mean)
    b <- zonal(slopeRaster * weightRaster, zoneRaster, fun = "sum", na.rm = TRUE) |>
      rename(b = slope) # ~ 11 min

    ## 4.2) sum of weights by zone (denominator of weighted mean)
    w <- zonal(weightRaster, zoneRaster, fun = "sum", na.rm = TRUE) # ~ 3.7 min

    ## 4.3) derive weighted mean
    a <- a |>
      left_join(b, by = "value") |>
      left_join(w, by = "value") |>
      mutate(wtd.mean.slope = b / w, .before = w)

    rm(w)

    ## 5) derive weighted sd
    a <- a |>
      left_join(
        zonal(
          x = weightRaster * (slopeRaster - classify(zoneRaster, rcl = a[, c("value", "wtd.mean.slope")]))^2,
          z = zoneRaster, fun = "sum", na.rm = TRUE
        ) |>
          rename(x = w),
        by = "value"
      ) |>
      mutate(wtd.sd = sqrt(x / w)) |>
      filter(!is.na(count)) |>
      select(!c(b, w, x))

    ##  6 b) write to file
    saveRDS(a, file = paste0("outputs/zoneStats_summary_", file.id, ".rds"))
  }

  return(invisible(NULL))
}

#' @importFrom dplyr filter mutate select
#' @importFrom ggplot2 aes geom_errorbar geom_hline geom_line ggplot scale_x_discrete xlab ylab
#' @importFrom stringr str_c str_detect str_sub
plotZoneStats <- function(file2plot = "outputs/zoneStats_summary_WBI_ecozone.rds",
                          weighted = TRUE, out = 1) {
  sumtab <- readRDS(file2plot) |>
    mutate(ageClass = factor(ageClass, levels = unique(ageClass), ordered = TRUE))
  catvar <- names(sumtab)[4]

  # y-axis label corresponding to examined time period
  tref <- cbind.data.frame(
    tp = paste0("t", 1:6),
    yr1 = c(1984, 1989, 1994, 1999, 2004, 2009),
    yr2 = c(1988, 1993, 1998, 2003, 2008, 2014)
  )

  tp <- ifelse(str_detect(file2plot, "_t"),
    paste0("(", str_c(filter(tref, tp == str_sub(file2plot, start = -6L, end = -5L)) |>
                        select(yr1, yr2) |> unlist(), collapse = "-"), ")"),
    "(1984-2014)"
  )

  if (weighted) {
    meanVar <- "wtd.mean.slope"
    sdVar <- "wtd.sd"
  } else {
    meanVar <- "mean.slope"
    sdVar <- "sd"
  }

  gp <- ggplot(sumtab, aes(x = ageClass, y = !!as.name(meanVar), group = !!as.name(catvar), color = !!as.name(catvar))) +
    xlab("stand age") +
    ylab(paste0("n-weighted mean AGB trend ", tp)) +
    scale_x_discrete(breaks = sumtab$ageClass, labels = sumtab$ageClass) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_line(linetype = 3) +
    geom_errorbar(aes(ymin = !!as.name(meanVar) - (!!as.name(sdVar) / sqrt(count)), ymax = !!as.name(meanVar) + (!!as.name(sdVar) / sqrt(count))),
      width = .2
    )

  if (out == 1) plot(gp) else return(gp)
}

#' @importFrom dplyr mutate relocate select
#' @importFrom ggplot2 aes facet_grid geom_hline geom_line ggplot ggtitle labs xlab ylab
#' @importFrom stringr str_sub
plotZoneStatsIntervals <- function(files2plot = file.path(paths$outputs, list.files(paths$outputs, pattern = "zoneStats_summary_WBI_ecozone_")),
                                   weighted = TRUE, xVar = "ageClass", groupVar = "tp", catVar = NULL, ptype = 1, plotResult = TRUE) {
  ## compile summary tables
  sumtab <- do.call(rbind, lapply(files2plot, function(x) {
    y <- readRDS(x) |>
      mutate(ageClass = factor(ageClass, levels = unique(ageClass), ordered = TRUE))
    y$tp <- factor(rep(str_sub(x, start = -6L, end = -5L), nrow(y)),
      levels = paste0("t", 1:6), ordered = TRUE
    )
    return(y |> relocate(tp, .before = value))
  }))

  ## validate arguments
  if (weighted) {
    meanVar <- "wtd.mean.slope"
    sdVar <- "wtd.sd"
  } else {
    meanVar <- "mean.slope"
    sdVar <- "sd"
  }

  if (is.null(catVar)) catVar <- names(sumtab)[5]
  if (!is.null(groupVar)) groupVar <- match.arg(groupVar, names(sumtab))
  xVar <- match.arg(xVar, names(sumtab))
  if (xVar == "ageClass") xlabel <- "stand age"
  if (xVar == "tp") xlabel <- "5-year time period (1984-2014)"

  if (ptype == 1) {
    gp <- ggplot(sumtab, aes(x = !!as.name(xVar), y = !!as.name(meanVar), group = !!as.name(groupVar), color = !!as.name(groupVar))) +
      labs(
        x = xlabel,
        y = "n-weighted mean AGB trend"
      ) +
      geom_hline(yintercept = 0, linetype = "dotted") +
      geom_line() +
      facet_grid(~ .data[[catVar]])
  } else {
    if (ptype == 2) {
      gp <- lapply(unique(sumtab[, catVar]), function(cat) {
        x <- filter(sumtab, !!as.name(catVar) == cat)

        return(ggplot(x, aes(x = !!as.name(xVar), y = !!as.name(meanVar), group = !!as.name(groupVar), color = !!as.name(groupVar))) +
          geom_line() +
          xlab(xlabel) +
          ylab("n-weighted mean AGB trend") +
          geom_hline(yintercept = 0, linetype = "dotted") +
          geom_line() +
          ggtitle(cat))
      })
      names(gp) <- unique(sumtab[, catVar])

      # pdf("all.pdf")
      # invisible(lapply(gp, print))
      # dev.off()
    }

    if (!plotResult) {
      return(gp)
    }
    if (ptype == 1) print(gp) else lapply(gp, print)
  }
}

## Test code using extracted cell values. Do results agree? Tested and yes :)
# rstack <- as.data.frame(c(slopeRaster, weightRaster, zoneRaster))
#
# group_by(rstack, zone) |>
#   filter(!is.na(slope) & !is.na(zone)) |>
#   dplyr::summarize(mean = format(mean(slope), scientific = TRUE),
#                    sd = format(sd(slope), scientific = TRUE),
#                    wmean = format(sum(w*slope) / sum(w), scientific = TRUE),
#                    meanXw = sum(w * slope),
#                    wmean2 = format(weighted.mean(slope, w, na.rm = TRUE), scientific = TRUE),
#                    wSD1 = format(weighted.sd(slope, w, na.rm = TRUE), scientific = TRUE),
#                    wSD2 = format(sdwt(slope, w), scientific = TRUE))

# cleanup -------------------------------------------------------------------------------------
unlink(paths$terra, recursive = TRUE)
