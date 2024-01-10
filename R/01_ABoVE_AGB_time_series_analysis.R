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
Require::Install(c("cowplot", "gridGraphics"), upgrade = FALSE)
Require::Require(c("dplyr", "ggplot2", "reproducible", "sf", "stringr", "terra",
                   "PredictiveEcology/AGBtrends (>= 0.0.4)"), upgrade = FALSE)

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
paths$tiles <- file.path(paths$outputs, "tiles") |>
  fs::dir_ls(regexp = "Bh", type = "directory") |>
  sort()

## set the max number of cores to use for parallel computations
no_cores <- AGBtrends::getNumCores()

file.remove(list.files(paths$terra, full.names = TRUE)) ## preemptive cleanup

terraOptions(
  tempdir = paths$terra,
  memmax = 25,
  memfrac = 0.8,
  progress = 1,
  verbose = TRUE,
  todisk = TRUE
)

## define time intervals (year ranges between 1984-2014)
timeint <- list(t1 = 1:5, t2 = 6:10, t3 = 11:15, t4 = 16:20, t5 = 21:25, t6 = 26:31)
timeint_all <- timeint |> unlist() |> unname() |> list(all = _)

## 1.1) Estimate cell-wise linear regression coefficients for undisrupted time series ---------
## aka "local" or "geographically weighted regression (GWR)"

f1 <- AGBtrends::gwr(paths$tiles, type = "slopes", cores = length(paths$tiles))
f2 <- AGBtrends::gwr(paths$tiles, type = "sample_size", cores = length(paths$tiles))

## 1.2) Combine tiled slope rasters into unified mosaics --------------------------------------

f3a <- AGBtrends::buildMosaics("slopes", intervals = timeint_all, src = paths$tiles, dst = paths$outputs)
f3b <- AGBtrends::buildMosaics("sample_size", intervals = timeint_all, src = paths$tiles, dst = paths$outputs)
f3 <- c(f3a, f3b)

## TODO: rebuild sample_size input tiles:
##  Warning messages:
##    1: In CPL_gdalbuildvrt(if (missing(source)) character(0) else source,  :
##       GDAL Message 1: gdalbuildvrt does not support heterogeneous band data type: expected Byte, got Float32.
##       Skipping /mnt/projects/CBM/2BT/ForProd/outputs/studyArea_WBI/tiles/Bh06v08/agb_sample_size_Bh06v08.tif


## 2.1) Calculate cell-specific slopes per 5-year time interval (n=6) -------------------------

### 2.1.1) calculate local slope coefficient for specified time interval ----------------------
f4 <- AGBtrends::gwrt(paths$tiles, type = "slopes", cores = no_cores, intervals = timeint)

### 2.1.2) stock number of non-NA values for subsequent weighted standard deviation -----------
f5 <- AGBtrends::gwrt(paths$tiles, type = "sample_size", cores = no_cores, intervals = timeint)

## 2.2) Combine tiled slope rasters into numerous unified mosaics -----------------------------

f6a <- AGBtrends::buildMosaics("slopes", intervals = timeint, src = paths$tiles, dst = paths$outputs)
f6b <- AGBtrends::buildMosaics("sample_size", intervals = timeint, src = paths$tiles, dst = paths$outputs)
f6 <- c(f6a, f6b)

## TODO: rebuild sample_size input tiles:
##  Warning messages:
##    1: In CPL_gdalbuildvrt(if (missing(source)) character(0) else source,  :
##       GDAL Message 1: gdalbuildvrt does not support heterogeneous band data type: expected Byte, got Float32.
##       Skipping /mnt/projects/CBM/2BT/ForProd/outputs/studyArea_WBI/tiles/Bh06v08/agb_sample_size_Bh06v08.tif

## Visual examination of results --------------------------------------------------------------

## verify hashes and file sizes
sapply(names(timeint), function(tp) {
  digest::digest(file.path(paths$outputs, paste0("AGB_slope_mosaic_", tp, ".tif")), algo = "xxhash64")
})
sapply(names(timeint), function(tp) {
  file.size(file.path(paths$outputs, paste0("AGB_slope_mosaic_", tp, ".tif")))
})

## TODO: finesse these plots further

plot_slope_mosaics <- function() {
  par(mfrow = c(3, 2))
  for (tp in names(timeint)) {
    plot(rast(file.path(paths$outputs, paste0("AGB_slope_mosaic_", tp, ".tif"))), main = tp)
  }
}

gg_slope_mosaics <- cowplot::plot_grid(plot_slope_mosaics)

ggsave(file.path(paths$outputs, "figures", "gg_slope_mosaics.png"),
       gg_slope_mosaics, height = 5, width = 10)

plot_slope_mosaic_hists <- function() {
  par(mfrow = c(2, 3))
  hist(rast(file.path(paths$outputs, paste0("AGB_slope_mosaic_t1.tif"))), main = "t1")
  hist(rast(file.path(paths$outputs, paste0("AGB_slope_mosaic_t2.tif"))), main = "t2")
  hist(rast(file.path(paths$outputs, paste0("AGB_slope_mosaic_t3.tif"))), main = "t3")
  hist(rast(file.path(paths$outputs, paste0("AGB_slope_mosaic_t4.tif"))), main = "t4")
  hist(rast(file.path(paths$outputs, paste0("AGB_slope_mosaic_t5.tif"))), main = "t5")
  hist(rast(file.path(paths$outputs, paste0("AGB_slope_mosaic_t6.tif"))), main = "t6")
}

gg_slope_mosaic_hists <- cowplot::plot_grid(plot_slope_mosaic_hists)

ggsave(file.path(paths$outputs, "figures", "gg_slope_mosaic_hists.png"),
       gg_slope_mosaic_hists, height = 5, width = 10)

# 3) Group slopes by age at time x ------------------------------------------------------------
##    (band argument determines reference layer/year),
##    effectively masking out pixels disturbed mid-time series

f7 <- AGBtrends::buildMosaics("age", intervals = timeint, src = paths$tiles, dst = paths$outputs)

# 4) Evaluate frequency distributions of forest landcover -----------------------------------

ftab <- readRDS(file.path(paths$outputs, "ABoVE_LandCover_freq_tables.rds"))

ftab <- do.call(rbind, lapply(1:length(ftab), function(i) {
  return(bind_cols(
    ecozone = cats(rast(ecoRast_tif))[[1]]$ECOZONE[i],
    group_by(ftab[[i]], value) |>
      summarize(count = mean(count, na.rm = TRUE)) |>
      left_join(reftab, by = "value") |>
      relocate(cat, .before = count)
  ))
}))

mutate(ftab, cat2 = ifelse(value %in% c(1:4), "Forested", "Non-Forested")) |>
  group_by(ecozone, cat2) |>
  summarize(catCount = sum(count)) |>
  group_by(ecozone) |>
  summarize(pcntForestedPixels = catCount[1] / sum(catCount))

# 5) rasterize study area by categorical zones of interest ------------------------------------
##   (basis of subsequent results comparison) WBI and ecozones by default

# targetCRS <- AGBtrends::Canada_Albers_Equal_Area_Conic
#
# source("modules/AGB_dataPrep/R/analysisZones.R")
# f_sA <- file.path(paths$outputs, "WBI_studyArea.gpkg")
# zoi <- createAnalysisZones(st_read(f_sA), targetCRS, "inputs")
# st_write(zoi, dsn = f_sA, delete_layer = TRUE)

no_cores <- AGBtrends::getNumCores(length(timeint))
cl <- parallelly::makeClusterPSOCK(
  no_cores,
  default_packages = c("AGBtrends", "sf", "terra"),
  rscript_libs = .libPaths(),
  autoStop = TRUE
)
parallel::clusterExport(cl, varlist = c("no_cores"), envir = environment())
parallel::clusterEvalQ(cl, terraOptions(tempdir = paths$terra, memfrac = 0.5 / no_cores))

parallel::parLapply(cl, seq(length(timeint)), function(i) {
  prepZones( # zoi = zoi,
    field = "ECOZONE", ## CAN ALTER FOR ECOREGION OR ECOPROVINCE, (which I did so results on file)
    ageClass = rast(file.path(paths$outputs, "mosaics", paste0("AGB_age_mosaic_classes_t", i, ".tif"))),
    file.id = paste0("WBI_ecozone_t", i),
    ow = TRUE
  )

  return(invisible(NULL))
}) # 17 min

parallel::stopCluster(cl)

# 6) Calculate comparative summary statistics by categorical ZOI for all time periods ---------

## Note in following that age at beginning of the 31 year time series (1984-2014)
## is identical to age at beginning of 't1' time interval (i.e. 1984-1988)
files <- list(
  list.files(paths$outputs, pattern = "slope_mosaic", full.names = TRUE),
  list.files(paths$outputs, pattern = "sample_size", full.names = TRUE),
  list.files(paths$outputs, pattern = "ZOIxageClass_WBI_ecozone", full.names = TRUE),
  list.files(paths$outputs, pattern = "ZOIxageClass_WBI_ecoregion", full.names = TRUE),
  list.files(paths$outputs, pattern = "ZOIxageClass_WBI_ecoprovince", full.names = TRUE)
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

no_cores <- AGBtrends::getNumCores(length(timeint) + 1)
cl <- parallelly::makeClusterPSOCK(no_cores,
  default_packages = c("AGBtrends", "dplyr", "sf", "terra"),
  rscript_libs = .libPaths(), autoStop = TRUE
)
parallel::clusterExport(cl, varlist = c("irast", "no_cores"), envir = environment())
parallel::clusterEvalQ(cl, terraOptions(tempdir = paths$terra, memfrac = 0.5 / no_cores))

system.time({
  parallel::parLapply(cl, seq(no_cores), function(i, svar = "ecozone", maskRaster = NULL) {
    ## TODO: use maskRaster file name to qualify file.id writeRaster tag
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

# 7) Diagnostic plots -------------------------------------------------------------------------

## TODO: move these to AGBtrends package

## Request 1: range, mode and mean of AGB values by ageClass, year 2000 -----------------------

## 1 a) range

## i) make mosaic for year 2000

agb_tifs <- fs::dir_ls(paths$tiles, regexp = "ragb", recurse = TRUE) ## 81 files

agb_mosaic <- file.path(paths$outputs, "mosaics", "agb_mosaic_2000.tif")
agb_mosaic_classes <- file.path(paths$outputs, "mosaics", "agb_mosaic_2000_classes.tif")

sf::gdal_utils(
  util = "buildvrt",
  source = agb_tifs,
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
       width = 8, height = 4)

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
       width = 8, height = 4)

## Request 2: cumulative delta AGB by ecozone -------------------------------------------------

agb_tifs <- fs::dir_ls(paths$tiles, regexp = "ragb", recurse = TRUE) ## 81 files

no_cores <- AGBtrends::getNumCores(20L) ## TODO: why 20 here; RAM limitations?
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

do.call(rbind, parallel::parLapply(cl, agb_tifs, function(r) {
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
    return(x * scale - shift)
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
         gg_ptab_ez, width = 8, height = 4)
})

# 8) Plot differences -------------------------------------------------------------------------

## TODO: use ggsave()

## 8 a) without disturbance mask --------------------------------------------------------------

## i=1 corresponds to 31-year time series, i=2 corresponds to time interval t1 (1984-1988), and so on and so forth
gg_71 <- plotZoneStats(
  file2plot = file.path(paths$outputs, "summaries", "zoneStats_summary_WBI_ecozone.rds")
)

ggsave(
  file.path(paths$outputs, "figures", "AGB_global_trends_WBI_ecozone_x_ageClass.png"),
  gg_71,
  width = 8,
  height = 4
)

## x Ecozone x ageClass
files2plot <- file.path(paths$outputs, "summaries") |>
  list.files(pattern = "zoneStats_summary_WBI_ecozone_", full.names = TRUE)

gg_72 <- plotZoneStatsIntervals(files2plot, weighted = TRUE, xVar = "tp", groupVar = "ageClass", ptype = 1)

ggsave(
  file.path(paths$outputs, "figures", paste0("AGB_temporal_trends_x_ECOZONE_x_ageClass_", Sys.Date(), ".png")),
  gg_72,
  width = 8,
  height = 4
)

## x ageClass x Ecozone
files2plot <- file.path(paths$outputs, "summaries") |>
  list.files(pattern = "zoneStats_summary_WBI_ecozone_", full.names = TRUE)

gg_73 <- plotZoneStatsIntervals(files2plot, weighted = TRUE, xVar = "tp",
                                catVar = "ageClass", groupVar = "ECOZONE",
                                ptype = 2, plotResult = FALSE) |>
  cowplot::plot_grid(plotlist = _)

ggsave(
  file.path(paths$outputs, "figures", paste0("AGB_temporal_trends_x_ageClass_x_ECOZONE_", Sys.Date(), ".png")),
  gg_73,
  width = 10,
  height = 5
)

## 8 b) with disturbance mask -----------------------------------------------------------------
gg_74 <- plotZoneStats(
  file2plot = file.path(paths$outputs, "summaries", "zoneStats_summary_WBI_distMask_ecozone.rds")
)

## TODO: verify & adjust output filename
# ggsave(
#   file.path(paths$outputs, "figures", paste0("AGB_temporal_trends_x_ECOZONE_distMask_", Sys.Date(), ".png")),
#   gg_74
# )

## x ageClass x Ecozone
files2plot <- file.path(paths$outputs, "summaries") |>
  list.files(pattern = "zoneStats_summary_WBI_distMask_ecozone_", full.names = TRUE)

gg_75 <- plotZoneStatsIntervals(files2plot, weighted = TRUE, xVar = "tp",
                                catVar = "ageClass", groupVar = "ECOZONE",
                                ptype = 2, plotResult = FALSE) |>
  cowplot::plot_grid(plotlist = _) ## TODO: why is `0-24` the ontly ageClass???

ggsave(
  file.path(paths$outputs, "figures", paste0("AGB_temporal_trends_x_ECOZONE_distMask_", Sys.Date(), ".png")),
  gg_75
)

files2plot = file.path(paths$outputs, "summaries") |>
  list.files(pattern = "WBI_distMask_ecozone", full.names = TRUE)

gg_76 <- plotZoneStatsIntervals(files2plot) ## TODO: why is `0-24` the ontly ageClass???

## TODO: verify & adjust output filename
# ggsave(
#   file.path(paths$outputs, "figures", paste0("AGB_temporal_trends_x_ECOZONE_distMask_", Sys.Date(), ".png")),
#   gg_76
# )

# 9) Test for significant differences between groups ------------------------------------------

## TODO

# cleanup -------------------------------------------------------------------------------------
unlink(paths$terra, recursive = TRUE)
