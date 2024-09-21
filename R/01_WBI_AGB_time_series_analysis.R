##  Updated: Sept 20, 2024
##
##  Authors: Alex M. Chubaty, PhD, FOR-CAST Research & Analytics
##           Tyler Rudolph, biologist M.Sc., CFS/NRCAN, Trade, Economics & Industry Branch
##
##  Description :
##   Recreate the ABoVE AGB time series analyses using simulated biomass, age, and disturbance
##   maps from the Western Boreal Initiative, which forecast vegetation and wildfire dynamics
##   from 2011-2100 under multiple climate change scenarios.
##
##   WBI simulations used 250m pixels, and output above ground biomass values in g/m^2;
##   ABoVE data product used 30m pixels, with AGB in Mg/ha with scaling factor of 0.01
##   (ABoVE values should be multiplied by the scale factor to obtain true values).
##
##   ABoVE: Mg/ha = 10e6 g / 10e4 m^2 = 10e2 g / m^2; scaled by 0.01 gives 10^4 g/m^2
##   therefore, WBI values need to be divided by 100 to be comparable to ABoVE.

# packages ------------------------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(googledrive)
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
allStudyAreas <- c("BC", "AB", "NT", "SK", "YT") ## MB doesn't overlap with ABoVE

## set the max number of cores to use for parallel computations
# options(parallelly.availableCores.fallback = 4L) ## set to limit the number of cores

pixelAreaHa <- 6.25 ## (250m * 250m) / (10e4 m^2/h) = 6.25 ha per pixel

## define time intervals (year ranges between 2011-2100)
timeint <- list(
  t1 = 1:5, t2 = 6:10, t3 = 11:15, t4 = 16:20, t5 = 21:25, t6 = 26:30,
  t7 = 31:35, t8 = 36:40, t9 = 41:45, t10 = 46:50, t11 = 51:55, t12 = 56:60,
  t13 = 61:65, t14 = 66:70, t15 = 71:75, t16 = 76:80, t17 = 81:85, t18 = 86:90
)
timeint_all <- timeint |> unlist() |> unname() |> list(all = _)

years <- (timeint_all |> unlist() |> unname()) + 2010
n_int <- length(timeint)

for (thisRep in allReps) {
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
  paths$figures <- file.path(paths$outputs, "figures") |> checkPath(create = TRUE)
  paths$mosaics <- file.path(paths$outputs, "mosaics") |> checkPath(create = TRUE)
  paths$summaries <- file.path(paths$outputs, "summaries") |> checkPath(create = TRUE)
  paths$terra <- file.path(paths$scratch, "terra", climateScenario, "01") |> checkPath(create = TRUE)

  file.remove(list.files(paths$terra, full.names = TRUE)) ## preemptive cleanup

  terraOptions(tempdir = paths$terra)

  ## 1.1) Estimate cell-wise linear regression coefficients for undisrupted time series ---------
  ## aka "local" or "geographically weighted regression (GWR)"

  agbdirs <- file.path(paths$outputs, "agb") |>
    list.dirs(recursive = FALSE)

  ## set the max number of cores to use for parallel computations
  no_cores <- min(
    as.integer(terra::free_RAM() / 1024^2 / 65), ## ~65 GB per thread
    AGBtrends::getNumCores(length(agbdirs))
  )

  cl <- parallelly::makeClusterPSOCK(
    no_cores,
    default_packages = c("AGBtrends", "terra"),
    rscript_libs = .libPaths(),
    autoStop = TRUE
  )

  ## TODO: more efficient to parallelize across provs rather than using app(); rework pkg funs
  f1 <- parallel::parLapply(cl, agbdirs, function(d) {
    AGBtrends::gwr(d, type = "slopes", cores = 1)
  }) |>
    unlist()

  f2 <- parallel::parLapply(cl, agbdirs, function(d) {
    AGBtrends::gwr(d, type = "sample_size", cores = 1)
  }) |>
    unlist()

  parallel::stopCluster(cl)
  terra::tmpFiles(remove = TRUE)
  gc()

  ## 1.2) Combine tiled slope rasters into unified mosaics --------------------------------------

  source("R/WBI_buildMosaics.R")

  f3a <- WBI_buildMosaics("slopes", intervals = timeint_all, src = agbdirs, dst = paths$mosaics)
  f3b <- WBI_buildMosaics("sample_size", intervals = timeint_all, src = agbdirs, dst = paths$mosaics)
  f3 <- c(f3a, f3b)

  ## 2.1) Calculate cell-specific slopes per 5-year time interval (n=18) -------------------------

  no_cores <- min(
    as.integer(terra::free_RAM() / 1024^2 / 20), ## <20 GB per thread
    AGBtrends::getNumCores(2 * length(agbdirs))
  )

  cl <- parallelly::makeClusterPSOCK(
    no_cores,
    default_packages = c("AGBtrends", "terra"),
    rscript_libs = .libPaths(),
    autoStop = TRUE
  )

  parallel::clusterExport(cl, c("agbdirs", "timeint"))

  ## TODO: currently more efficient to parallelize across provs rather than using app(); rework pkg funs
  ## and it's better still to do both sets of calcs in parallel
  f45 <- parallel::parLapply(cl, seq(2 * length(agbdirs)), function(i) {
    if (i <= length(agbdirs)) {
      ### 2.1.1) calculate local slope coefficient for specified time interval
      d <- agbdirs[i]
      AGBtrends::gwrt(d, type = "slopes", cores = 1, intervals = timeint)
    } else {
      ### 2.1.2) stock number of non-NA values for subsequent weighted standard deviation
      d <- agbdirs[i - length(agbdirs)]
      AGBtrends::gwrt(d, type = "sample_size", cores = 1, intervals = timeint)
    }
  }) |>
    unlist() ## ~7 mins per interval for NT; 18 intervals; ~2.1h total time

  f4 <- grep("slopes", f45, value = TRUE)
  f5 <- grep("sample_size", f45, value = TRUE)

  parallel::stopCluster(cl)

  ## 2.2) Combine tiled slope rasters into numerous unified mosaics -----------------------------

  source("R/WBI_buildMosaics.R")

  f6a <- WBI_buildMosaics("slopes", intervals = timeint, src = agbdirs, dst = paths$mosaics)
  f6b <- WBI_buildMosaics("sample_size", intervals = timeint, src = agbdirs, dst = paths$mosaics)
  f6 <- c(f6a, f6b)

  ## Visual examination of results --------------------------------------------------------------

  ## verify hashes and file sizes
  sapply(names(timeint), function(tp) {
    file.path(paths$mosaics, paste0("agb_slopes_mosaic_", tp, ".tif")) |>
      digest::digest(algo = "xxhash64")
  })
  sapply(names(timeint), function(tp) {
    file.path(paths$mosaics, paste0("agb_slopes_mosaic_", tp, ".tif")) |>
      file.size()
  })

  ## TODO: finesse these plots further

  plot_slope_mosaics <- function() {
    par(mfrow = c(3, 6))
    for (tp in names(timeint)) {
      file.path(paths$mosaics, paste0("agb_slopes_mosaic_", tp, ".tif")) |>
        rast() |>
        plot(main = tp)
    }
  }

  gg_slope_mosaics <- cowplot::plot_grid(plot_slope_mosaics)

  ggsave(file.path(paths$figures, "gg_slopes_mosaics.png"),
         gg_slope_mosaics, height = 12, width = 16)

  plot_slope_mosaic_hists <- function() {
    par(mfrow = c(3, 6))
    for (i in seq_len(n_int)) {
      file.path(paths$mosaics, paste0("agb_slopes_mosaic_t", i, ".tif")) |>
        rast() |>
        hist(main = paste0("t", i), maxcell = 1e+08)
    }
  }

  gg_slope_mosaic_hists <- cowplot::plot_grid(plot_slope_mosaic_hists)

  ggsave(file.path(paths$figures, "gg_slopes_mosaic_hists.png"),
         gg_slope_mosaic_hists, height = 12, width = 16)

  # 3) Group slopes by age at time x ------------------------------------------------------------
  ##    (band argument determines reference layer/year),
  ##    effectively masking out pixels disturbed mid-time series

  agedirs <- file.path(paths$outputs, "age") |>
    list.dirs(recursive = FALSE)

  source("R/WBI_buildMosaics.R")

  f7 <- WBI_buildMosaics(type = "age", intervals = timeint, src = agedirs, dst = paths$mosaics) ## ~570 GB

  gc()

  # 4) Evaluate frequency distributions of forest landcover -----------------------------------

  if (FALSE) { ## TODO
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
  }

  # 5) rasterize study area by categorical zones of interest ------------------------------------
  ##   (basis of subsequent results comparison) WBI and ecozones by default

  no_cores <- AGBtrends::getNumCores(n_int)
  cl <- parallelly::makeClusterPSOCK(
    no_cores,
    default_packages = c("AGBtrends", "sf", "terra"),
    rscript_libs = .libPaths(),
    autoStop = TRUE
  )
  parallel::clusterExport(cl, varlist = c("no_cores", "paths"), envir = environment())
  parallel::clusterEvalQ(cl, {
    terraOptions(tempdir = paths$terra, memfrac = 0.5 / no_cores)

    invisible(NULL)
  })

  for (eco in c("ECOZONE", "ECOREGION", "ECOPROVINCE")) {
    parallel::parLapply(cl, seq(n_int), function(i, eco) {
      targetCRS <- file.path(paths$outputs, "agb", "AB", "ragb_AB.tif") |>
        rast() |>
        crs()

      zoi <- file.path("outputs", "studyArea_WBI", "WBI_studyArea.gpkg") |>
        st_read(quiet = TRUE) |>
        st_transform(targetCRS)

      prepZones(
        zoi = zoi,
        field = eco,
        ageClass = rast(file.path(paths$mosaics, paste0("agb_age_mosaic_classes_t", i, ".tif"))),
        fileID = paste0("WBI_", tolower(eco), "_t", i),
        destinationPath = paths$outputs,
        overwrite = TRUE
      )

      return(invisible(NULL))
    },  eco = eco) # ~110 GB
  }

  parallel::stopCluster(cl)

  # 6) Calculate comparative summary statistics by categorical ZOI for all time periods ---------

  ## Note in following that age at beginning of the 90 year time series (2011-2100)
  ## is identical to age at beginning of 't1' time interval (i.e. 2011-2015)
  files <- list(
    list.files(file.path(paths$mosaics), pattern = "slopes_mosaic", full.names = TRUE),
    list.files(file.path(paths$mosaics), pattern = "sample_size", full.names = TRUE),
    list.files(paths$outputs, pattern = "ZOIxageClass_WBI_ecozone", full.names = TRUE),
    list.files(paths$outputs, pattern = "ZOIxageClass_WBI_ecoregion", full.names = TRUE),
    list.files(paths$outputs, pattern = "ZOIxageClass_WBI_ecoprovince", full.names = TRUE)
  ) |>
    lapply(grep, pattern = "[.]tif$", value = TRUE) ## only tifs; not aux files

  irast <- list(
    slope = files[[1]],
    w = files[[2]],
    ## these last values '1' (below) refer to ageClass at 'time 0' (i.e. 2011) used
    ## for the complete time series slope raster mosaic stats assessment.
    ecozone = files[[3]][c(seq(n_int), 1)],
    ecoregion = files[[4]][c(seq(n_int), 1)],
    ecoprovince = files[[5]][c(seq(n_int), 1)]
  )

  no_cores <- AGBtrends::getNumCores(n_int + 1) ## ~36 GB per thread

  cl <- parallelly::makeClusterPSOCK(no_cores,
                                     default_packages = c("AGBtrends", "dplyr", "sf", "terra"),
                                     rscript_libs = .libPaths(), autoStop = TRUE
  )

  parallel::clusterExport(cl, varlist = c("irast", "no_cores", "paths"), envir = environment())

  parallel::clusterEvalQ(cl, {
    terraOptions(tempdir = paths$terra, memfrac = 0.5 / no_cores, todisk = TRUE)

    invisible(NULL)
  })

  parallel::parLapply(cl, seq(no_cores), function(i, svar = "ecozone") {
    ## maskRaster arg can be either e.g. was it disturbed? or e.g. is it forested? or both?
    ## e.g., for pixels disturbed over course of time series (according to ABoVE) use:
    ##   rast(file.path(paths$outputs, "mosaics", "binary_disturbed_mosaic.tif")) |>
    ##     classify(cbind(1, NA))

    ## run without mask
    if (i == no_cores) {
      file.id <- paste0("WBI_", svar)
    } else {
      file.id <- paste0("WBI_", svar, "_t", i)
    }

    zoneStats(
      slopeRaster = rast(irast$slope[i]),
      weightRaster = rast(irast$w[i]),
      zoneRaster = rast(irast[[svar]][i]),
      maskRaster = NULL,
      fileID = file.id,
      destinationPath = paths$summaries
    )

    gc()

    ## run with disturbance mask
    if (i == no_cores) {
      file.id <- paste0("WBI_distMask_", svar)
    } else {
      file.id <- paste0("WBI_distMask_", svar, "_t", i)
    }

    msk <- file.path(paths$mosaics, "binary_disturbed_mosaic.tif") |>
      rast()
    stopifnot(nlyr(msk) == 1)
    msk <- classify(msk, cbind(1, NA))

    zoneStats(
      slopeRaster = rast(irast$slope[i]),
      weightRaster = rast(irast$w[i]),
      zoneRaster = rast(irast[[svar]][i]),
      maskRaster = msk,
      fileID = file.id,
      destinationPath = paths$summaries
    )

    gc()

    return(invisible(NULL))
  })

  parallel::stopCluster(cl)
  terra::tmpFiles(remove = TRUE)
  gc()

  # 7) Diagnostic plots -------------------------------------------------------------------------

  ## TODO: move these to AGBtrends package

  ## Request 1: range, mode and mean of AGB values by ageClass ----------------------------------

  ## 1 a) range

  ## i) make mosaic for year 2100

  myear <- 2100 ## TODO: which other years should we look at?
  agb_tifs <- fs::dir_ls(agbdirs, regexp = "ragb", recurse = TRUE)
  agb_vrt <- file.path(paths$terra, paste0("agb_", myear, ".vrt"))

  agb_mosaic <- file.path(paths$mosaics, paste0("agb_mosaic_", myear, ".tif"))
  agb_mosaic_classes <- file.path(paths$mosaics, paste0("agb_mosaic_", myear, "_classes.tif"))

  sf::gdal_utils(
    util = "buildvrt",
    source = agb_tifs,
    destination = agb_vrt,
    options = c("-b", as.character(length(years)))
  )

  sf::gdal_utils(
    util = "warp",
    source = agb_vrt,
    destination = agb_mosaic
  )

  ## i) classify into Mg/ha bins similar to Wang et al.
  classify(rast(agb_mosaic) * 0.01, ## convert to Mg/ha
           rcl = c(0, 50, 100, 150, 250), ## Mg/ha bins
           include.lowest = TRUE, brackets = TRUE, right = FALSE,
           filename = agb_mosaic_classes,
           overwrite = TRUE
  )

  ## ii) compute sum of AGB (in Mg) by age class (t18 = 2100)
  agbSum <- zonal(rast(agb_mosaic) * pixelAreaHa,
                  rast(file.path(paths$mosaics, "agb_age_mosaic_classes_t18.tif")),
                  fun = "sum", na.rm = TRUE
  )

  ## iv) compute sum of AGB (in Mg) by AGB class as in Wang et al.
  agbClass <- zonal(rast(agb_mosaic) * pixelAreaHa,
                    rast(agb_mosaic_classes),
                    fun = "sum", na.rm = TRUE
  )

  ## iii) visualize AGB (in Tg * 0.01) by age class
  gg_agb_age_class <- ggplot(data = agbSum, aes(x = ageClass, y = I(agb_mosaic_2100 * 1e-6 * 0.01))) +
    scale_x_discrete(name = "Stand Age Class", labels = c("0-24", "25-49", "50-79", "80-124", ">= 125")) +
    scale_y_continuous(name = "AGB (Tg)") +
    geom_bar(stat = "identity")

  ggsave(file.path(paths$figures, "AGB_distribution_x_ageClass.png"), gg_agb_age_class,
         width = 16, height = 12)

  ## iv) visualize AGB (in Tg * 0.01) by AGB class as per Wang et al. (2021)
  gg_agb_agb_class <- ggplot(
    data = agbClass |>
      rename(agbClass = agb_mosaic_2100, agb = agb_mosaic_2100.1) |>
      mutate(agbClass = factor(agbClass, levels = agbClass)),
    aes(x = agbClass, y = I(agb * 1e-6 * 0.01))
  ) +
    scale_x_discrete(
      name = "AGB Class (Mg ha-1)",
      labels = c("0-50", "50-100", "100-150", ">150")
    ) +
    scale_y_continuous(name = "AGB Stock (Tg)") +
    geom_bar(stat = "identity")

  ggsave(file.path(paths$figures, "AGB_distribution_x_AGBClass.png"), gg_agb_agb_class,
         width = 16, height = 12)

  ## Request 2: cumulative delta AGB by ecozone -------------------------------------------------

  agb_tifs <- fs::dir_ls(agbdirs, regexp = "ragb", recurse = TRUE)

  no_cores <- AGBtrends::getNumCores(length(agb_tifs))
  cl <- parallelly::makeClusterPSOCK(no_cores,
                                     default_packages = c("dplyr", "sf", "terra"),
                                     rscript_libs = .libPaths(),
                                     autoStop = TRUE)
  parallel::clusterExport(cl, varlist = c("no_cores", "paths", "pixelAreaHa", "years"))
  parallel::clusterEvalQ(cl, {
    terraOptions(
      tempdir = paths$terra,
      memmax = 25,
      memfrac = 0.6 / no_cores,
      progress = 1,
      verbose = TRUE
    )

    invisible(NULL)
  })

  fout <- file.path(paths$outputs, "sum_AGB_x_Ecozone.rds")

  if (file.exists(fout)) unlink(fout)

  do.call(rbind, parallel::parLapply(cl, agb_tifs, function(r) {
    targetCRS <- file.path(paths$outputs, "agb", "AB", "ragb_AB.tif") |>
      rast() |>
      crs()

    ez <- file.path("outputs", "studyArea_WBI", "WBI_studyArea.gpkg") |>
      st_read(quiet = TRUE) |>
      st_transform(targetCRS) |>
      select("ECOZONE") |>
      rasterize(rast(r, lyr = as.character(years[1])), field = "ECOZONE")

    yrs <- names(rast(r))

    x <- do.call(rbind, lapply(yrs, function(yr) {
      (rast(r, lyr = yr) * pixelAreaHa) |>
        zonal(ez, fun = "sum", na.rm = TRUE) |>
        mutate(Year = yr, .before = all_of(yr)) |>
        dplyr::rename(AGB = yr)
    }))

    return(x)
  })) |>
    saveRDS(file = fout)

  parallel::stopCluster(cl)

  ptab <- readRDS(fout) |>
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

  lapply(seq(length(ecozones)), function(i) {
    df <- filter(ptab, ECOZONE == ecozones[i]) |> ungroup()

    ## https://finchstudio.io/blog/ggplot-dual-y-axes/
    ## scale and shift variables calculated based on desired mins and maxes
    max_first <- max(df$dAGB, na.rm = TRUE) # Specify max of first y axis
    max_second <- max(df$AGB, na.rm = TRUE) # Specify max of second y axis
    min_first <- min(df$dAGB, na.rm = TRUE) # Specify min of first y axis
    min_second <- min(df$AGB, na.rm = TRUE) # Specify min of second y axis

    ## scale and shift variables calculated based on desired mins and maxes
    scale <- (max_second - min_second) / (max_first - min_first)
    shift <- min_first - min_second

    ## Function to scale secondary axis
    scale_function <- function(x, scale, shift) {
      return(x * scale - shift)
    }

    ## Function to scale secondary variable values
    inv_scale_function <- function(x, scale, shift) {
      return((x + shift) / scale)
    }

    gg_ptab_ez <- ggplot(
      data = df |>
        bind_rows(data.frame(ECOZONE = ecozones[i], AGB = NA, Year = years[1], dAGB = 0)),
      aes(x = Year, y = dAGB, weight = c(rep(1, nrow(df)), 100), color = "dAGB (Tg)")
    ) +
      ggtitle(ecozones[i]) +
      geom_smooth(method = "loess") +
      geom_point(aes(y = inv_scale_function(c(df$AGB, NA), scale, shift),
                     color = "Total AGB (Tg)"), pch = 20, cex = 0.75) +
      scale_x_continuous(breaks = seq.int(from = head(years, 1), to = tail(years, 1), by = 5),
                         labels = identity) +
      scale_y_continuous("Cumulative AGB change (Tg)",
                         limits = c(min_first, max_first * 1.1),
                         sec.axis = sec_axis(~ scale_function(., scale, shift),
                                             name = "Total AGB (Tg)")
      ) +
      geom_hline(yintercept = 0, lty = "dashed") +
      labs(color = "Units")

    ggsave(file.path(paths$figures, paste0("AGB_distribution_x_Year_", ecozones[i], ".png")),
           gg_ptab_ez, width = 16, height = 8)
  })

  # 8) Plot differences -------------------------------------------------------------------------
  ## 8 a) without disturbance mask --------------------------------------------------------------

  t_ref <- tref(timeint, years)

  ## i=1 corresponds to 31-year time series, i=2 corresponds to time interval t1 (1984-1988), and so on and so forth
  gg_71 <- plotZoneStats(
    files2plot = file.path(paths$summaries, "zoneStats_summary_WBI_ecozone.rds"),
    tref = t_ref
  )

  ggsave(
    file.path(paths$figures, "AGB_global_trends_WBI_ecozone_x_ageClass.png"),
    gg_71,
    width = 16,
    height = 12
  )

  ## x Ecozone x ageClass
  f2p <- file.path(paths$summaries) |>
    list.files(pattern = "zoneStats_summary_WBI_ecozone_", full.names = TRUE)

  gg_72 <- plotZoneStatsIntervals(
    files2plot = f2p,
    tref = t_ref,
    weighted = TRUE,
    xVar = "tp",
    groupVar = "ageClass",
    ptype = 1
  )

  ggsave(
    file.path(paths$figures, paste0("AGB_temporal_trends_x_ECOZONE_x_ageClass_", Sys.Date(), ".png")),
    gg_72,
    width = 16,
    height = 12
  )

  ## x ageClass x Ecozone
  f2p <- file.path(paths$summaries) |>
    list.files(pattern = "zoneStats_summary_WBI_ecozone_", full.names = TRUE)

  gg_73 <- plotZoneStatsIntervals(
    files2plot = f2p,
    tref = t_ref,
    weighted = TRUE,
    xVar = "tp",
    catVar = "ageClass",
    groupVar = "ECOZONE",
    ptype = 2,
    plotResult = FALSE
  ) |>
    cowplot::plot_grid(plotlist = _)

  ggsave(
    file.path(paths$figures, paste0("AGB_temporal_trends_x_ageClass_x_ECOZONE_", Sys.Date(), ".png")),
    gg_73,
    width = 10,
    height = 5
  )

  ## 8 b) with disturbance mask -----------------------------------------------------------------
  gg_74 <- plotZoneStats(
    files2plot = file.path(paths$summaries, "zoneStats_summary_WBI_distMask_ecozone.rds"),
    tref = t_ref
  )

  ## TODO: verify & adjust output filename
  # ggsave(
  #   file.path(paths$figures, paste0("AGB_temporal_trends_x_ECOZONE_distMask_", Sys.Date(), ".png")),
  #   gg_74
  # )

  ## x ageClass x Ecozone
  f2p <- file.path(paths$summaries) |>
    list.files(pattern = "zoneStats_summary_WBI_distMask_ecozone_", full.names = TRUE)

  gg_75 <- plotZoneStatsIntervals(
    files2plot = f2p,
    tref = t_ref,
    weighted = TRUE,
    xVar = "tp",
    catVar = "ageClass",
    groupVar = "ECOZONE",
    ptype = 2,
    plotResult = FALSE
  ) |>
    cowplot::plot_grid(plotlist = _) ## NOTE: only age class `0-24` here

  ggsave(
    file.path(paths$figures, paste0("AGB_temporal_trends_x_ECOZONE_distMask_", Sys.Date(), ".png")),
    gg_75
  )

  f2p <- file.path(paths$summaries) |>
    list.files(pattern = "WBI_distMask_ecozone", full.names = TRUE)

  gg_76 <- plotZoneStatsIntervals(
    files2plot,
    tref = t_ref
  ) ## NOTE: only age class `0-24` here

  ## TODO: verify & adjust output filename
  # ggsave(
  #   file.path(paths$figures, paste0("AGB_temporal_trends_x_ECOZONE_distMask_", Sys.Date(), ".png")),
  #   gg_76
  # )

  # 9) Test for significant differences between groups ------------------------------------------

  ## TODO

  # 10) upload results --------------------------------------------------------------------------

  drive_outputs <- as_id("19XHTS6V09ARYbVZ6hM1SQiUSLjWuO8-M")

  gid <- drive_ls(drive_outputs) |>
    filter(name == basename(paths$outputs)) |>
    as_id()

  if (length(gid) == 0) {
    gid <- drive_mkdir(name = basename(paths$outputs), path = drive_outputs) |>
      as_id()
  }

  files2upload <- c(
    f1, f2, f3, f4, f5, f6, f7,
    list.files(file.path(paths$figures), full.names = TRUE)
  )

  purrr::walk(files2upload, drive_put, path = gid)
}

# cleanup -------------------------------------------------------------------------------------
terra::tmpFiles(orphan = TRUE, remove = TRUE)
