## Date Created: Nov 17, 2022
## Auteurs: Tyler Rudolph, biologist M.Sc., CFS/NRCAN, Trade, Economics & Industry Branch
##          Alex M. Chubaty, PhD, FOR-CAST Research & Analytics
##
## Name of script : "00_ABoVE_data_prep.R"
##
## Description :
##   Estimate cell-specific stand age for raster tiles from the ABoVE AGB time series
##   using a combination of ABoVE Disturbance Agents and the CaNFIR stand age mosaic (kNN 2020)

# package installation and loading ------------------------------------------------------------
Require::Require(c("dplyr", "reproducible", "sf", "stringr", "terra",
                   "PredictiveEcology/AGBtrends (>= 0.0.3)"), upgrade = FALSE)

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
# options(parallelly.availableCores.fallback = 4L) ## set to limit the number of cores
no_cores <- AGBtrends::getNumCores() ## use up to half the number cores or fallback

terraOptions(tempdir = paths$terra, todisk = TRUE)

## define time intervals (year ranges between 1984-2014)
timeint <- list(t1 = 1:5, t2 = 6:10, t3 = 11:15, t4 = 16:20, t5 = 21:25, t6 = 26:31)
timeint_all <- timeint |> unlist() |> unname() |> list(all = _)

# 1) data import ------------------------------------------------------------------------------

agb_gpkg <- file.path(paths$outputs, "ABoVE_AGB_study_area.gpkg")
dstagnt_gpkg <- file.path(paths$outputs, "ABoVE_DistAgents_study_area.gpkg")
studyArea_gpkg <- file.path(paths$outputs, "WBI_studyArea.gpkg")

years <- 1984:2014

## 1.1) Import ABoVE product tiles ------------------------------------------------------------
agb_tiles <- st_read(agb_gpkg, "tileset")
dist_tiles <- st_read(dstagnt_gpkg, "tileset")

## 1.2) Index 30m ABG input rasters -----------------------------------------------------------
agbdsn <- file.path(paths$inputs, "ABoVE_AGB_30m", "data")
agbfiles <- list.files(agbdsn, pattern = "AGB_B", full.names = TRUE)

## 1.3) Index ABoVE disturbance history raster masks ------------------------------------------
distdsn <- file.path(paths$inputs, "ABoVE_ForestDisturbance_Agents", "data")
distfiles <- list.files(distdsn, pattern = ".tif", full.names = TRUE)

## 1.4) Import WBI study area -----------------------------------------------------------------
wbi <- st_read(studyArea_gpkg)

## 1.5) Reduce input file selection to tiles contained within WBI study area & align filename indices
agbfiles <- agbfiles[apply(relate(x = vect(agb_tiles), y = vect(wbi), relation = "intersects"), 1, any)]
tilenames <- file.path(paths$outputs, "tiles", str_sub(str_remove_all(basename(agbfiles), "ABoVE_AGB_"), end = -5L))
distfiles <- distfiles[match(basename(tilenames), str_sub(basename(distfiles), end = 7L))]

# 2) Estimate cell-specific stand age ---------------------------------------------------------
## using a combination of ABoVE Disturbance Agents and the CaNFIR stand age mosaic (kNN 2020)

ptime <- system.time({
  age_mosaic <- file.path(paths$inputs, "CaNFIR", "mosaic_age.tif")
  outfiles <- AGBtrends::ABovE_CaNFIR_standAge(agbfiles, distfiles, age_mosaic, tilenames) ## TODO: test
}) # 1.6 hours

## Create mask limiting input values to cells disturbed during time series --------------------
## (according to ABoVE)

no_cores <- AGBtrends::getNumCores(25L) ## TODO: why 25 here; RAM limitation?
mempercore <- as.integer(terra::free_RAM() / 1024^2 / cores) ## TODO: propagate mempercore throughout
cl <- parallelly::makeClusterPSOCK(
  no_cores,
  default_packages = c("stringr", "terra"),
  rscript_libs = .libPaths(),
  autoStop = TRUE
)

parallel::clusterExport(cl, varlist = "distdsn", envir = environment())
parallel::clusterEvalQ(cl, {
  terraOptions(
    tempdir = paths$terra,
    memmax = min(25, mempercore),
    memfrac = 0.8 / no_cores,
    progress = 1,
    verbose = TRUE
  )
})

parallel::parLapply(cl, distfiles, function(f) {
  writeRaster(
    any(rast(f)),
    filename = file.path(paths$outputs, "binary_disturbed", paste0(str_sub(basename(f), end = 7L), "_disturbed.tif")),
    overwrite = TRUE
  )
})

parallel::stopCluster(cl)

## build raster mosaic
dist_src <- list.files(file.path(paths$outputs, "binary_disturbed"), full.names = TRUE)
dist_mosaic <- file.path(paths$outputs, "mosaics")

f0a <- AGBtrends::buildMosaics("binary_disturbed", intervals = NULL, src = dist_src, dst = dist_mosaic)

# 3) Evaluate frequency of forested land cover types by ecozone -------------------------------
## * FOR IMPLEMENTATION AS MASK *

tileIDs <- dir(file.path(paths$outputs, "tiles"))

## Create land cover mosaic
lc_src <- list.files(file.path(paths$inputs, "ABoVE_LandCover"), full.names = TRUE)
lc_dst <- file.path(paths$outputs, "mosaics")

## TODO: confirm we want simplified version
f0b <- AGBtrends::buildMosaics("LandCover_Simplified", intervals = timeint_all, src = lc_src, dst = lc_dst)

## Rasterize ecozones to fit
ecoRast_tif <- file.path(paths$terra, "ecoRast.tif")

ecoRast <- terra::rasterize(
  x = wbi,
  y = rast(file.path(paths$outputs, "mosaics", "landcover_simplified_mosaic.tif")),
  field = "ECOZONE",
  filename = ecoRast_tif
)

no_cores <- AGBtrends::getNumCores(6L) ## TODO: why 6 here; number of ecoRast_tif levels
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
})

system.time({
  ftab <- parLapply(cl, cats(rast(ecoRast_tif))[[1]]$value, function(ezone) {
    return(freq(mask(
      rast(file.path(paths$outputs, "landcover_simplified_mosaic.tif")),
      rast(ecoRast_tif),
      maskvalues = ezone,
      inverse = TRUE
    )))
  })
})

stopCluster(cl)

names(ftab) <- cats(rast(ecoRast_tif))[[1]]$ECOZONE
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

# cleanup -------------------------------------------------------------------------------------
unlink(paths$terra, recursive = TRUE)
