# packages ------------------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(reproducible)
library(sf)
library(stringr)
library(terra)

# global parameters for project setup ---------------------------------------------------------

projName <- workflowtools::findProjectName()
studyAreaName <- "AGB_WBI"
user <- Sys.info()[["user"]]

## specify which WBI runs to use for AGB analyses

climateGCM <- "CanESM5"
climateSSP <- 370
climateScenario <- paste0(climateGCM, "_SSP", climateSSP)

allReps <- sprintf("%02d", 1:5)
thisRep <- allReps[1]

allStudyAreas <- c("BC", "AB", "NT", "SK", "YT") ## MB doesn't overlap with ABoVE (see #5 below)

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
paths$terra <- checkPath(file.path(paths$scratch, "terra"), create = TRUE)

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

# 1) locate WBI simulated AGB rasters ---------------------------------------------------------

agbrasters <- file.path(paths$inputs, allOutputDirs) |>
  rep(length(years)) |>
  sort() |>
  file.path(paste0("simulatedBiomassMap_", years, "_year", years, ".tif"))

stopifnot(all(file.exists(agbrasters)))

# 2) locate WBI simulated burn (disturbance) maps ---------------------------------------------

distrasters <- file.path(paths$inputs, allOutputDirs) |>
  rep(length(years)) |>
  sort() |>
  file.path(paste0("burnMap_", years, "_year", years, ".tif"))

stopifnot(all(file.exists(agerasters)))

# 3) locate WBI simulated time since fire maps  -----------------------------------------------

agerasters <- file.path(paths$inputs, allOutputDirs, "postprocess") |>
  rep(length(years)) |>
  sort() |>
  file.path(paste0("timeSinceFire_", years, ".tif"))

stopifnot(all(file.exists(distrasters)))

# 4) create or locate WBI simulated landcover maps --------------------------------------------

## TODO: do we need LCC maps from WBI? use Tati's LCC fun

# 5) OPTIONAL: Visually compare available tiles between ABoVE products and WBI study area -----

## previously-created ABoVE and WBI study area polygons
agb_gpkg <- file.path("outputs", "studyArea_WBI", "ABoVE_AGB_study_area.gpkg")
dstagnt_gpkg <- file.path("outputs", "studyArea_WBI", "ABoVE_DistAgents_study_area.gpkg")
studyArea_gpkg <- file.path("outputs", "studyArea_WBI", "WBI_studyArea.gpkg")
wbi_provs_gpkg <- file.path(paths$outputs, "WBI_studyArea_provs.gpkg")

agb_sa <- st_read(agb_gpkg, "study_area", quiet = TRUE) |>
  st_transform(targetCRS)
agb_tiles <- st_read(agb_gpkg, "tileset", quiet = TRUE) |>
  st_transform(targetCRS)

dist_sa <- st_read(dstagnt_gpkg, "study_area", quiet = TRUE) |>
  st_transform(targetCRS)
dist_tiles <- st_read(dstagnt_gpkg, "tileset", quiet = TRUE) |>
  st_transform(targetCRS)

if (FALSE) {
  ## only needs to be run once
  bcrzip <- "https://www.birdscanada.org/download/gislab/bcr_terrestrial_shape.zip"

  bcrWB <- Cache(
    prepInputs,
    url = bcrzip,
    cachePath = paths$cache,
    destinationPath = "inputs", ## use AGB inputs dir, not WBI outputs
    targetCRS = targetCRS,
    fun = "sf::st_read"
  ) |>
    filter(BCR %in% c(4, 6:8))

  provsWB <- Cache(
    prepInputs,
    url = "https://www12.statcan.gc.ca/census-recensement/2011/geo/bound-limit/files-fichiers/2016/lpr_000b16a_e.zip",
    cachePath = paths$cache,
    destinationPath = "inputs", ## use AGB inputs dir, not WBI outputs
    targetFile = "lpr_000b16a_e.shp",
    alsoExtract = "similar",
    targetCRS = targetCRS,
    fun = "sf::st_read"
  ) |>
    filter(PREABBR %in% c("B.C.", "Alta.", "Sask.", "Man.", "Y.T.", "N.W.T.", "Nvt.")) |>
    st_cast("MULTIPOLYGON")

  wbi_provs <- bcrWB |>
    st_crop(provsWB[provsWB$PREABBR != "Nvt.", ]) |>
    st_intersection(provsWB) |>
    group_by(PRENAME) |>
    summarize() |>
    st_buffer(0.0) |>
    st_make_valid()
  wbi_provs$PREABBR <- c("AB", "BC", "MB", "NT", "NT", "SK", "YT") ## NU is part of NT rasters
  st_write(wbi_provs, dsn = wbi_provs_gpkg, driver = "GPKG", delete_layer = TRUE)

  rm(bcrzip, bcrWB, provsWB)
}

wbi <- st_read(studyArea_gpkg, quiet = TRUE) |>
  st_transform(targetCRS)
wbi_provs <- st_read(wbi_provs_gpkg, quiet = TRUE) |>
  st_transform(targetCRS)

## comparison plot
plot_studyArea_tiles <- function() {
  plot(wbi_provs |> st_geometry())
  plot(filter(dist_tiles, apply(relate(x = vect(dist_tiles), y = vect(wbi_provs), relation = "intersects"), 1, any)) |>
         st_geometry(), border = "red", add = TRUE)
  plot(filter(agb_tiles, apply(relate(x = vect(agb_tiles), y = vect(wbi_provs), relation = "intersects"), 1, any)) |>
         st_geometry(),
    border = "lightblue", add = TRUE
  )
} ## no grey tile overlap with MB, NU

gg_tiles <- cowplot::plot_grid(plot_studyArea_tiles)

ggsave(
  filename = file.path(paths$outputs, "figures", paste0("ABoVE_tiles_", studyAreaName, ".png")),
  plot = gg_tiles,
  create.dir = TRUE,
  height = 8, width = 12
)

# cleanup -------------------------------------------------------------------------------------
terra::tmpFiles(orphan = TRUE, remove = TRUE)
unlink(paths$terra, recursive = TRUE)
