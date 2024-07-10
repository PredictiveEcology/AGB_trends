# project basics ------------------------------------------------------------------------------

if (file.exists("~/.Renviron")) readRenviron("~/.Renviron") ## GITHUB_PAT
if (file.exists("AGB_trends.Renviron")) readRenviron("AGB_trends.Renviron") ## database credentials

.ncores <- min(parallel::detectCores() / 2, 32L)
.nodename <- Sys.info()[["nodename"]]
.user <- Sys.info()[["user"]]

if (exists(".mode", .GlobalEnv)) {
  stopifnot(.mode %in% c("development", "testing", "production"))
} else {
  .mode <- "development"
}

## packages, paths and options --------------------------------------------------------------------------

library(data.table)
library(plyr)
library(pryr)
library(googledrive)
library(httr)
library(sf)
library(terra)

library(reproducible)
library(SpaDES.config)
library(SpaDES.core)

prjDir <- SpaDES.config::findProjectPath()

stopifnot(identical(prjDir, getwd()))

options(
  Ncpus = .ncores,
  repos = c(CRAN = "https://cloud.r-project.org")
)

workflowtools::check_project_packages(prjDir)

# configure project ---------------------------------------------------------------------------

## TODO

config <- list(
  cacheDBtype = switch(.user,
                       #achubaty = "postgresql",
                       "sqlite"),
  googleUser = switch(.user,
                      achubaty = "achubaty@for-cast.ca",
                      trudolph = "tylerdrudolph@gmail.com")
)

# project paths -------------------------------------------------------------------------------

prjPaths <- list(
  cachePath = "cache",
  inputPath = "inputs",
  modulePath = "modules",
  outputPath = "outputs",
  scratchPath = file.path(tempdir(), "scratch", basename(prjDir))
)
prjPaths$rasterPath <- checkPath(file.path(prjPaths$scratchPath, "raster"), create = TRUE)
prjPaths$terraPath <- checkPath(file.path(prjPaths$scratchPath, "terra"), create = TRUE)

do.call(SpaDES.core::setPaths, prjPaths) ## set paths for simulation

# project options -----------------------------------------------------------------------------

# opts <- SpaDES.config::setProjectOptions(config)

opts <- options(
  reproducible.conn = SpaDES.config::dbConnCache(config$cacheDBtype),
  reproducible.destinationPath = prjPaths$inputPath
)

terra::terraOptions(tempdir = prjPaths$terraPath, todisk = TRUE)

quickPlot::dev.useRSGD(useRSGD = quickPlot::isRstudioServer())

# authenticate google user ---------------------------------------------------------------------

## allow authentication via email/oauth (interactive) or by token (non-interactive)
SpaDES.config::authGoogle(tryToken = "forprod", tryEmail = config$googleUser)

# simulation ----------------------------------------------------------------------------------

## define study area

sf_use_s2(FALSE)

targetCRS <- AGBtrends::Canada_Albers_Equal_Area_Conic

myStudyArea <- switch(
  .mode,
  production = {
    bcrWBI <- Cache(prepInputs,
                    url = "https://www.birdscanada.org/download/gislab/bcr_terrestrial_shape.zip",
                    destinationPath = file.path(prjPaths$inputPath, "WBI"),
                    fun = "sf::st_read") |>
      st_transform(targetCRS) |>
      filter(BCR %in% c(4, 6:8))

    provsWBI <- geodata::gadm(country = "CAN", level = 1, path = file.path(prjPaths$inputPath, "WBI")) |>
      st_as_sf() |>
      st_transform(targetCRS) |>
      filter(NAME_1 %in% c("British Columbia", "Alberta", "Saskatchewan", "Manitoba",
                           "Yukon", "Northwest Territories", "Nunavut"))

    studyAreaWBI <- Cache(postProcess, provsWBI, studyArea = bcrWBI, useSAcrs = TRUE, filename2 = NULL) |>
      st_union() |>
      st_buffer(0)

    gpkgFile <- file.path(prjPaths$outputPath, "WBI_studyArea.gpkg")
    if (!file.exists(gpkgFile)) {
      st_write(studyAreaWBI, dsn = gpkgFile, driver = "GPKG", delete_layer = TRUE)
    }

    studyAreaWBI
  },
  development = {
    ## random study area for development
    studyAreaRnd1 <- terra::vect(cbind(x = -113.530, y = 61.530), crs = "EPSG:4326") |>
      SpaDES.tools::randomStudyArea(size = 3e10, seed = 42) |>
      st_as_sf() |>
      st_transform(targetCRS)

    studyAreaRnd1
  },
  testing = {
    ## random study area for testing (larger + diff location than devel)
    studyAreaRnd2 <- terra::vect(cbind(x = -126.95, y = 61.30), crs = "EPSG:4326") |>
      SpaDES.tools::randomStudyArea(size = 6e10, seed = pi) |>
      st_as_sf() |>
      st_transform(targetCRS)

    studyAreaRnd2
  }
)

## start the simulation
mySim <- simInitAndSpades(
  times = list(start = 0, end = 1),
  params = list(
    .GlobalEnv = list(
      analysisZonesType = "ecozone",
      .plots = c("screen", "png", "raw"),
      .studyAreaName = ifelse(.mode == "production", "WBI", "test"),
      .useParallel = TRUE
    ),
    AGB_dataPrep = list(
      ## TODO
    ),
    AGB_analyses = list(
      ## TODO
    )
  ),
  objects = list(
    studyArea = myStudyArea
  ),
  modules = list(
    "AGB_dataPrep",
    "AGB_analyses"
  ),
  paths = prjPaths
)
