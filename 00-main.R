# project basics ------------------------------------------------------------------------------

if (file.exists("~/.Renviron")) readRenviron("~/.Renviron") ## GITHUB_PAT
if (file.exists("BC_HRV.Renviron")) readRenviron("BC_HRV.Renviron") ## database credentials

.ncores <- min(parallel::detectCores() / 2, 32L)
.nodename <- Sys.info()[["nodename"]]
.user <- Sys.info()[["user"]]

######

if (exists(".mode", .GlobalEnv)) {
  stopifnot(.mode %in% c("development", "testing", "production"))
} else {
  .mode <- "development"
}

#####

prjDir <- switch(.user,
                 achubaty  = "~/GitHub/AGB_trends",
                 trudolph = "~/GitHub/AGB_trends",
                 "~/GitHub/AGB_trends")

stopifnot(identical(normalizePath(prjDir), getwd())) ## ensure we're working in the project directory

## set new temp dir in scratch directory (existing /tmp too small for large callr ops)
## see https://github.com/r-lib/callr/issues/172
if (grepl("for-cast[.]ca", .nodename) && !grepl("larix", .nodename)) {
  oldTmpDir <- tempdir()
  newTmpDir <- file.path("/mnt/scratch", .user, basename(prjDir), "tmp")
  if (!dir.exists(newTmpDir)) dir.create(newTmpDir, recursive = TRUE)
  newTmpDir <- tools::file_path_as_absolute(newTmpDir)
  Sys.setenv(TMPDIR = newTmpDir)
  unlink(oldTmpDir, recursive = TRUE)
  tempdir(check = TRUE)
}

options(
  Ncpus = .ncores,
  repos = c(CRAN = "https://cloud.r-project.org")
)

# install and load packages -------------------------------------------------------------------

pkgDir <- file.path(tools::R_user_dir(basename(prjDir), "data"), "packages",
                    version$platform, getRversion()[, 1:2])
dir.create(pkgDir, recursive = TRUE, showWarnings = FALSE)
.libPaths(pkgDir, include.site = FALSE)
message("Using libPaths:\n", paste(.libPaths(), collapse = "\n"))

if (!"remotes" %in% rownames(installed.packages(lib.loc = .libPaths()[1]))) {
  install.packages("remotes")
}

Require.version <- "PredictiveEcology/Require@v0.3.1" ## use CRAN version
if (!"Require" %in% rownames(installed.packages(lib.loc = .libPaths()[1])) ||
    packageVersion("Require", lib.loc = .libPaths()[1]) != "0.3.1") {
  remotes::install_github(Require.version)
}

library(Require)

setLinuxBinaryRepo()

Require(c(
  "PredictiveEcology/SpaDES.project@transition (>= 0.0.7.9003)", ## TODO: use development once merged
  "PredictiveEcology/SpaDES.config@development (>= 0.0.2.9050)",
  "PredictiveEcology/SpaDES.tools@development"
), standAlone = TRUE, upgrade = FALSE)

modulePkgs <- unname(unlist(packagesInModules(modulePath = file.path(prjDir, "modules"))))
otherPkgs <- c(
  "archive", "details", "DBI", "s-u/fastshp", "httpuv", "logging", "RPostgres", "RSQLite"
)

Install(unique(c(modulePkgs, otherPkgs)), standAlone = TRUE, upgrade = FALSE)

## NOTE: always load packages LAST, after installation above;
##       ensure plyr loaded before dplyr or there will be problems
Require(c("plyr", "dplyr",
          "data.table", "googledrive", "httr", "pryr", "sessioninfo", "sf", "terra",
          "SpaDES.core"),
        standAlone = TRUE, upgrade = FALSE)

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
  scratchPath = if (.user == "achubaty") {
    if (.nodename == "larix.for-cast.ca") {
      file.path(tempdir(), "scratch", basename(prjDir))
    } else {
      file.path("/mnt/scratch", .user, basename(prjDir))
    }
  } else if (.user == "trudolph") {
    if (grepl("for-cast[.]ca", .nodename)) {
      file.path("/mnt/scratch", .user, basename(prjDir))
    } else {
      "scratch"
    }
  } else {
    ## default for unknown user
    file.path(tempdir(), "scratch", baseName(prjDir))
  }
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
    AGB_dataPrep = list(
      analysisZonesType = "ecozone",
      .plots = c("screen", "png", "raw"),
      .studyAreaName = if (.mode == "production") "WBI" else "test",
      .useParallel = TRUE
    ),
    AGB_analyses = list(
      analysisZonesType = "ecozone",
      .plots = c("screen", "png", "raw"),
      .studyAreaName = if (.mode == "production") "WBI" else "test",
      .useParallel = TRUE
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
