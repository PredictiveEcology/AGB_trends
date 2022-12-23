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
  repos = c(CRAN = "https://cran.rstudio.com"),
  Require.RPackageCache = "default" ## will use default package cache directory: `RequirePkgCacheDir()`
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

Require.version <- "PredictiveEcology/Require@development"
if (!"Require" %in% rownames(installed.packages(lib.loc = .libPaths()[1])) ||
    packageVersion("Require", lib.loc = .libPaths()[1]) < "0.2.5") {
  remotes::install_github(Require.version)
}

library(Require)

## temporarily until new Rcpp release on CRAN in early 2023 ----------------------------------------
options("Require.otherPkgs" = setdiff(getOption("Require.otherPkgs"), "Rcpp")) ## remove Rcpp from "forced source"
RcppVersionNeeded <- package_version("1.0.9.3")

RcppVersionAvail <- if (!"Rcpp" %in% rownames(installed.packages(lib.loc = .libPaths()[1]))) {
  package_version(data.table::as.data.table(available.packages())[Package == "Rcpp", Version])
} else {
  package_version(packageVersion("Rcpp", lib.loc = .libPaths()[1]))
}

if (RcppVersionAvail < RcppVersionNeeded) {
  Require(paste0("Rcpp (>= ", RcppVersionNeeded, ")"),  repos = "https://rcppcore.github.io/drat",
          require = FALSE, verbose = 1)
}
##

setLinuxBinaryRepo()

Require(c("PredictiveEcology/SpaDES.project@transition (>= 0.0.7.9003)", ## TODO: use development once merged
          "PredictiveEcology/SpaDES.config@development (>= 0.0.2.9050)",
          "PredictiveEcology/SpaDES.tools@development"),
        upgrade = FALSE, standAlone = TRUE)

modulePkgs <- unname(unlist(packagesInModules(modulePath = file.path(prjDir, "modules"))))
otherPkgs <- c(
  "archive", "details", "DBI", "s-u/fastshp", "logging", "RPostgres", "RSQLite"
)

Require(unique(c(modulePkgs, otherPkgs)), require = FALSE, standAlone = TRUE, upgrade = FALSE)

## NOTE: always load packages LAST, after installation above;
##       ensure plyr loaded before dplyr or there will be problems
Require(c("data.table", "plyr", "pryr", "SpaDES.core",
          "dplyr", "googledrive", "httr", "sessioninfo", "sf", "terra"),
        upgrade = FALSE, standAlone = TRUE)

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
      file.path("/mnt/scratch/achubaty", basename(prjDir))
    }
  } else if (.user == "trudolph") {
    if (grepl("for-cast[.]ca", .nodename)) {
      file.path("/mnt/scratch/trudolph", basename(prjDir))
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

# project options -----------------------------------------------------------------------------

# opts <- SpaDES.config::setProjectOptions(config)

opts <- options(
  reproducible.conn = SpaDES.config::dbConnCache(config$cacheDBtype),
  reproducible.destinationPath = prjPaths$inputPath
)

terra::terraOptions(tempdir = prjPaths$terraPath, todisc = TRUE)

quickPlot::dev.useRSGD(useRSGD = quickPlot::isRstudioServer())

# authenticate google user ---------------------------------------------------------------------

## allow authentication via email/oauth (interactive) or by token (non-interactive)
SpaDES.config::authGoogle(tryToken = "AGB_trends", tryEmail = config$googleUser)

# simulation ----------------------------------------------------------------------------------

do.call(SpaDES.core::setPaths, prjPaths) ## set paths for simulation

## define study area

## "Canada_Albers_Equal_Area_Conic" - no recognized EPSG code, using wkt:
targetCRS <- paste0("PROJCRS[\"Canada_Albers_Equal_Area_Conic\",\n",
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

myStudyArea <- switch(
  .mode,
  production = {
    bcrWBI <- Cache(prepInputs,
                    url = "https://www.birdscanada.org/download/gislab/bcr_terrestrial_shape.zip",
                    destinationPath = file.path(prjPaths$inputPath, "WBI"),
                    fun = "sf::st_read") %>%
      filter(BCR %in% c(4, 6:8))

    provsWBI <- geodata::gadm(country = "CAN", level = 1, path = file.path(prjPaths$inputPath, "WBI")) %>%
      st_as_sf() %>%
      st_transform(targetCRS) %>%
      filter(NAME_1 %in% c("British Columbia", "Alberta", "Saskatchewan", "Manitoba",
                           "Yukon", "Northwest Territories", "Nunavut"))

    studyAreaWBI <- Cache(postProcess, provsWBI, studyArea = bcrWBI, useSAcrs = TRUE, filename2 = NULL)

    gpkgFile <- file.path(prjPaths$outputPath, "WBI_studyArea.gpkg")
    if (!file.exists(gpkgFile)) {
      st_write(studyAreaWBI, dsn = gpkgFile, driver = "GPKG", delete_layer = TRUE)
    }

    studyAreaWBI
  },
  development = {
    ## random study area for development
    ctr1 <- sp::SpatialPoints(matrix(c(-113.530, 61.530), ncol = 2))
    crs(ctr1) <- "EPSG:4326"
    studyAreaRnd1 <- SpaDES.tools::randomStudyArea(ctr1, size = 3e10, seed = 42) %>%
      st_as_sf() %>%
      st_transform(targetCRS)

    studyAreaRnd1
  },
  testing = {
    ## random study area for testing (larger + diff location than devel)
    ctr2 <- sp::SpatialPoints(matrix(c(-126.95, 61.30), ncol = 2))
    crs(ctr2) <- "EPSG:4326"
    studyAreaRnd2 <- SpaDES.tools::randomStudyArea(ctr2, size = 6e10, seed = pi) %>%
      st_as_sf() %>%
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
      .studyAreaName = if (.mode == "production") "WBI" else "test"
    )
  ),
  objects = list(
    studyArea = myStudyArea
  ),
  modules = list("AGB_dataPrep"),
  paths = prjPaths
)
