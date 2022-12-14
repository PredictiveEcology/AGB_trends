# project basics ------------------------------------------------------------------------------

if (file.exists("~/.Renviron")) readRenviron("~/.Renviron") ## GITHUB_PAT
if (file.exists("BC_HRV.Renviron")) readRenviron("BC_HRV.Renviron") ## database credentials

.ncores <- min(parallel::detectCores() / 2, 32L)
.nodename <- Sys.info()[["nodename"]]
.user <- Sys.info()[["user"]]

prjDir <- switch(.user,
                 achubaty  = "~/GitHub/AGB_trends",
                 trudolph = "~/GitHub/AGB_trends",
                 "~/GitHub/AGB_trends")

stopifnot(identical(normalizePath(prjDir), getwd())) ## ensure we're working in the project directory

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
          "PredictiveEcology/SpaDES.config@development (>= 0.0.2.9050)"),
        upgrade = FALSE, standAlone = TRUE)

modulePkgs <- unname(unlist(packagesInModules(modulePath = file.path(prjDir, "modules"))))
otherPkgs <- c(
  "archive", "details", "DBI", "s-u/fastshp", "logging", "RPostgres", "RSQLite"
)

Require(unique(c(modulePkgs, otherPkgs)), require = FALSE, standAlone = TRUE, upgrade = FALSE)

## NOTE: always load packages LAST, after installation above;
##       ensure plyr loaded before dplyr or there will be problems
Require(c("data.table", "plyr", "pryr", "SpaDES.core",
          "dplyr", "googledrive", "httr", "sf", "sessioninfo", "terra"),
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
  st_write(studyAreaWBI, dsn = gpkgFile, driver = "GPKG")
}

## start the simulation
mySim <- simInitAndSpades(
  times = list(start = 0, end = 1),
  params = list(
    AGB_dataPrep = list(
      analysisZonesType = "ecozone",
      .plots = c("screen", "png", "raw")
    )
  ),
  objects = list(
    studyArea = studyAreaWBI
  ),
  modules = list("AGB_dataPrep"),
  paths = prjPaths
)
