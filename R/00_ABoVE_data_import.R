## Date Created: Nov 17, 2022
## Auteurs: Tyler Rudolph, biologist M.Sc., CFS/NRCAN, Trade, Economics & Industry Branch
##          Alex M. Chubaty, PhD, FOR-CAST Research & Analytics
##
## Name of script : "00_ABoVE_data_import.R"
##
## Description :
##   Import 1) above-ground biomass (AGB) and 2) disturbance (year/type) time series datasets
##   created through the Arctic-Boreal Vulnerability Experiment (ABoVE) research project
##   <https://daac.ornl.gov/cgi-bin/dataset_lister.pl?p=34>;
##   3) Import CaNFIR kNN stand age estimation (2020)
##   4) Create spatial reference polygons corresponding to individual ABoVE tiles

# package installation and loading ------------------------------------------------------------
Require::Install(c("cowplot", "gridGraphics"), upgrade = FALSE)
Require::Require(
  c("dplyr", "ggplot2", "googledrive", "reproducible", "sf", "stringr", "terra",
    "PredictiveEcology/AGBtrends (>= 0.0.4)"),
  upgrade = FALSE
)

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

## set the max number of cores to use for parallel computations
no_cores <- AGBtrends::getNumCores()

terraOptions(tempdir = paths$terra, todisk = TRUE)

auth_json <- list.files(pattern = "forprod-.*[.]json")
if (length(auth_json) == 0) {
  stop("Google Service Account token file (.json) not found.")
}

drive_auth(path = auth_json)

## ABoVE default CRS (Canada_Albers_Equal_Area_Conic)
targetCRS <- AGBtrends::Canada_Albers_Equal_Area_Conic |> crs()

# 1) download tiled ABoVE AGB rasters ---------------------------------------------------------
agbtiles <- drive_ls(as_id("1CdJ4t_Ja0qcQgk5DIkFi8QSYQVcP9EV4"))

## 1 a) setup parallel processing (one thread per tile) ---------------------------------------
cl <- parallelly::makeClusterPSOCK(
  no_cores,
  default_packages = c("googledrive"),
  rscript_libs = .libPaths(),
  autoStop = TRUE
)
parallel::clusterExport(cl, varlist = c("agbtiles", "auth_json", "paths"))

## 1 b) import raster tiles -------------------------------------------------------------------
parLapply(cl, seq(nrow(agbtiles)), function(m) {
  drive_auth(path = auth_json)
  retry(quote({
    httr::with_config(config = httr::config(http_version = 2), {
      drive_download(
        file = as_id(agbtiles[m, ]),
        path = file.path(paths$inputs, "ABoVE_AGB_30m", "data", agbtiles$name[m]),
        overwrite = TRUE
      ) ## TODO: use check files instead + check sums
    })
  }))
})

stopCluster(cl)

# 2) download tiled ABoVE Forest Disturbance Agents -------------------------------------------
distiles <- drive_ls(as_id("1CNalAGmw9fO0-TuMMqLt3HrMNTj6cnER"))

## 2 a) setup parallel processing (one thread per tile) ---------------------------------------
cl <- makeCluster(no_cores)
clusterExport(cl, varlist = c("distiles", "auth_json", "paths"))

## 2 b) import raster tiles -------------------------------------------------------------------
parLapply(cl, seq(nrow(distiles)), function(m) {
  library(googledrive)
  drive_auth(path = auth_json)
  retry(quote({
    httr::with_config(config = httr::config(http_version = 2), {
      drive_download(
        file = as_id(distiles[m, ]),
        path = file.path(paths$inputs, "ABoVE_ForestDisturbance_Agents", "data", distiles$name[m]),
        overwrite = TRUE
      )
    })
  }))
})

parallel::stopCluster(cl)

# 3) download CaNFIR stand age product + supporting file types --------------------------------

## create destination folder if non-existent
checkPath(file.path(paths$inputs, "CaNFIR"), create = TRUE)

## authenticate w/ Google Cloud
preProcess(
  url = "https://drive.google.com/file/d/1Thvxc9I8d1DzE7_hMovKESaCQGuHEYKy",
  targetFile = "mosaic_age.tif",
  fun = "terra::rast",
  alsoExtract = "similar",
  destinationPath = file.path(paths$inputs, "CaNFIR"),
  cacheRepo = paths$cache,
  userTags = c("bigRasters", "CaNFIR"), ## keywords to use in the cache db to make it easy to find later
  overwrite = TRUE
)

# 4) Import and pre-process ABoVE Landcover time series product (1984 - 2014) -----------------

## 4a) Download ABoVE Landcover time series data and metadata ---------------------------------

## Set up access to EarthData webserver (in bash)
## Set up your ~/.netrc file as listed here:
##   https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+cURL+And+Wget

# cd ~
# touch .netrc
# echo "machine urs.earthdata.nasa.gov login USERNAME password PASSWORD" > .netrc
# chmod 0600 .netrc
# touch .urs_cookies

# library(httr)
# netrc_path <- "~/.netrc"
# cookie_path <- "~/.urs_cookies"
# downloaded_file_path <- "~/test"
# set_config(config(
#   followlocation = 1,
#   netrc = 1,
#   netrc_file = netrc_path,
#   cookie = cookie_path,
#   cookiefile = cookie_path,
#   cookiejar = cookie_path
# ))

tileIDs <- file.path(paths$inputs, "ABoVE_AGB_30m", "data") |>
  dir() |>
  stringr::str_sub(start = -11L, end = -5L) |>
  unique()

url_prefix_data <- "https://daac.ornl.gov/daacdata/above/Annual_Landcover_ABoVE/data/"

## Get metadata
httr::GET(
  url = "https://daac.ornl.gov/daacdata/above/Annual_Landcover_ABoVE/comp/Annual_Landcover_ABoVE.pdf",
  httr::write_disk(file.path(paths$inputs, "ABoVE_Landcover", "Annual_Landcover_ABoVE.pdf"), overwrite = TRUE)
)

httr::GET(
  url = paste0(url_prefix_data, "accuracy_assess_1984-2014.csv"),
  httr::write_disk(file.path(paths$inputs, "ABoVE_Landcover", "accuracy_assess_1984-2014.csv"), overwrite = TRUE)
)

httr::GET(
  url = paste0(url_prefix_data, "accuracy_summary_1984-2014.csv"),
  httr::write_disk(file.path(paths$inputs, "ABoVE_Landcover", "accuracy_summary_1984-2014.csv"))
)

## Download individually tiled rasters

lc_files <- paste0(url_prefix_data, "ABoVE_LandCover_", tileIDs, ".tif")
lcs_files <- paste0(url_prefix_data, "ABoVE_LandCover_Simplified_", tileIDs, ".tif")

lc_dir <- file.path(paths$inputs, "ABoVE_LandCover")
lcs_dir <- file.path(paths$inputs, "ABoVE_LandCover") ## using same dir, since same dir used at url

lapply(lc_files, function(f) {
  httr::GET(
    url = f,
    httr::write_disk(file.path(lc_dir, basename(f)), overwrite = TRUE)
  )
})

lapply(lcs_files, function(f) {
  httr::GET(
    url = f,
    httr::write_disk(file.path(lcs_dir, basename(f)), overwrite = TRUE)
  )
})

# 5) Create polygons from ABoVE raster tiles --------------------------------------------------

## 5 a) Create vector polygons corresponding to ABoVE AGB raster tiles ------------------------
dsn <- file.path(paths$inputs, "ABoVE_AGB_30m", "data")
agb_tifs <- list.files(dsn, pattern = "AGB_B")
agb_gpkg <- file.path(paths$outputs, "ABoVE_AGB_study_area.gpkg")

## individual tile bounding boxes
boxes <- do.call(rbind, lapply(agb_tifs, function(x) as.vector(ext(rast(file.path(dsn, x))))))

## wide area bounding box
Abox <- c(apply(boxes, 2, min)[c(1, 3)], apply(boxes, 2, max)[c(2, 4)])[c(1, 3, 2, 4)] |>
  st_bbox() |>
  st_as_sfc() |>
  st_as_sf() |>
  st_set_crs(targetCRS)

st_write(
  merge(Abox, data.frame(description = "ABoVE_AGB_study_area")),
  dsn = agb_gpkg,
  layer = "study_area",
  driver = "GPKG", delete_layer = TRUE
)

## append individual tiles
vtiles <- do.call(rbind, lapply(1:nrow(boxes), function(m) {
  c(boxes[m, 1], boxes[m, 3], boxes[m, 2], boxes[m, 4]) |>
    st_bbox() |>
    st_as_sfc() |>
    st_as_sf(data.frame(tile_name = agb_tifs[m]))
})) |>
  st_set_crs(targetCRS)

st_write(
  vtiles,
  dsn = agb_gpkg,
  layer = "tileset", driver = "GPKG", delete_layer = TRUE
)

## 5 b) Create vector polygons corresponding to ABoVE disturbance history raster tiles --------
dsn <- file.path(paths$inputs, "ABoVE_ForestDisturbance_Agents", "data")
dstagnt_tifs <- list.files(dsn, pattern = ".tif")
dstagnt_gpkg <- file.path(paths$outputs, "ABoVE_DistAgents_study_area.gpkg")

## individual tile bounding boxes
boxes <- do.call(rbind, lapply(dstagnt_tifs, function(x) as.vector(ext(rast(file.path(dsn, x))))))

## wide area bounding box
Abox <- c(apply(boxes, 2, min)[c(1, 3)], apply(boxes, 2, max)[c(2, 4)])[c(1, 3, 2, 4)] |>
  st_bbox() |>
  st_as_sfc() |>
  st_as_sf() |>
  st_set_crs(targetCRS)

st_write(
  merge(Abox, data.frame(description = "ABoVE_ForestDisturbance_Agents_study_area")),
  dsn = dstagnt_gpkg,
  layer = "study_area", driver = "GPKG", delete_layer = TRUE
)

## append individual tiles
vtiles <- do.call(rbind, lapply(1:nrow(boxes), function(m) {
  c(boxes[m, 1], boxes[m, 3], boxes[m, 2], boxes[m, 4]) |>
    st_bbox() |>
    st_as_sfc() |>
    st_as_sf(data.frame(tile_name = dstagnt_tifs[m]))
})) |>
  st_set_crs(targetCRS)

st_write(vtiles, dsn = dstagnt_gpkg, layer = "tileset", driver = "GPKG", delete_layer = TRUE)

# 6) Import WBI study area --------------------------------------------------------------------

bcrzip <- "https://www.birdscanada.org/download/gislab/bcr_terrestrial_shape.zip"

bcrWB <- Cache(
  prepInputs,
  url = bcrzip,
  cacheRepo = paths$cache,
  destinationPath = paths$inputs,
  targetCRS = targetCRS,
  fun = "sf::st_read"
) |>
  filter(BCR %in% c(4, 6:8))

provsWB <- Cache(
  prepInputs,
  url = "https://www12.statcan.gc.ca/census-recensement/2011/geo/bound-limit/files-fichiers/2016/lpr_000b16a_e.zip",
  cacheRepo = paths$cache,
  destinationPath = paths$inputs,
  targetFile = "lpr_000b16a_e.shp",
  alsoExtract = "similar",
  targetCRS = targetCRS,
  fun = "sf::st_read"
) |>
  filter(PREABBR %in% c("B.C.", "Alta.", "Sask.", "Man.", "Y.T.", "N.W.T.", "Nvt.")) |>
  st_cast("MULTIPOLYGON")

## crop WBI to western Cdn provinces and save to file
studyArea_gpkg <- file.path(paths$outputs, "WBI_studyArea.gpkg")
studyArea <- bcrWB |>
  st_crop(provsWB[provsWB$PREABBR != "Nvt.", ]) |>
  st_intersection(provsWB) |>
  group_by(BCR, Label) |>
  summarize() |>
  AGBtrends::createAnalysisZones(targetCRS, paths$inputs) |>
  st_write(dsn = studyArea_gpkg, driver = "GPKG", delete_layer = TRUE)

## OPTIONAL: Visually compare available tiles between ABoVE products and WBI study area -------

## 1) AGB tiles
agb_sa <- st_read(agb_gpkg, "study_area")
agb_tiles <- st_read(agb_gpkg, "tileset")

## 2) Disturbance history tiles
dist_sa <- st_read(dstagnt_gpkg, "study_area")
dist_tiles <- st_read(dstagnt_gpkg, "tileset")

## 3) WBI study area
wbi <- st_read(studyArea_gpkg)

## 4) Comparison plot

plot_studyArea_tiles <- function() {
  plot(wbi |> st_geometry())
  plot(filter(dist_tiles, apply(relate(x = vect(dist_tiles), y = vect(wbi), relation = "intersects"), 1, any)) |> st_geometry(), border = "red", add = TRUE)
  plot(filter(agb_tiles, apply(relate(x = vect(agb_tiles), y = vect(wbi), relation = "intersects"), 1, any)) |> st_geometry(),
    border = "lightblue", add = TRUE
  )
}

gg_tiles <- cowplot::plot_grid(plot_studyArea_tiles)

ggsave(file.path(paths$outputs, "figures", paste0("ABoVE_tiles_", studyAreaName, ".png")), gg_tiles, width = 12, height = 8)

# cleanup -------------------------------------------------------------------------------------
unlink(paths$terra, recursive = TRUE)
