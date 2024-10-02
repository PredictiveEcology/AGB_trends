WBI_standAge <- function(agbfiles, agefiles, studyareas, years, cl = NULL) {
  cores <- getNumCores(length(studyareas))

  if (is.null(cl)) {
    cl <- parallelly::makeClusterPSOCK(
      cores,
      default_packages = c("sf", "terra"),
      rscript_libs = .libPaths(),
      autoStop = TRUE
    )
    parallel::clusterExport(
      cl,
      varlist = c("agbfiles", "agefiles", "cores", "studyareas", "years"),
      envir = environment()
    )
    parallel::clusterEvalQ(cl, {
      terraOptions(
        memmax = min(25, as.integer(terra::free_RAM() / 1024^2)),
        memfrac = 0.8 / cores,
        progress = 0,
        verbose = TRUE,
        print = FALSE
      )
    })
  }

  outFiles <- parallel::parLapply(cl, seq(length(studyareas)), function(i) {
    sadir <- studyareas[i]
    sadir_out <- file.path(sadir) |>
      sub("inputs", "outputs", x = _) |>
      sub("/([A-Z][A-Z])_(.*)", "/AGB_WBI/\\2/\\1/", x = _)
    sadir_out_agb <- file.path(dirname(sadir_out), "agb", basename(sadir_out)) |>
      reproducible::checkPath(create = TRUE)
    sadir_out_age <- file.path(dirname(sadir_out), "age", basename(sadir_out)) |>
      reproducible::checkPath(create = TRUE)

    sa <- substr(basename(sadir), 1, 2)

    ## crop/mask to ABoVE tile areas
    wbi_rast_crs <- rast(agbfiles[1]) |> crs()

    agb_gpkg <- file.path("outputs", "studyArea_WBI", "ABoVE_AGB_study_area.gpkg")
    agb_tiles <- st_read(agb_gpkg, "tileset", quiet = TRUE) |>
      st_buffer(1e3) |>
      st_buffer(-1e3) |>
      st_union() |>
      st_make_valid() |>
      st_transform(wbi_rast_crs)

    wbi_provs <- file.path("outputs_wbi", "AGB_WBI", "WBI_studyArea_provs.gpkg") |>
      st_read(quiet = TRUE) |>
      st_transform(wbi_rast_crs)
    wbi_provs$PREABBR <- c("AB", "BC", "MB", "NT", "NT", "SK", "YT") ## NU is part of NT rasters
    prov <- wbi_provs[wbi_provs$PREABBR == sa, ] |>
      st_intersection(agb_tiles) |>
      st_union() |>
      vect()

    ## Import AGB rasters (0 = NA)
    ragb <- grep(sadir, agbfiles, value = TRUE) |>
      sort() |>
      rast() |>
      crop(prov) |>
      mask(prov)
    names(ragb) <- as.character(years)
    ragb <- writeRaster(ragb, file.path(sadir_out_agb, paste0("ragb_", sa, ".tif")), overwrite = TRUE)

    ## NOTE: age maps here are 'time since fire'
    ## - stand age as of year associated with raster layer (negative values become NA);
    rage <- grep(sadir, agefiles, value = TRUE) |>
      sort() |>
      rast() |>
      crop(prov) |>
      mask(prov) |>
      classify(rcl = cbind(-100L, 0L, NA), include.lowest = TRUE)
    names(rage) <- years
    rage <- writeRaster(rage, file.path(sadir_out_age, paste0("rage_", sa, ".tif")), overwrite = TRUE)

    gc()

    return(list(sources(ragb), sources(rage)))
  }) |>
    unlist(recursive = TRUE)

  parallel::stopCluster(cl)

  return(invisible(outFiles))
}
