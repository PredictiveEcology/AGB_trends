WBI_prepRasters <- function(type, studyAreas, years, paths, dst, cl = NULL) {
  if (is.null(cl)) {
    no_cores <- AGBtrends::getNumCores(length(studyAreas))
    mempercore <- as.integer(terra::free_RAM() / 1024^2 / no_cores) ## TODO: propagate this; find better ram detection fun
    cl <- parallelly::makeClusterPSOCK(
      no_cores,
      default_packages = c("sf", "stringr", "terra"),
      rscript_libs = .libPaths(),
      autoStop = TRUE
    )
  }

  parallel::clusterExport(cl, varlist = c("dst", "mempercore", "no_cores", "paths", "years"),
                          envir = environment())
  parallel::clusterEvalQ(cl, {
    terraOptions(
      tempdir = paths$terra,
      memmax = min(25, mempercore),
      memfrac = 0.8 / no_cores,
      progress = 1,
      verbose = TRUE
    )

    invisible(NULL)
  })

  out <- parallel::parLapply(cl, studyAreas, function(p) {
    wbi_rast_crs <- file.path(dirname(dst), "agb", p, paste0("ragb_", p, ".tif")) |>
      rast() |>
      crs()

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
    prov <- wbi_provs[wbi_provs$PREABBR == p, ] |>
      st_intersection(agb_tiles) |>
      st_union() |>
      vect()

    if (type == "disturbed") {
      allYears <- list.files(dst, full.names = TRUE, pattern = paste0(p, "_")) |>
        grep("_(binary_disturbed|disturbed_all)", x = _, invert = TRUE, value = TRUE) |>
        sort()

      stk <- rast(allYears) |>
        crop(prov) |>
        mask(prov)
      names(stk) <- years

      f_stk <- file.path(dst, paste0(p, "_", type, "_all.tif"))
      writeRaster(stk, f_stk, overwrite = TRUE)

      stk <- writeRaster(
        any(stk),
        filename = gsub("disturbed_all", "binary_disturbed", f_stk),
        overwrite = TRUE
      )

      try(unlink(allYears))
    } else {
      ## 'agb' or 'age'
      dst_p <- file.path(dst, p) |> reproducible::checkPath(create = TRUE)
      f <- file.path(dst_p, paste0("r", type, "_", p, ".tif"))

      stk <- rast(f) |>
        crop(prov) |>
        mask(prov)
      names(stk) <- years ## should already be the case, but enforce it anyway
      stk <- writeRaster(stk, f, overwrite = TRUE)
    }

    sources(stk)
  }) |>
    unlist()

  # parallel::stopCluster(cl)
  terra::tmpFiles(remove = TRUE)

  return(out)
}
