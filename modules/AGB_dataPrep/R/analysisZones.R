createAnalysisZones <- function(studyArea, targetCRS, destinationPath) {
  if (!is(studyArea, "Spatial")) {
    studyArea <- as_Spatial(studyArea)
  }

  urlList <- list(
    ecodistrict = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/district/ecodistrict_shp.zip",
    ecoregion = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/region/ecoregion_shp.zip",
    ecoprovince = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/province/ecoprovince_shp.zip",
    ecozone = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/zone/ecozone_shp.zip"
  )

  eco <- lapply(urlList, function(url) {
    suppressWarnings({
      prepInputs(url = url,
                 destinationPath = destinationPath,
                 studyArea = studyArea,
                 fun = "sf::st_read",
                 overwrite = TRUE) |>
        st_transform(targetCRS) |>
        st_cast("POLYGON")
    })
  })

  eco[[1]] <- select(eco[[1]], ECODISTRIC, geometry) |>
    dplyr::rename(ECODISTRICT = ECODISTRIC)
  eco[[2]] <- select(eco[[2]], REGION_NAM, geometry) |>
    dplyr::rename(ECOREGION = REGION_NAM)
  eco[[3]] <- select(eco[[3]], ECOPROVINC, geometry) |>
    dplyr::rename(ECOPROVINCE = ECOPROVINC)
  eco[[4]] <- select(eco[[4]], ZONE_NAME, geometry) |>
    dplyr::rename(ECOZONE = ZONE_NAME) |>
    mutate(ECOZONE = tools::toTitleCase(tolower(ECOZONE)))

  ## intersect them all, removing slivers, lines, points, etc.
  analysisZones <- eco[[4]] |>
    st_intersection(eco[[3]]) |>
    st_intersection(eco[[2]]) |>
    st_intersection(eco[[1]]) |>
    st_buffer(0)
  analysisZones <- analysisZones[which(!is.na(st_dimension(analysisZones))), ]
  rownames(analysisZones) <- 1:nrow(analysisZones)

  return(analysisZones)
}
