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
                 overwrite = TRUE) %>%
        st_transform(targetCRS) %>%
        st_cast("POLYGON")
    })
  })

  eco[[1]] <- select(eco[[1]], ECODISTRIC, geometry) %>%
    mutate(ECODISTRICT = ECODISTRIC, ECODISTRIC = NULL, .before = "geometry")
  eco[[2]] <- select(eco[[2]], REGION_NAM, geometry) %>%
    mutate(ECOREGION = REGION_NAM, REGION_NAM = NULL, .before = "geometry")
  eco[[3]] <- select(eco[[3]], ECOPROVINC, geometry) %>%
    mutate(ECOPROVINCE = ECOPROVINC, ECOPROVINC = NULL, .before = "geometry")
  eco[[4]] <- select(eco[[4]], ZONE_NAME, geometry) %>%
    mutate(ECOZONE = tools::toTitleCase(tolower(ZONE_NAME)), ZONE_NAME = NULL, .before = "geometry")

  ## intersect them all, removing slivers, lines, points, etc.
  analysisZones <- eco[[4]] %>%
    st_intersection(eco[[3]]) %>%
    st_intersection(eco[[2]]) %>%
    st_intersection(eco[[1]]) %>%
    st_buffer(0)
  analysisZones <- analysisZones[, -which(grepl("[.]before", colnames(analysisZones)))]
  analysisZones <- analysisZones[which(!is.na(st_dimension(analysisZones))), ]
  rownames(analysisZones) <- 1:nrow(analysisZones)

  return(analysisZones)
}
