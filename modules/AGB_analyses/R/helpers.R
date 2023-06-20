#' Derive slope of numerical vector across a time series
#'
#' @param x vector of timeseries values
#'
ts_slope <- function(x) {
  y <- which(!is.na(x))
  if (length(y) >= 2) {
    x <- unlist(x[!is.na(x)])
    if (length(unique(x)) == 1) {
      return(0)
    } else {
      return(sum((x - mean(x)) * (y - mean(y))) / sum((x - mean(x))^2))
    }
  } else {
    return(NA)
  }
}

#' number of non-NA values in a sample
#'
#' @param x vector of timeseries values
#'
ts_nsamp <- function(x) {
  x <- sum(!is.na(x))
  if (x > 1) {
    return(x)
  } else {
    return(NA)
  }
}

#' `prepZones`
#'
#' Create unique categorical zones for comparative analysis
#' (e.g. study area ZOI x age class for desired time period)
#'
#' @param zoi TODO
#' @param field TODO
#' @param ageClass TODO
#' @param file.id TODO
#' @param cropRaster TODO
#' @param overwrite logical. should the file be overwritten if it exists?
#'
#' @return TODO
#'
#' @export
prepZones <- function(zoi = st_read("outputs/WBI_studyArea.gpkg", quiet = TRUE),
                      field = "ECOZONE",
                      ageClass = rast("inputs/clean/AGB_age_mosaic_classes_t1.tif"),
                      file.id = "WBI_ecozone_t1",
                      cropRaster = NULL,
                      overwrite = TRUE) {
  if (!is.null(cropRaster)) {
    ageClass <- crop(ageClass, cropRaster)
  }

  ## QC of field to rasterize
  # if(!inherits(class(pull(zoi, field)), 'numeric')) zoi <- mutate(zoi, !!field := as.integer(factor(!!as.name(field))))

  ## define comparative zones of interest (ZOI x age class for specified time period) within study area
  rasterize(
    x = vect(zoi),
    y = ageClass,
    field = field,
    filename = paste0("outputs/ZOI_", file.id, ".tif"),
    overwrite = overwrite
  )

  # safeguard
  maxchar <- max(nchar(levels(rast(paste0("outputs/ZOI_", file.id, ".tif")))[[1]]$value))
  rcoef <- as.numeric(str_c(c(1, rep(0, maxchar)), collapse = ""))

  ## combine study area zoi and age class mosaic into unique combined raster categories (when age class from 100 - 500 [n=5])
  zoneCat <- rast(paste0("outputs/ZOI_", file.id, ".tif"))
  jointRast <- as.int(rcoef * ageClass + zoneCat)

  ## re-assign factor levels
  ftab <- data.frame(value = unique(jointRast)) |>
    rename(value = ageClass) |>
    mutate(
      ageClass = levels(ageClass)[[1]]$ageClass[match(as.integer(str_sub(value, start = 1L, end = 1L)), levels(ageClass)[[1]]$value)],
      !!field := levels(zoneCat)[[1]][, field][match(as.integer(str_sub(value, start = -maxchar)), levels(zoneCat)[[1]]$value)]
    )

  ftab <- bind_cols(ftab, zoneCat = paste0(ftab[, 3], " (", ftab[, 2], " yrs)")) |>
    relocate(zoneCat, .after = value)

  levels(jointRast) <- ftab

  writeRaster(jointRast, filename = paste0("outputs/ZOIxageClass_", file.id, ".tif"), overwrite = overwrite)

  return(invisible(NULL))
}

#' Calculate summary statistics by categorical zone of interest
#'
#' (weighted mean & sd of slopes by age class, time period & study area zones of interest)
#' - WBI and ecozones by default for every time period
#'
#' @param slopeRaster TODO
#' @param weightRaster TODO
#' @param zoneRaster TODO
#' @param cropRaster TODO
#' @param maskRaster TODO
#' @param file.id TODO: e.g. one of `"WBI_ecozones"`
#'
#' @return TODO
#'
#' @export
zoneStats <- function(slopeRaster = NULL, weightRaster = NULL, zoneRaster = NULL,
                      cropRaster = NULL, maskRaster = NULL, file.id = NULL, path = NULL) {
  names(weightRaster) <- "w"
  names(slopeRaster) <- "slope"

  ## crop to smaller (e.g. test) area, if one provided
  if (!is.null(cropRaster)) {
    slopeRaster <- crop(slopeRaster, cropRaster)
    weightRaster <- crop(weightRaster, cropRaster)
    zoneRaster <- crop(zoneRaster, cropRaster)
    if (!is.null(maskRaster)) maskRaster <- crop(maskRaster, cropRaster)
  }

  ## apply mask, if provided
  if (!is.null(maskRaster)) {
    slopeRaster <- mask(slopeRaster, maskRaster)
    weightRaster <- mask(weightRaster, maskRaster)
    zoneRaster <- mask(zoneRaster, maskRaster)
  } # 13 min

  a <- cats(zoneRaster)[[1]]
  levels(zoneRaster) <- NULL
  names(zoneRaster) <- "value"

  ## 1) count of slope observations by zone
  a <- a |>
    left_join(zonal(slopeRaster, zoneRaster, fun = "notNA"), by = "value") |>
    rename(count = slope)

  if (any(!is.na(a$count))) {
    ## 2) geometric (unweighted) mean by zone
    meanRast <- zonal(slopeRaster, zoneRaster, fun = "mean", as.raster = TRUE, na.rm = TRUE)

    a <- a |>
      left_join(zonal(meanRast, zoneRaster, fun = "min", na.rm = TRUE), by = "value") |>
      rename(mean.slope = slope)

    ## 3) sd by zone
    a <- a |>
      left_join(
        zonal(
        x = (slopeRaster - meanRast)^2,
        z = zoneRaster, fun = "sum", na.rm = TRUE
      ), by = "value") |>
      mutate(sd = sqrt(slope / count)) |>
      select(!slope)

    ## 4.1) sum of wi*xi by zone (numerator of weighted mean)
    b <- zonal(slopeRaster * weightRaster, zoneRaster, fun = "sum", na.rm = TRUE) |>
      rename(b = slope)

    ## 4.2) sum of weights by zone (denominator of weighted mean)
    w <- zonal(weightRaster, zoneRaster, fun = "sum", na.rm = TRUE)

    ## 4.3) derive weighted mean
    a <- a |>
      left_join(., b, by = "value") |>
      left_join(., w, by = "value") |>
      mutate(wtd.mean.slope = b / w, .before = w)

    rm(w)

    ## 5) derive weighted sd
    a <- a |>
      left_join(
        zonal(
          x = weightRaster * (slopeRaster - classify(zoneRaster, rcl = a[, c("value", "wtd.mean.slope")]))^2,
          z = zoneRaster, fun = "sum", na.rm = TRUE
        ) |> rename(x = w),
        by = "value"
      ) |>
      mutate(wtd.sd = sqrt(x / w)) |>
      filter(!is.na(count)) |>
      select(!c(b, w, x))

    ##  6 b) write to file
    saveRDS(a, file = file.path(path, paste0("zoneStats_summary_", file.id, ".rds")))
  }

  return(invisible(NULL))
}

#' Plot summary statistics
#'
#' @param file2plot TODO
#' @param weighted TODO
#' @param out TODO
#'
#' @return TODO
#'
#' @export
plotZoneStats <- function(file2plot = NULL, weighted = TRUE, out = 1) {
  sumtab <- readRDS(file2plot) |>
    mutate(ageClass = factor(ageClass, levels = unique(ageClass), ordered = TRUE))
  catvar <- names(sumtab)[4]

  # y-axis label corresponding to examined time period
  tref <- cbind.data.frame(
    tp = paste0("t", 1:6),
    yr1 = c(1984, 1989, 1994, 1999, 2004, 2009),
    yr2 = c(1988, 1993, 1998, 2003, 2008, 2014)
  )

  tp <- ifelse(str_detect(file2plot, "_t"),
    paste0("(", str_c(filter(tref, tp == str_sub(file2plot, start = -6L, end = -5L)) |>
                        select(yr1, yr2) |> unlist(), collapse = "-"), ")"),
    "(1984-2014)"
  )

  if (weighted) meanVar <- "wtd.mean.slope" else meanVar <- "mean.slope"
  if (weighted) sdVar <- "wtd.sd" else sdVar <- "sd"

  gp <- ggplot(sumtab, aes(x = ageClass, y = !!as.name(meanVar), group = !!as.name(catvar), color = !!as.name(catvar))) +
    xlab("stand age") +
    ylab(paste0("n-weighted mean AGB trend ", tp)) +
    scale_x_discrete(breaks = sumtab$ageClass, labels = sumtab$ageClass) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_line(linetype = 3) +
    geom_errorbar(aes(ymin = !!as.name(meanVar) - (!!as.name(sdVar) / sqrt(count)), ymax = !!as.name(meanVar) + (!!as.name(sdVar) / sqrt(count))),
      width = .2
    )

  if (out == 1) plot(gp) else return(gp)
}

plotZoneStatsIntervals <- function(files2plot = NULL, weighted = TRUE, xVar = "ageClass",
                                   groupVar = "tp", catVar = NULL, ptype = 1, plotResult = TRUE) {
  ## compile summary tables
  sumtab <- do.call(rbind, lapply(files2plot, function(x) {
    y <- readRDS(x) |>
      mutate(ageClass = factor(ageClass, levels = unique(ageClass), ordered = TRUE))
    y$tp <- factor(rep(str_sub(x, start = -6L, end = -5L), nrow(y)),
      levels = paste0("t", 1:6), ordered = TRUE
    )
    return(y |> relocate(tp, .before = value))
  }))

  ## validate arguments
  if (weighted) meanVar <- "wtd.mean.slope" else meanVar <- "mean.slope"
  if (weighted) sdVar <- "wtd.sd" else sdVar <- "sd"

  if (is.null(catVar)) catVar <- names(sumtab)[5]
  if (!is.null(groupVar)) groupVar <- match.arg(groupVar, names(sumtab))
  xVar <- match.arg(xVar, names(sumtab))
  if (xVar == "ageClass") xlabel <- "stand age"
  if (xVar == "tp") xlabel <- "5-year time period (1984-2014)"

  if (ptype == 1) {
    gp <- ggplot(sumtab, aes(x = !!as.name(xVar), y = !!as.name(meanVar), group = !!as.name(groupVar), color = !!as.name(groupVar))) +
      labs(
        x = xlabel,
        y = "n-weighted mean AGB trend"
      ) +
      geom_hline(yintercept = 0, linetype = "dotted") +
      geom_line() +
      facet_grid(~ .data[[catVar]])
  } else {
    if (ptype == 2) {
      gp <- lapply(unique(sumtab[, catVar]), function(cat) {
        x <- filter(sumtab, !!as.name(catVar) == cat)

        return(ggplot(x, aes(x = !!as.name(xVar), y = !!as.name(meanVar), group = !!as.name(groupVar), color = !!as.name(groupVar))) +
          geom_line() +
          xlab(xlabel) +
          ylab("n-weighted mean AGB trend") +
          geom_hline(yintercept = 0, linetype = "dotted") +
          geom_line() +
          ggtitle(cat))
      })
      names(gp) <- unique(sumtab[, catVar])
    }

    if (!plotResult) {
      return(gp)
    }
    if (ptype == 1) print(gp) else lapply(gp, print)
  }
}
