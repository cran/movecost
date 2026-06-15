#' Internal utilities for movecost 3.0
#'
#' These helpers are not exported. They centralise input coercion,
#' validation, unit handling, and time formatting so that every
#' user-facing function behaves consistently.
#'
#' @name movecost-internals
#' @keywords internal
NULL

#' Coerce any supported vector-input type to 'sf'
#'
#' Accepts 'sf', 'sfc', 'SpatVector', and legacy 'sp' objects
#' (SpatialPointsDataFrame, SpatialPoints, SpatialLinesDataFrame,
#' SpatialPolygonsDataFrame, ...), so that scripts written for
#' movecost <= 2.x keep working when fed to the new functions.
#'
#' @param x vector data in any supported class.
#' @param arg_name name of the argument (used in error messages).
#' @return an 'sf' object.
#' @keywords internal
#' @noRd
.as_sf <- function(x, arg_name = "input") {
  if (inherits(x, "sf")) return(x)
  if (inherits(x, "sfc")) return(sf::st_as_sf(x))
  if (inherits(x, "SpatVector")) return(sf::st_as_sf(x))
  if (inherits(x, c("SpatialPointsDataFrame", "SpatialPoints",
                    "SpatialLinesDataFrame", "SpatialLines",
                    "SpatialPolygonsDataFrame", "SpatialPolygons"))) {
    return(sf::st_as_sf(x))
  }
  stop("'", arg_name, "' must be an sf, sfc, SpatVector, or sp object; got: ",
       paste(class(x), collapse = "/"), call. = FALSE)
}

#' Coerce any supported raster-input type to 'SpatRaster'
#'
#' Accepts 'SpatRaster', legacy 'RasterLayer', or a file path.
#'
#' @param x raster data in any supported class.
#' @return a 'SpatRaster'.
#' @keywords internal
#' @noRd
.as_spatraster <- function(x) {
  if (inherits(x, "SpatRaster")) return(x)
  if (inherits(x, "RasterLayer")) return(terra::rast(x))
  if (is.character(x) && length(x) == 1L) return(terra::rast(x))
  stop("'dtm' must be a SpatRaster, RasterLayer, or a file path; got: ",
       paste(class(x), collapse = "/"), call. = FALSE)
}

#' Validate that a raster uses a projected coordinate system
#' @keywords internal
#' @noRd
.assert_projected <- function(r) {
  ll <- terra::is.lonlat(r)
  if (is.na(ll)) {
    stop("The input DTM has no coordinate reference system.\n",
         "Assign a projected CRS first, e.g.: terra::crs(dtm) <- 'EPSG:32633'.",
         call. = FALSE)
  }
  if (ll) {
    stop("The input DTM uses a geographic (lon/lat) coordinate system.\n",
         "movecost requires a projected coordinate system (planar units in metres).\n",
         "Reproject first, e.g.: terra::project(dtm, 'EPSG:32633').", call. = FALSE)
  }
  invisible(TRUE)
}

#' Align a vector layer to the DTM's CRS
#'
#' The same CRS can be encoded differently by different file formats or
#' libraries (WKT variants, legacy proj4 strings), in which case a strict
#' equality test produces false alarms. This helper therefore TRANSFORMS the
#' vector layer to the DTM's CRS whenever the two CRS objects are not equal:
#' if they are equivalent, the transformation is an identity; if they truly
#' differ, the layer is properly reprojected (with a message). A layer with
#' no CRS at all is assumed to be in the DTM's CRS, with a warning.
#'
#' @param r the reference SpatRaster (the DTM).
#' @param v an sf object.
#' @param arg_name name of the argument (used in messages).
#' @return the sf object, in the DTM's CRS.
#' @keywords internal
#' @noRd
.align_crs <- function(r, v, arg_name = "input") {
  crs_r <- sf::st_crs(terra::crs(r))
  crs_v <- sf::st_crs(v)
  if (is.na(crs_r)) return(v)   # DTM CRS unknown: nothing to align to
  if (is.na(crs_v)) {
    warning("'", arg_name, "' has no coordinate reference system; ",
            "assuming the DTM's CRS.", call. = FALSE)
    sf::st_crs(v) <- crs_r
    return(v)
  }
  if (crs_r != crs_v) {
    v <- sf::st_transform(v, crs_r)
    message("Note: '", arg_name,
            "' was aligned (transformed) to the DTM's CRS.")
  }
  v
}

#' Convert seconds to the user-selected time unit
#' @keywords internal
#' @noRd
.time_convert <- function(x, time = c("h", "m")) {
  time <- match.arg(time)
  if (time == "h") x / 3600 else x / 60
}

#' Format a cost expressed in decimal hours or minutes as HH:MM:SS
#'
#' Replaces the dependency on 'chron' used in movecost <= 2.x.
#'
#' @param x numeric vector of decimal time values.
#' @param time unit of 'x': "h" (hours) or "m" (minutes).
#' @return character vector in "HH:MM:SS" format ('NA' kept as NA).
#' @keywords internal
#' @noRd
.fmt_hms <- function(x, time = c("h", "m")) {
  time <- match.arg(time)
  secs <- if (time == "h") x * 3600 else x * 60
  out <- rep(NA_character_, length(x))
  ok <- is.finite(secs)
  s <- round(secs[ok])
  out[ok] <- sprintf("%02d:%02d:%02d", s %/% 3600L, (s %% 3600L) %/% 60L, s %% 60L)
  out
}

#' Default contour interval: one tenth of the finite value range
#'
#' signif() (not round()) so that small cost ranges never collapse to a
#' zero interval.
#'
#' @keywords internal
#' @noRd
.default_breaks <- function(vals) {
  v <- vals[is.finite(vals)]
  signif((max(v) - min(v)) / 10, 2)
}

#' Sequence of contour levels over the finite value range
#' @keywords internal
#' @noRd
.contour_levels <- function(vals, breaks) {
  v <- vals[is.finite(vals)]
  seq(min(v), max(v), by = breaks)
}

#' Cells of a SpatRaster touched by an sf geometry
#'
#' Used for barriers and for path-density counts.
#'
#' @keywords internal
#' @noRd
.cells_touched <- function(r, geom_sf) {
  v <- terra::vect(geom_sf)
  ex <- terra::extract(r, v, cells = TRUE, touches = TRUE)
  unique(ex$cell[!is.na(ex$cell)])
}

#' Null-coalescing helper
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a
