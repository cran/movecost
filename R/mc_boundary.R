#' Cost-limit boundaries around origin locations
#'
#' Delineates, around each origin location, the area reachable within a given
#' cost limit (for instance, the 30-minute walking catchment), reusing the
#' graph stored in a \code{\link{mc_surface}} object. The accumulated cost
#' surfaces for all origins are obtained with a single multi-source graph
#' query, so the function scales efficiently with the number of origins
#' (movecost <= 2.x recomputed the whole conductance matrix for each call).\cr
#'
#' Unlike movecost <= 2.x, which derived the boundary from contour lines of a
#' combined minimum-cost surface, version 3.0 polygonises each origin's
#' reachable area directly (cells with accumulated cost <= \code{limit}), which
#' yields closed polygons by construction and allows area and perimeter to be
#' computed exactly. Overlapping catchments of different origins are preserved
#' as separate (overlapping) polygons.
#'
#' @param surface a \code{movecost_surface} object created by
#'   \code{\link{mc_surface}}.
#' @param origin origin location(s) (sf points, SpatVector, or legacy
#'   SpatialPointsDataFrame).
#' @param limit cost limit defining the boundary, in the unit of the cost
#'   function (for time-based functions, in the unit selected via
#'   \code{time}); if NULL (default), one tenth of the range of accumulated
#'   values across all origins is used.
#' @param time unit for time-based cost functions: "h" (hours, default) or
#'   "m" (minutes).
#'
#' @return An object of class \code{movecost_boundary}, a list with components
#' \itemize{
#'   \item \code{boundaries}: \code{sf} polygons (one feature per origin) with
#'     \code{origin_id}, \code{limit}, \code{area} (m^2) and \code{perimeter}
#'     (m) attributes;
#'   \item \code{accum}: the minimum accumulated cost surface across origins
#'     (\code{SpatRaster}), useful for plotting;
#'   \item \code{origin}: the origin location(s) (\code{sf});
#'   \item metadata used by the plot method.
#' }
#' Plot it at any time with \code{plot()} (see
#' \code{\link{plot.movecost_boundary}}).
#'
#' @examples
#' dtm <- mc_volc()
#' destin <- mc_destin_loc()
#'
#' surf <- mc_surface(dtm, funct = "t", move = 8)
#' # 5-minute walking catchments around three locations
#' bnd <- mc_boundary(surf, origin = destin[1:3, ], limit = 5, time = "m")
#'
#' bnd$boundaries     # polygons with area and perimeter
#' plot(bnd)
#'
#' @seealso \code{\link{mc_surface}}, \code{\link{plot.movecost_boundary}}
#' @export
mc_boundary <- function(surface, origin, limit = NULL, time = "h") {
  stopifnot(inherits(surface, "movecost_surface"))
  ori <- .points_to_cells(surface, origin, "origin")

  vals <- .accum_matrix(surface, ori$cells, mode = "out")
  vals <- .convert_units(vals, surface$funct, time)

  if (is.null(limit)) {
    limit <- .default_breaks(as.vector(vals))
    message("No 'limit' supplied; using 1/10 of the accumulated cost range: ",
            limit)
  }

  crs <- sf::st_crs(terra::crs(surface$dtm))
  dtm_template <- terra::rast(surface$dtm)
  polys <- vector("list", length(ori$cells))
  for (i in seq_along(ori$cells)) {
    r <- dtm_template
    terra::values(r) <- as.numeric(vals[i, ] <= limit)
    r <- terra::mask(r, surface$dtm)
    g <- sf::st_as_sf(terra::as.polygons(r))   # aggregated by cell value
    g <- g[g[[1]] == 1, ]                      # keep the reachable area
    polys[[i]] <- sf::st_sf(origin_id = i, limit = limit,
                            geometry = sf::st_geometry(g))
  }
  boundaries <- do.call(rbind, polys)
  boundaries$area <- sf::st_area(boundaries)
  boundaries$perimeter <- sf::st_length(sf::st_cast(boundaries, "MULTILINESTRING"))

  vmin <- if (nrow(vals) > 1L) apply(vals, 2, min) else vals[1, ]
  vmin[!is.finite(vmin)] <- NA_real_
  accum <- terra::rast(surface$dtm)
  terra::values(accum) <- vmin
  accum <- terra::mask(accum, surface$dtm)
  names(accum) <- "accumulated_cost"

  cf <- .cf_get(surface$funct)
  out <- list(boundaries = boundaries,
              accum = accum,
              origin = ori$sf,
              limit = limit,
              time = time,
              funct = surface$funct,
              funct.label = surface$funct.label,
              cost.unit = if (cf$time_output) time else cf$unit,
              dtm = surface$dtm)
  class(out) <- "movecost_boundary"
  out
}

#' @rdname movecost-methods
#' @export
print.movecost_boundary <- function(x, ...) {
  cat("movecost cost-limit boundaries\n")
  cat("  cost function : ", x$funct, "\n", sep = "")
  cat("  cost limit    : ", x$limit, " (", x$cost.unit, ")\n", sep = "")
  cat("  origins       : ", nrow(x$origin), "\n", sep = "")
  print(sf::st_drop_geometry(x$boundaries))
  cat("Use plot() to visualise; mc_export() to write to disk.\n")
  invisible(x)
}
