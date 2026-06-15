#' Cost allocation (Thiessen-like polygons based on movement cost)
#'
#' Partitions the landscape into zones of pertinence, assigning each DTM cell
#' to the origin location from which it can be reached at the smallest
#' accumulated cost. The result is the cost-based analogue of Voronoi/Thiessen
#' polygons.\cr
#'
#' All the accumulated cost surfaces are obtained with a single multi-source
#' graph query on the \code{\link{mc_surface}} object: with n origins,
#' movecost <= 2.x performed n full conductance-matrix constructions plus n
#' accumulation runs, whereas this implementation performs none of the former
#' and one batched run of the latter, making allocation over many sites orders
#' of magnitude faster.
#'
#' @param surface a \code{movecost_surface} object created by
#'   \code{\link{mc_surface}}.
#' @param origin origin location(s) (sf points, SpatVector, or legacy
#'   SpatialPointsDataFrame); at least two.
#' @param time unit for time-based cost functions: "h" (hours, default) or
#'   "m" (minutes).
#' @param breaks optional isoline interval for the minimum-cost surface; if
#'   NULL (default), one tenth of its range. Isolines are returned for optional
#'   display in the plot method.
#'
#' @return An object of class \code{movecost_alloc}, a list with components
#' \itemize{
#'   \item \code{alloc}: allocation zones (\code{SpatRaster}; cell value = index
#'     of the prevailing origin);
#'   \item \code{zones}: the same zones as polygons (\code{sf}) with an
#'     \code{origin_id} attribute;
#'   \item \code{accum.min}: minimum accumulated cost across origins
#'     (\code{SpatRaster});
#'   \item \code{isolines}: isolines of \code{accum.min} (\code{sf});
#'   \item \code{origin}: the origin locations (\code{sf});
#'   \item metadata used by the plot method.
#' }
#' Plot it at any time with \code{plot()} (see
#' \code{\link{plot.movecost_alloc}}).
#'
#' @examples
#' dtm <- mc_volc()
#' destin <- mc_destin_loc()
#'
#' surf <- mc_surface(dtm, funct = "t", move = 8)
#' al <- mc_alloc(surf, origin = destin)
#'
#' plot(al)                    # allocation zones
#' plot(al, isolines = TRUE)   # with cost isolines overlaid
#'
#' @seealso \code{\link{mc_surface}}, \code{\link{plot.movecost_alloc}}
#' @export
mc_alloc <- function(surface, origin, time = "h", breaks = NULL) {
  stopifnot(inherits(surface, "movecost_surface"))
  ori <- .points_to_cells(surface, origin, "origin")
  if (length(ori$cells) < 2L) {
    stop("Cost allocation requires at least two origin locations.", call. = FALSE)
  }

  vals <- .accum_matrix(surface, ori$cells, mode = "out")
  vals <- .convert_units(vals, surface$funct, time)

  finite_any <- apply(vals, 2, function(col) any(is.finite(col)))
  zone <- rep(NA_real_, ncol(vals))
  zone[finite_any] <- apply(vals[, finite_any, drop = FALSE], 2, which.min)
  vmin <- rep(NA_real_, ncol(vals))
  vmin[finite_any] <- apply(vals[, finite_any, drop = FALSE], 2, min)

  alloc <- terra::rast(surface$dtm)
  terra::values(alloc) <- zone
  alloc <- terra::mask(alloc, surface$dtm)
  names(alloc) <- "allocation"

  accum.min <- terra::rast(surface$dtm)
  terra::values(accum.min) <- vmin
  accum.min <- terra::mask(accum.min, surface$dtm)
  names(accum.min) <- "min_accumulated_cost"

  zones <- sf::st_as_sf(terra::as.polygons(alloc))   # aggregated by zone id
  names(zones)[1] <- "origin_id"

  if (is.null(breaks)) breaks <- .default_breaks(vmin)
  iso <- tryCatch(
    sf::st_as_sf(terra::as.contour(accum.min,
                                   levels = .contour_levels(vmin, breaks))),
    error = function(e) NULL
  )

  cf <- .cf_get(surface$funct)
  out <- list(alloc = alloc,
              zones = zones,
              accum.min = accum.min,
              isolines = iso,
              origin = ori$sf,
              breaks = breaks,
              time = time,
              funct = surface$funct,
              funct.label = surface$funct.label,
              cost.unit = if (cf$time_output) time else cf$unit,
              dtm = surface$dtm)
  class(out) <- "movecost_alloc"
  out
}

#' @rdname movecost-methods
#' @export
print.movecost_alloc <- function(x, ...) {
  cat("movecost cost allocation\n")
  cat("  cost function : ", x$funct, "\n", sep = "")
  cat("  origins       : ", nrow(x$origin), "\n", sep = "")
  cat("  zones         : ", nrow(x$zones), " polygons\n", sep = "")
  cat("Use plot() to visualise; mc_export() to write to disk.\n")
  invisible(x)
}
