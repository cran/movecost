#' Accumulated cost surface and isolines around one or more origins
#'
#' Computes the anisotropic accumulated cost of movement around the input
#' origin location(s), reusing the graph stored in a
#' \code{\link{mc_surface}} object (nothing is recomputed). If more than one
#' origin is supplied, the returned surface stores, for each cell, the minimum
#' accumulated cost across origins (as in movecost <= 2.x).\cr
#'
#' Isolines (contour lines) are derived from the accumulated surface at an
#' interval set by \code{breaks}; if no value is supplied, the interval defaults
#' to one tenth of the range of the accumulated values.\cr
#'
#' For cost functions expressing cost as walking time, the accumulated values
#' are returned in the unit selected via \code{time} ("h" hours, default, or
#' "m" minutes); for the energy functions, values are expressed as accumulated
#' energy (kJ for the Pandolf, Pandolf-corrected, Van Leusen, and
#' Llobera-Sluckin functions; J/kg for the Minetti, Herzog, and Ardigo
#' functions; kcal for Hare's function); abstract-cost functions keep their
#' native unit (see the \code{cost_unit} column of
#' \code{\link{mc_cost_functions}}).
#'
#' @param surface a \code{movecost_surface} object created by
#'   \code{\link{mc_surface}}.
#' @param origin origin location(s) (sf points, SpatVector, or legacy
#'   SpatialPointsDataFrame).
#' @param time unit for time-based cost functions: "h" (hours, default) or
#'   "m" (minutes).
#' @param breaks isoline interval; if NULL (default), one tenth of the range of
#'   accumulated values.
#'
#' @return An object of class \code{movecost_accum}, a list with components
#' \itemize{
#'   \item \code{accum}: accumulated cost surface (\code{SpatRaster}), masked to
#'     the DTM (NoData cells, e.g. sea, are NA);
#'   \item \code{isolines}: contour lines (\code{sf} lines) with a \code{level}
#'     attribute;
#'   \item \code{origin}: the origin location(s) (\code{sf});
#'   \item \code{breaks}, \code{levels}, \code{time}, \code{funct},
#'     \code{funct.label}, \code{cost.unit}: metadata used by the plot method.
#' }
#' Plot it at any time, and as many times as desired, with \code{plot()} (see
#' \code{\link{plot.movecost_accum}}); no computation is re-run.
#'
#' @examples
#' dtm <- mc_volc()
#' start <- mc_volc_loc()
#'
#' surf <- mc_surface(dtm, funct = "t", move = 8)
#' acc <- mc_accum(surf, origin = start, breaks = 0.05)
#'
#' \donttest{
#' plot(acc)                 # raster + isolines
#' plot(acc, type = "contour")  # isolines only
#' }
#'
#' @seealso \code{\link{mc_surface}}, \code{\link{plot.movecost_accum}},
#'   \code{\link{mc_export}}
#' @export
mc_accum <- function(surface, origin, time = "h", breaks = NULL) {
  stopifnot(inherits(surface, "movecost_surface"))
  ori <- .points_to_cells(surface, origin, "origin")

  vals <- .accum_matrix(surface, ori$cells, mode = "out")
  vals <- if (nrow(vals) > 1L) apply(vals, 2, min) else vals[1, ]
  vals <- .convert_units(vals, surface$funct, time)
  vals[!is.finite(vals)] <- NA_real_

  accum <- terra::rast(surface$dtm)
  terra::values(accum) <- vals
  accum <- terra::mask(accum, surface$dtm)
  names(accum) <- "accumulated_cost"

  if (is.null(breaks) || breaks <= 0) breaks <- .default_breaks(vals)
  levels <- tryCatch(.contour_levels(vals, breaks), error = function(e) NULL)
  iso <- if (is.null(levels)) NULL else tryCatch(
    sf::st_as_sf(terra::as.contour(accum, levels = levels)),
    error = function(e) NULL
  )

  cf <- .cf_get(surface$funct)
  out <- list(accum = accum,
              isolines = iso,
              origin = ori$sf,
              breaks = breaks,
              levels = levels,
              time = time,
              funct = surface$funct,
              funct.label = surface$funct.label,
              cost.unit = if (cf$time_output) time else cf$unit,
              dtm = surface$dtm)
  class(out) <- "movecost_accum"
  out
}

#' @rdname movecost-methods
#' @export
print.movecost_accum <- function(x, ...) {
  rng <- range(terra::values(x$accum, mat = FALSE), na.rm = TRUE)
  cat("movecost accumulated cost surface\n")
  cat("  cost function : ", x$funct, "\n", sep = "")
  cat("  origins       : ", nrow(x$origin), "\n", sep = "")
  cat("  cost range    : ", signif(rng[1], 4), " - ", signif(rng[2], 4),
      " (", x$cost.unit, ")\n", sep = "")
  cat("  isoline interval: ", x$breaks, "\n", sep = "")
  cat("Use plot() to visualise; mc_export() to write to disk.\n")
  invisible(x)
}
