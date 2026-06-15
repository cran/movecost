#' Least-cost paths from an origin to one or more destinations
#'
#' Computes least-cost path(s) (LCPs) from an origin to one or more destination
#' locations, reusing the graph stored in a \code{\link{mc_surface}} object.
#' Because terrain costs are anisotropic, the path back from a destination to
#' the origin can differ from the outward path; set \code{return.base = TRUE}
#' to also compute the return paths.\cr
#'
#' Each returned path carries its accumulated \code{cost} (in the unit of the
#' surface's cost function; for time-based functions, in the unit selected via
#' \code{time}) and its \code{length}. For time-based functions a
#' \code{cost_hms} column additionally expresses the cost in sexagesimal format
#' (HH:MM:SS). The same attributes are appended to the returned copy of the
#' destination locations, mirroring the \code{dest.loc.w.cost} output of
#' movecost <= 2.x.
#'
#' @param surface a \code{movecost_surface} object created by
#'   \code{\link{mc_surface}}.
#' @param origin single origin location (sf points, SpatVector, or legacy
#'   SpatialPointsDataFrame); if more than one feature is supplied, the first
#'   is used.
#' @param destin destination location(s) (same accepted classes).
#' @param return.base TRUE or FALSE (default): also compute the least-cost
#'   paths from each destination back to the origin.
#' @param time unit for time-based cost functions: "h" (hours, default) or
#'   "m" (minutes).
#'
#' @return An object of class \code{movecost_paths}, a list with components
#' \itemize{
#'   \item \code{paths}: \code{sf} lines with \code{cost}, \code{length} and
#'     (for time-based functions) \code{cost_hms} attributes;
#'   \item \code{paths.back}: \code{sf} lines of the return paths (NULL unless
#'     \code{return.base = TRUE});
#'   \item \code{destin}: copy of the destination locations with \code{cost}
#'     (and \code{cost_hms}) appended;
#'   \item \code{origin}: the origin (\code{sf});
#'   \item metadata used by the plot method.
#' }
#' Plot it at any time with \code{plot()} (see
#' \code{\link{plot.movecost_paths}}); no computation is re-run.
#'
#' @examples
#' dtm <- mc_volc()
#' start <- mc_volc_loc()
#' destin <- mc_destin_loc()
#'
#' surf <- mc_surface(dtm, funct = "t", move = 8)
#' lcp <- mc_paths(surf, origin = start, destin = destin)
#'
#' lcp$paths          # sf object: costs, lengths, geometries
#' plot(lcp)
#'
#' # paths back to the origin (anisotropy makes them differ):
#' lcp2 <- mc_paths(surf, start, destin[2, ], return.base = TRUE)
#' plot(lcp2)
#'
#' @seealso \code{\link{mc_surface}}, \code{\link{plot.movecost_paths}},
#'   \code{\link{mc_rank}}, \code{\link{mc_network}}
#' @export
mc_paths <- function(surface, origin, destin, return.base = FALSE, time = "h") {
  stopifnot(inherits(surface, "movecost_surface"))
  ori <- .points_to_cells(surface, origin, "origin")
  dst <- .points_to_cells(surface, destin, "destin")
  if (length(ori$cells) > 1L) {
    warning("Multiple origin features supplied; using the first one.", call. = FALSE)
    ori$cells <- ori$cells[1]
    ori$sf <- ori$sf[1, ]
  }

  paths <- .lcp_sf(surface, ori$cells, dst$cells, time)

  destin_out <- dst$sf
  destin_out$cost <- paths$cost
  cf <- .cf_get(surface$funct)
  if (cf$time_output) {
    destin_out$cost_hms <- .fmt_hms(paths$cost, time)
    paths$cost_hms <- destin_out$cost_hms
  }

  paths.back <- NULL
  if (isTRUE(return.base)) {
    backs <- lapply(dst$cells, function(dc) .lcp_sf(surface, dc, ori$cells, time))
    paths.back <- do.call(rbind, backs)
    if (cf$time_output) paths.back$cost_hms <- .fmt_hms(paths.back$cost, time)
  }

  out <- list(paths = paths,
              paths.back = paths.back,
              origin = ori$sf,
              destin = destin_out,
              time = time,
              funct = surface$funct,
              funct.label = surface$funct.label,
              cost.unit = if (cf$time_output) time else cf$unit,
              dtm = surface$dtm,
              barrier = surface$barrier)
  class(out) <- "movecost_paths"
  out
}

#' Internal: least-cost paths between cells as an sf object
#'
#' Runs Dijkstra once from the source cell, extracts the vertex paths to all
#' target cells, converts them to LINESTRINGs, and attaches accumulated cost
#' (converted to user units) and geometric length.
#'
#' @keywords internal
#' @noRd
.lcp_sf <- function(surface, from_cell, to_cells, time = "h") {
  sp <- igraph::shortest_paths(surface$graph,
                               from = from_cell,
                               to = to_cells,
                               mode = "out",
                               output = "vpath",
                               algorithm = "dijkstra")
  costs_native <- .accum_matrix(surface, from_cell, mode = "out")[1, to_cells]
  costs <- .convert_units(costs_native, surface$funct, time)

  crs <- sf::st_crs(terra::crs(surface$dtm))
  geoms <- lapply(seq_along(sp$vpath), function(i) {
    v <- as.integer(sp$vpath[[i]])
    if (length(v) < 2L) {
      # unreachable destination, or destination == origin: empty linestring
      return(sf::st_linestring())
    }
    sf::st_linestring(terra::xyFromCell(surface$dtm, v))
  })
  out <- sf::st_sf(id = seq_along(geoms),
                   cost = as.numeric(costs),
                   geometry = sf::st_sfc(geoms, crs = crs))
  out$length <- sf::st_length(out)
  if (any(sf::st_is_empty(out) & is.infinite(out$cost))) {
    warning("One or more destinations are unreachable from the origin ",
            "(e.g. separated by a barrier or NoData cells); ",
            "their paths are empty and their cost is Inf.", call. = FALSE)
  }
  out
}

#' @rdname movecost-methods
#' @export
print.movecost_paths <- function(x, ...) {
  cat("movecost least-cost paths\n")
  cat("  cost function : ", x$funct, "\n", sep = "")
  cat("  paths         : ", nrow(x$paths),
      if (!is.null(x$paths.back)) paste0(" (+", nrow(x$paths.back), " back to origin)"),
      "\n", sep = "")
  cat("  cost unit     : ", x$cost.unit, "\n", sep = "")
  print(sf::st_drop_geometry(x$paths))
  cat("Use plot() to visualise; mc_export() to write to disk.\n")
  invisible(x)
}
