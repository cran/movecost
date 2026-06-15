#' Ranked (sub-optimal) least-cost paths between an origin and a destination
#'
#' Computes the optimal least-cost path and up to \code{k - 1} progressively
#' sub-optimal alternatives between two locations. Alternatives are obtained by
#' iterative penalisation: after each path is found, the conductance of the
#' graph edges incident to its cells is multiplied by \code{penalty} (0.01 by
#' default), and Dijkstra's algorithm is run again. This discourages, without
#' strictly forbidding, the reuse of earlier routes, and therefore returns
#' spatially distinct alternatives rather than trivial one-cell deviations
#' (which is what an exact k-shortest-path algorithm would return on a raster
#' graph). The procedure is the same in spirit as movecost <= 2.x's
#' \code{moverank()}, but operates directly on the graph edge weights, so the
#' conductance matrix is never rebuilt and any number of ranks can be
#' requested (movecost <= 2.x was limited to 6).\cr
#'
#' Each returned path carries its \code{rank} (1 = optimal), its accumulated
#' \code{cost} computed on the \emph{original} (unpenalised) graph - so costs
#' of different ranks are directly comparable - and its \code{length}.
#'
#' @param surface a \code{movecost_surface} object created by
#'   \code{\link{mc_surface}}.
#' @param origin single origin location (sf points, SpatVector, or legacy
#'   SpatialPointsDataFrame).
#' @param destin single destination location (same accepted classes).
#' @param k number of paths to return, optimal one included (default 3).
#' @param penalty conductance multiplier applied to edges incident to already
#'   found paths (default 0.01; smaller values push alternatives further away
#'   from earlier routes).
#' @param time unit for time-based cost functions: "h" (hours, default) or
#'   "m" (minutes).
#'
#' @return An object of class \code{movecost_rank}, a list with components
#' \itemize{
#'   \item \code{paths}: \code{sf} lines with \code{rank}, \code{cost},
#'     (for time-based functions) \code{cost_hms}, and \code{length}
#'     attributes;
#'   \item \code{corridor}: the least-cost corridor (\code{SpatRaster}) between
#'     the two locations, stored so that \code{plot(x, type = "corridor")} can
#'     draw the ranked paths over it without further computation;
#'   \item \code{origin}, \code{destin}: the input locations (\code{sf});
#'   \item metadata used by the plot method.
#' }
#' Plot it at any time with \code{plot()}: \code{plot(x)} draws the paths over
#' the terrain, \code{plot(x, type = "corridor")} draws them over the
#' least-cost corridor, and \code{plot(x, type = "chart")} draws a bubble chart
#' of length and cost by rank (see \code{\link{plot.movecost_rank}}).
#'
#' @examples
#' dtm <- mc_volc()
#' start <- mc_volc_loc()
#' destin <- mc_destin_loc()
#'
#' surf <- mc_surface(dtm, funct = "t", move = 8)
#' rk <- mc_rank(surf, origin = start, destin = destin[2, ], k = 3)
#'
#' rk$paths        # ranks, costs, lengths
#' plot(rk)
#'
#' @seealso \code{\link{mc_surface}}, \code{\link{plot.movecost_rank}},
#'   \code{\link{mc_paths}}
#' @export
mc_rank <- function(surface, origin, destin, k = 3, penalty = 0.01, time = "h") {
  stopifnot(inherits(surface, "movecost_surface"))
  if (k < 1) stop("'k' must be at least 1.", call. = FALSE)
  if (penalty <= 0 || penalty >= 1) {
    stop("'penalty' must be strictly between 0 and 1.", call. = FALSE)
  }
  ori <- .points_to_cells(surface, origin, "origin")
  dst <- .points_to_cells(surface, destin, "destin")
  o_cell <- ori$cells[1]
  d_cell <- dst$cells[1]

  g_work <- surface$graph       # working copy: igraph objects copy on modify
  crs <- sf::st_crs(terra::crs(surface$dtm))
  cf <- .cf_get(surface$funct)

  geoms <- vector("list", k)
  costs <- numeric(k)
  for (r in seq_len(k)) {
    sp <- igraph::shortest_paths(g_work, from = o_cell, to = d_cell,
                                 mode = "out", output = "vpath",
                                 algorithm = "dijkstra")
    v <- as.integer(sp$vpath[[1]])
    if (length(v) < 2L) {
      warning("Only ", r - 1L, " path(s) could be found; ",
              "the destination became unreachable at rank ", r, ".",
              call. = FALSE)
      geoms <- geoms[seq_len(r - 1L)]
      costs <- costs[seq_len(r - 1L)]
      break
    }
    geoms[[r]] <- sf::st_linestring(terra::xyFromCell(surface$dtm, v))
    # cost evaluated on the ORIGINAL graph, so ranks are comparable
    ep <- igraph::E(surface$graph, path = v)
    costs[r] <- sum(igraph::edge_attr(surface$graph, "weight", ep))
    # penalise edges incident to the path's vertices for the next iteration
    if (r < k) {
      eids <- igraph::incident_edges(g_work, v, mode = "all")
      eids <- unique(as.integer(unlist(eids)))
      w <- igraph::edge_attr(g_work, "weight", eids)
      g_work <- igraph::set_edge_attr(g_work, "weight", index = eids,
                                      value = w / penalty)
    }
  }

  costs <- .convert_units(costs, surface$funct, time)
  paths <- sf::st_sf(rank = seq_along(geoms),
                     cost = costs,
                     geometry = sf::st_sfc(geoms, crs = crs))
  if (cf$time_output) paths$cost_hms <- .fmt_hms(paths$cost, time)
  paths$length <- sf::st_length(paths)

  # least-cost corridor (classic "reach" form) between the two locations, so
  # that plot(type = "corridor") can overlay the ranked paths on it without
  # the user having to compute it separately (two distance calls; inexpensive)
  cv <- .accum_matrix(surface, o_cell, mode = "out")[1, ] +
        .accum_matrix(surface, d_cell, mode = "out")[1, ]
  cv <- .convert_units(cv, surface$funct, time)
  cv[!is.finite(cv)] <- NA_real_
  corridor <- terra::rast(surface$dtm)
  terra::values(corridor) <- cv
  corridor <- terra::mask(corridor, surface$dtm)
  names(corridor) <- "corridor_cost"

  out <- list(paths = paths,
              corridor = corridor,
              origin = ori$sf[1, ],
              destin = dst$sf[1, ],
              k = k,
              penalty = penalty,
              time = time,
              funct = surface$funct,
              funct.label = surface$funct.label,
              cost.unit = if (cf$time_output) time else cf$unit,
              dtm = surface$dtm)
  class(out) <- "movecost_rank"
  out
}

#' @rdname movecost-methods
#' @export
print.movecost_rank <- function(x, ...) {
  cat("movecost ranked least-cost paths\n")
  cat("  cost function : ", x$funct, "\n", sep = "")
  cat("  paths returned: ", nrow(x$paths), " (penalty = ", x$penalty, ")\n", sep = "")
  print(sf::st_drop_geometry(x$paths))
  cat("Use plot() to visualise; mc_export() to write to disk.\n")
  invisible(x)
}
