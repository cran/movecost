#' Least-cost network between multiple locations
#'
#' Computes least-cost paths between multiple locations, either between all
#' unique pairs (\code{type = "allpairs"}) or from each location to its closest
#' (cost-wise) neighbour (\code{type = "neigh"}), together with the full
#' origin-to-origin cost matrix and, optionally, a path-density raster.\cr
#'
#' The graph stored in the \code{\link{mc_surface}} object is reused
#' throughout: the cost matrix requires a single multi-source query and the
#' paths one Dijkstra run per location, whereas movecost <= 2.x rebuilt the
#' conductance matrix before any of this could start. For n locations this
#' reduces the dominant cost from O(n) full matrix constructions to zero.\cr
#'
#' The path-density raster (\code{density = TRUE}) counts, for each cell, how
#' many least-cost paths traverse it; the percentage version expresses the
#' count against the total number of paths. Density is computed exactly, from
#' the graph vertices of each path (no rasterisation approximation).
#'
#' @param surface a \code{movecost_surface} object created by
#'   \code{\link{mc_surface}}.
#' @param nodes locations to be connected (sf points, SpatVector, or legacy
#'   SpatialPointsDataFrame); at least two.
#' @param type "allpairs" (default) for LCPs between all unique pairs of
#'   locations (in both directions, since costs are anisotropic), or "neigh"
#'   for LCPs from each location to its cost-wise nearest neighbour.
#' @param density TRUE or FALSE (default): also compute the path-density
#'   rasters.
#' @param time unit for time-based cost functions: "h" (hours, default) or
#'   "m" (minutes).
#'
#' @return An object of class \code{movecost_network}, a list with components
#' \itemize{
#'   \item \code{paths}: \code{sf} lines with \code{from}, \code{to},
#'     \code{cost}, and \code{length} attributes;
#'   \item \code{cost.matrix}: full asymmetric cost matrix between locations
#'     (in the unit of the cost function; for time-based functions, in the
#'     unit selected via \code{time});
#'   \item \code{density.count}, \code{density.perc}: path-density rasters
#'     (\code{SpatRaster}; NULL unless \code{density = TRUE});
#'   \item \code{nodes}: the input locations (\code{sf}), with a \code{node}
#'     column holding the index (1..n) that matches the row and column names
#'     of \code{cost.matrix} and the node labels drawn by \code{plot()};
#'   \item metadata used by the plot method.
#' }
#' Plot it at any time with \code{plot()} (see
#' \code{\link{plot.movecost_network}}).
#'
#' @examples
#' dtm <- mc_volc()
#' destin <- mc_destin_loc()
#'
#' surf <- mc_surface(dtm, funct = "t", move = 8)
#' nw <- mc_network(surf, nodes = destin, type = "allpairs")
#'
#' nw$cost.matrix              # asymmetric origin-to-origin costs
#' plot(nw)
#'
#' # network connecting each site to its cost-wise nearest neighbour
#' nw2 <- mc_network(surf, nodes = destin, type = "neigh")
#' plot(nw2)
#'
#' @seealso \code{\link{mc_surface}}, \code{\link{plot.movecost_network}}
#' @export
mc_network <- function(surface, nodes, type = c("allpairs", "neigh"),
                       density = FALSE, time = "h") {
  stopifnot(inherits(surface, "movecost_surface"))
  type <- match.arg(type)
  nds <- .points_to_cells(surface, nodes, "nodes")
  n <- length(nds$cells)
  if (n < 2L) stop("A network requires at least two locations.", call. = FALSE)

  # full asymmetric cost matrix from one batched query
  cm <- igraph::distances(surface$graph, v = nds$cells, to = nds$cells,
                          mode = "out", algorithm = "dijkstra")
  cm <- .convert_units(cm, surface$funct, time)
  dimnames(cm) <- list(seq_len(n), seq_len(n))

  # which pairs to trace
  if (type == "allpairs") {
    pairs <- which(matrix(TRUE, n, n) & !diag(n), arr.ind = TRUE)
  } else {
    cm_na <- cm; diag(cm_na) <- NA
    nearest <- apply(cm_na, 1, which.min)
    pairs <- cbind(seq_len(n), nearest)
  }

  crs <- sf::st_crs(terra::crs(surface$dtm))
  paths_list <- vector("list", nrow(pairs))
  vcells_all <- vector("list", nrow(pairs))
  for (i in seq_len(n)) {
    sel <- which(pairs[, 1] == i)
    if (!length(sel)) next
    sp <- igraph::shortest_paths(surface$graph,
                                 from = nds$cells[i],
                                 to = nds$cells[pairs[sel, 2]],
                                 mode = "out", output = "vpath",
                                 algorithm = "dijkstra")
    for (k in seq_along(sel)) {
      v <- as.integer(sp$vpath[[k]])
      vcells_all[[sel[k]]] <- v
      geom <- if (length(v) < 2L) sf::st_linestring()
              else sf::st_linestring(terra::xyFromCell(surface$dtm, v))
      paths_list[[sel[k]]] <- sf::st_sf(
        from = pairs[sel[k], 1], to = pairs[sel[k], 2],
        cost = cm[pairs[sel[k], 1], pairs[sel[k], 2]],
        geometry = sf::st_sfc(geom, crs = crs))
    }
  }
  paths <- do.call(rbind, paths_list)
  paths$length <- sf::st_length(paths)

  density.count <- density.perc <- NULL
  if (isTRUE(density)) {
    counts <- table(unlist(vcells_all))
    dvals <- rep(NA_real_, terra::ncell(surface$dtm))
    dvals[as.integer(names(counts))] <- as.numeric(counts)
    density.count <- terra::rast(surface$dtm)
    terra::values(density.count) <- dvals
    density.count <- terra::mask(density.count, surface$dtm)
    names(density.count) <- "path_count"
    density.perc <- density.count / nrow(paths) * 100
    names(density.perc) <- "path_percentage"
  }

  # number the nodes so their index is explicit in the data too, matching the
  # row/column names of cost.matrix (and the labels drawn by plot())
  nodes_out <- nds$sf
  nodes_out$node <- seq_len(n)

  cf <- .cf_get(surface$funct)
  out <- list(paths = paths,
              cost.matrix = cm,
              density.count = density.count,
              density.perc = density.perc,
              nodes = nodes_out,
              type = type,
              time = time,
              funct = surface$funct,
              funct.label = surface$funct.label,
              cost.unit = if (cf$time_output) time else cf$unit,
              dtm = surface$dtm)
  class(out) <- "movecost_network"
  out
}

#' @rdname movecost-methods
#' @export
print.movecost_network <- function(x, ...) {
  cat("movecost least-cost network (type: ", x$type, ")\n", sep = "")
  cat("  cost function : ", x$funct, "\n", sep = "")
  cat("  locations     : ", nrow(x$nodes), "; paths: ", nrow(x$paths), "\n", sep = "")
  cat("  cost unit     : ", x$cost.unit, "\n", sep = "")
  cat("Cost matrix (asymmetric):\n")
  print(round(x$cost.matrix, 3))
  cat("Use plot() to visualise; mc_export() to write to disk.\n")
  invisible(x)
}
