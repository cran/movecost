#' Compare least-cost paths produced by different cost functions
#'
#' Computes, for the same origin and destination(s), the least-cost paths
#' produced by two or more cost functions, so that route choice and cost can be
#' compared across cost models. One cost surface is built per function (this is
#' unavoidable, since each function defines its own graph weights), but each is
#' built once and all the remaining computation reuses it.\cr
#'
#' Note that costs produced by functions of different families (walking time,
#' energy, abstract) are expressed in different units and are comparable in
#' terms of route geometry, not of cost magnitude; the cost table returned
#' reports the unit alongside each value.
#'
#' @param dtm Digital Terrain Model (SpatRaster, RasterLayer, or file path),
#'   as in \code{\link{mc_surface}}.
#' @param origin single origin location (sf points, SpatVector, or legacy
#'   SpatialPointsDataFrame).
#' @param destin destination location(s) (same accepted classes).
#' @param functs character vector of cost-function codes to compare (see
#'   \code{\link{mc_cost_functions}}), e.g. \code{c("t", "ks", "hrz")}.
#' @param move number of directions in which cells are connected (4, 8, or 16;
#'   default 16).
#' @param return.base TRUE or FALSE (default): also compute the paths back to
#'   the origin for each cost function.
#' @param time unit for time-based cost functions: "h" (hours, default) or
#'   "m" (minutes).
#' @param ... further parameters passed to \code{\link{mc_surface}}
#'   (\code{N}, \code{W}, \code{L}, \code{V}, \code{sl.crit}, \code{barrier},
#'   \code{field}, \code{cogn.slp}, \code{topo.dist}).
#'
#' @return An object of class \code{movecost_comp}, a list with components
#' \itemize{
#'   \item \code{paths}: \code{sf} lines for all functions, with \code{funct},
#'     \code{destination}, \code{cost}, \code{cost_unit}, and \code{length}
#'     attributes;
#'   \item \code{paths.back}: return paths (NULL unless
#'     \code{return.base = TRUE});
#'   \item \code{cost.table}: data frame summarising cost and length per
#'     function and destination;
#'   \item \code{origin}, \code{destin}: the input locations (\code{sf});
#'   \item metadata used by the plot method.
#' }
#' Plot it at any time with \code{plot()} (see
#' \code{\link{plot.movecost_comp}}); \code{plot(x, type = "chart")} compares
#' the path length and the cost across functions, as boxplots of their
#' distribution across destinations when three or more destinations are used
#' (the \code{add.chart} output of movecost <= 2.x), or as bars otherwise.
#'
#' @examples
#' dtm <- mc_volc()
#' start <- mc_volc_loc()
#' destin <- mc_destin_loc()
#'
#' cmp <- mc_comp(dtm, origin = start, destin = destin[2, ],
#'                functs = c("t", "ks", "mp"), move = 8)
#'
#' cmp$cost.table
#' plot(cmp)                  # routes on the terrain
#' plot(cmp, type = "chart")  # comparative bar chart
#'
#' @seealso \code{\link{mc_surface}}, \code{\link{mc_paths}},
#'   \code{\link{plot.movecost_comp}}
#' @export
mc_comp <- function(dtm, origin, destin, functs, move = 16,
                    return.base = FALSE, time = "h", ...) {
  if (length(functs) < 2L) {
    stop("Provide at least two cost-function codes to compare.", call. = FALSE)
  }
  invisible(lapply(functs, .cf_get))   # validate all codes upfront

  paths_all <- list()
  backs_all <- list()
  dtm <- .as_spatraster(dtm)
  for (f in functs) {
    surf <- mc_surface(dtm = dtm, funct = f, move = move, ...)
    lcp <- mc_paths(surf, origin = origin, destin = destin,
                    return.base = return.base, time = time)
    cf <- .cf_get(f)
    p <- lcp$paths
    p$cost_hms <- NULL    # present only for time-based functions; drop so that
                          # paths of different families can be row-bound
    p$funct <- f
    p$destination <- seq_len(nrow(p))
    p$cost_unit <- if (cf$time_output) time else cf$unit
    paths_all[[f]] <- p
    if (isTRUE(return.base)) {
      pb <- lcp$paths.back
      pb$cost_hms <- NULL
      pb$funct <- f
      pb$destination <- seq_len(nrow(pb))
      pb$cost_unit <- if (cf$time_output) time else cf$unit
      backs_all[[f]] <- pb
    }
  }
  paths <- do.call(rbind, paths_all)
  rownames(paths) <- NULL
  paths.back <- if (isTRUE(return.base)) do.call(rbind, backs_all) else NULL

  cost.table <- sf::st_drop_geometry(paths)[, c("funct", "destination", "cost",
                                                "cost_unit", "length")]

  ori_sf <- .as_sf(origin, "origin")
  dst_sf <- .as_sf(destin, "destin")
  out <- list(paths = paths,
              paths.back = paths.back,
              cost.table = cost.table,
              origin = ori_sf,
              destin = dst_sf,
              functs = functs,
              time = time,
              dtm = dtm)
  class(out) <- "movecost_comp"
  out
}

#' @rdname movecost-methods
#' @export
print.movecost_comp <- function(x, ...) {
  cat("movecost cost-function comparison\n")
  cat("  functions     : ", paste(x$functs, collapse = ", "), "\n", sep = "")
  cat("  destinations  : ", nrow(x$destin), "\n", sep = "")
  print(x$cost.table)
  cat("Use plot() to visualise routes; plot(x, type = 'chart') for a bar chart.\n")
  invisible(x)
}
