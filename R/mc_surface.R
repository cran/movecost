#' Create a reusable cost surface (the computational core of movecost 3.0)
#'
#' \code{mc_surface()} builds, once, everything that any subsequent analysis needs:
#' the directed graph of movement costs between adjacent DTM cells, derived from a
#' user-selected cost function. The returned object is then fed to
#' \code{\link{mc_accum}}, \code{\link{mc_paths}}, \code{\link{mc_corridor}},
#' \code{\link{mc_boundary}}, \code{\link{mc_alloc}}, \code{\link{mc_network}},
#' and \code{\link{mc_rank}}, none of which recomputes the graph.\cr
#'
#' \strong{Why this design}: in movecost <= 2.x every analysis function rebuilt
#' the conductance matrix from scratch at each call (a network between n sites
#' implied n or more redundant constructions). Building the graph once and
#' reusing it makes multi-analysis workflows and large location sets dramatically
#' faster, and makes the package's internal structure modular: one function, one
#' responsibility.\cr
#'
#' \strong{Cost model}: for each pair of adjacent cells (4, 8, or 16 neighbours,
#' set via \code{move}), the terrain slope is computed as rise over run between
#' cell centres, signed in the direction of travel (positive uphill, negative
#' downhill), so movement costs are anisotropic. The user-selected cost function
#' (see \code{\link{mc_cost_functions}} for all 26 implemented functions, their
#' formulas, units, and references) maps the slope to either a walking speed or
#' the reciprocal of a per-metre cost; the edge cost is then the distance between
#' cell centres divided by speed (yielding seconds) or multiplied by the
#' per-metre cost. Accumulated costs and least-cost paths are obtained with
#' Dijkstra's algorithm via the \emph{igraph} library.\cr
#'
#' \strong{Cognitive slope}: if \code{cogn.slp = TRUE}, positive slopes are
#' multiplied by 1.99 and negative slopes by 2.31 before the cost function is
#' applied, following Pingel T.J. (2013), Modeling Slope as a Contributor to
#' Route Selection in Mountainous Areas, Cartography and Geographic Information
#' Science, 37(2), 137-148.\cr
#'
#' \strong{Topographic (surface) distance}: if \code{topo.dist = TRUE}, costs are
#' corrected by the factor \eqn{sqrt(1 + slope^2)}, accounting for the actual
#' distance travelled over sloped terrain rather than its planar projection.
#' The correction is clearly appropriate for time-based functions; for metabolic
#' functions, which were empirically derived against planar distance, it should
#' be used with caution. In gently undulating landscapes the difference is
#' negligible; it matters in rugged terrain with slopes beyond 20-30 percent.\cr
#'
#' \strong{Barriers}: areas where movement is inhibited can be supplied via
#' \code{barrier} (any sf, SpatVector, or legacy sp lines/polygons object).
#' Edges incident to cells touched by the barrier have their conductance
#' multiplied by \code{field}: 0 (default) makes the barrier impassable, while
#' intermediate values (e.g. 0.01) penalise but do not forbid crossing. Note
#' that with \code{move = 16}, knight-move connections can "jump" a thin linear
#' barrier, exactly as in movecost <= 2.x; use \code{move = 8} if that must not
#' happen.\cr
#'
#' \strong{Terrain factor raster}: the terrain factor N (see
#' \code{\link{mc_cost_functions}}) can vary across the study area. Supply a
#' raster via \code{terrain} in either of two forms:
#' \itemize{
#'   \item a \emph{numeric} raster whose cells hold terrain factors directly
#'     (e.g. 1.00 on roads, 1.79 on bogs);
#'   \item a \emph{categorical} raster of terrain classes (integer codes, or a
#'     terra factor raster with category labels), together with a lookup table
#'     \code{terrain.table} pairing each class with its terrain factor.
#'     \code{\link{mc_terrain_factors}} returns documented reference values
#'     from which the table can be built.
#' }
#' The scalar \code{N} keeps its role and acts as the \emph{fallback} value:
#' cells with no terrain value (NA in the raster, cells outside its extent, or
#' classes missing from the table) take \code{N}, so a raster can map only the
#' terrain types of interest (roads, bogs, ...) and leave the rest to a single
#' general value (e.g. \code{N = 1.67}, off-path). A message reports how many
#' cells were filled in this way; a warning lists classes missing from the
#' table. Terrain factors must be positive; \code{Inf} makes a cell
#' impassable.\cr
#' Each edge between two adjacent cells takes the \emph{mean} of their terrain
#' factors. Half of an edge lies in each cell, and for every implemented
#' function the edge cost is linear in N (affine for the Pandolf and Van Leusen
#' functions, where N multiplies only the terrain-dependent term of the
#' metabolic equation); the mean therefore reproduces exactly the cost of
#' crossing half of the edge on each terrain. Note that N is applied
#' \emph{inside} the cost function, as with the scalar N, and not as a
#' multiplier of the final cost: for the Pandolf-type functions the two are
#' different. With \code{move = 16}, knight-move edges pass close to two
#' further cells whose terrain is not considered (the same simplification
#' made for elevation).\cr
#' If the raster's coordinate system or grid differs from the DTM's, it is
#' aligned to the DTM with nearest-neighbour resampling
#' (\code{terra::project(method = "near")}), which never invents values and is
#' therefore safe for class codes; a warning is issued, as the data are
#' changed. For full control, align the raster beforehand.\cr
#' A few cost functions do not take a terrain factor (see the \code{uses_N}
#' column of \code{\link{mc_cost_functions}()}) and ignore both \code{N} and
#' \code{terrain}, with a warning: the off-path functions (\code{"tofp"},
#' \code{"icmoffp"}, \code{"icfoffp"}), which already embed off-path walking,
#' and Alberti's function (\code{"alb"}), which models the walking pace of a
#' flock rather than of a human walker. To model terrain with a raster, use the
#' on-path version of an off-path function (\code{"t"}, \code{"icmonp"},
#' \code{"icfonp"}), so that the terrain effect is not counted twice.\cr
#' A caution: published terrain factors were obtained in different studies,
#' under specific conditions and for specific cost models (see
#' \code{\link{mc_terrain_factors}}); a terrain raster makes their use easy,
#' but does not make them more accurate.\cr
#'
#' \strong{Irregular DTMs}: cells with NoData (for instance, sea around a
#' coastline) are simply excluded from the graph, so accumulated costs and
#' least-cost paths can never traverse them. The \code{irregular.dtm} workaround
#' of movecost <= 2.x is therefore no longer needed and has no equivalent here.\cr
#'
#' \strong{Coordinate system}: the DTM must use a projected coordinate system
#' with planar units in metres; geographic (lon/lat) rasters are rejected with
#' an informative error.
#'
#' @param dtm Digital Terrain Model: a terra \code{SpatRaster}, a legacy
#'   \code{RasterLayer}, or a file path readable by \code{terra::rast()}. If
#'   NULL, elevation data are downloaded for the area defined by
#'   \code{studyplot} via \code{\link{mc_dtm}} (requires the \emph{elevatr}
#'   package).
#' @param funct character: code of the cost function to be used (default
#'   \code{"t"}, Tobler's on-path hiking function). Run
#'   \code{\link{mc_cost_functions}()} for the complete annotated list.
#' @param move number of directions in which cells are connected: 4 (rook's
#'   case), 8 (queen's case), or 16 (knight and one-cell queen moves; default).
#' @param studyplot polygon (sf, SpatVector, or sp) defining the area for which
#'   online elevation data are downloaded when \code{dtm} is NULL.
#' @param z zoom level (resolution) of the downloaded elevation data (0-15;
#'   default 9); only used when \code{dtm} is NULL. See \code{\link{mc_dtm}}.
#' @param barrier lines or polygons (sf, SpatVector, or sp) representing areas
#'   where movement is inhibited (see Details).
#' @param field conductance multiplier applied to edges touching the barrier
#'   (0 by default, i.e. movement inhibited).
#' @param cogn.slp TRUE or FALSE (default): use the 'cognitive slope' in place
#'   of the real slope (see Details).
#' @param topo.dist TRUE or FALSE (default): correct costs for topographic
#'   (surface) distance rather than planar distance (see Details).
#' @param N terrain factor representing ease of movement (1 by default; see
#'   \code{\link{mc_cost_functions}} for reference values): a single positive
#'   number, applied to the whole study area. If \code{terrain} is supplied,
#'   \code{N} is used only for the cells with no terrain value (see Details).
#' @param sl.crit critical slope in percent (default 10), used only by the
#'   wheeled-vehicle cost function.
#' @param W walker's body weight in kg (default 70), used by the Pandolf and
#'   Van Leusen functions.
#' @param L carried load in kg (default 0), used by the Pandolf and Van Leusen
#'   functions.
#' @param V walking speed in m/s (default 1.2), used by the Pandolf, Van
#'   Leusen, and Ardigo functions; if 0, speed is derived cell-by-cell from the
#'   Tobler on-path function.
#' @param terrain optional raster of terrain factors or terrain classes
#'   (SpatRaster, RasterLayer, or a file path readable by
#'   \code{terra::rast()}), making the terrain factor vary across the study
#'   area (see Details). NULL (default) applies the scalar \code{N}
#'   everywhere.
#' @param terrain.table optional data frame translating the classes of a
#'   categorical \code{terrain} raster into terrain factors, with columns
#'   \code{class} (numeric codes, or the category labels of a terra factor
#'   raster) and \code{N} (terrain factor); a two-column data frame without
#'   these names is read as class, then N. Required for a terra factor raster;
#'   if NULL, the cell values of a numeric raster are taken as terrain factors.
#'
#' @return An object of class \code{movecost_surface}: a list with components
#' \itemize{
#'   \item \code{dtm}: the input DTM (\code{SpatRaster});
#'   \item \code{graph}: the directed igraph cost graph (one vertex per cell);
#'   \item \code{cost.raster}: \code{SpatRaster} storing, for each cell, the
#'     mean value of the cost function over the edges leaving the cell (walking
#'     speed in km/h for speed-based functions, conductance otherwise) - the
#'     analogue of the \code{cost.surface} component of movecost <= 2.x;
#'   \item \code{funct}, \code{params}, \code{move}: metadata describing how
#'     the graph was built (used by downstream functions and plot methods);
#'   \item \code{terrain}: \code{SpatRaster} of the terrain factor actually
#'     used for each cell (after alignment, class lookup, and filling with
#'     \code{N}), or NULL if no terrain raster was supplied; draw it with
#'     \code{plot(x, type = "terrain")}.
#' }
#' The object has \code{print()} and \code{plot()} methods; \code{plot()}
#' returns a ggplot2 object (see \code{\link{plot.movecost_surface}}).\cr
#'
#' \strong{Saving to disk}: because SpatRaster objects hold external pointers,
#' use \code{\link{mc_save}} and \code{\link{mc_load}} (not bare
#' \code{saveRDS()}) to store movecost objects across sessions.
#'
#' @examples
#' # load the package sample DTM (Maungawhau volcano) and locations
#' dtm <- mc_volc()
#' start <- mc_volc_loc()
#' destin <- mc_destin_loc()
#'
#' # build the cost surface once (Tobler's on-path function, 8 directions)
#' surf <- mc_surface(dtm, funct = "t", move = 8)
#'
#' # ...and reuse it across analyses without recomputation:
#' acc <- mc_accum(surf, origin = start, breaks = 0.05)
#' lcp <- mc_paths(surf, origin = start, destin = destin)
#'
#' # visualise on demand (no re-running of any computation):
#' plot(acc)
#' plot(lcp)
#'
#' # a terrain factor varying across the area: a (fictitious) categorical
#' # raster with class 1 below 120 m and class 2 (rock) above, translated
#' # into terrain factors via a lookup table
#' terr <- terra::classify(dtm, rbind(c(-Inf, 120, 1), c(120, Inf, 2)))
#' tf <- mc_terrain_factors("time")
#' tab <- data.frame(class = 1:2,
#'                   N = tf$N[match(c("Lawn grass", "Rock"), tf$terrain)])
#' surf_t <- mc_surface(dtm, funct = "t", move = 8,
#'                      terrain = terr, terrain.table = tab)
#' surf_t
#' plot(surf_t, type = "terrain")
#'
#' @references Alberti (2019) <doi:10.1016/j.softx.2019.100331>.
#' @seealso \code{\link{mc_cost_functions}}, \code{\link{mc_accum}},
#'   \code{\link{mc_paths}}, \code{\link{mc_corridor}}, \code{\link{mc_boundary}},
#'   \code{\link{mc_alloc}}, \code{\link{mc_network}}, \code{\link{mc_rank}},
#'   \code{\link{mc_dtm}}, \code{\link{mc_terrain_factors}},
#'   \code{\link{mc_save}};
#'   visualisation: \code{\link{plot.movecost_surface}}
#' @export
mc_surface <- function(dtm = NULL, funct = "t", move = 16,
                       studyplot = NULL, z = 9,
                       barrier = NULL, field = 0,
                       cogn.slp = FALSE, topo.dist = FALSE,
                       N = 1, sl.crit = 10, W = 70, L = 0, V = 1.2,
                       terrain = NULL, terrain.table = NULL) {

  # -- acquire and validate the DTM ------------------------------------------
  if (is.null(dtm)) {
    if (is.null(studyplot)) {
      stop("Either 'dtm' or 'studyplot' must be provided.", call. = FALSE)
    }
    dtm <- mc_dtm(studyplot = studyplot, z = z)
  }
  dtm <- .as_spatraster(dtm)
  if (terra::nlyr(dtm) > 1) {
    warning("The DTM has ", terra::nlyr(dtm),
            " layers; only the first is used.", call. = FALSE)
    dtm <- dtm[[1]]
  }
  .assert_projected(dtm)
  if (!move %in% c(4, 8, 16)) {
    stop("'move' must be 4, 8, or 16.", call. = FALSE)
  }
  if (!is.numeric(field) || field < 0 || field > 1) {
    stop("'field' must be a conductance multiplier between 0 and 1.",
         call. = FALSE)
  }

  if (!is.numeric(N) || length(N) != 1L || !is.finite(N) || N <= 0) {
    stop("'N' must be a single positive number (the terrain factor).",
         call. = FALSE)
  }
  if (is.null(terrain) && !is.null(terrain.table)) {
    stop("'terrain.table' was supplied without a 'terrain' raster.",
         call. = FALSE)
  }

  cf <- .cf_get(funct)
  params <- list(N = N, sl.crit = sl.crit, W = W, L = L, V = V,
                 cogn.slp = cogn.slp, topo.dist = topo.dist, field = field,
                 terrain = FALSE)

  # functions that do not take a terrain factor (the off-path functions, which
  # already embed off-path walking, and Alberti's, which models a flock's
  # walking pace) ignore N and 'terrain', with a warning saying why
  if (!isTRUE(cf$uses_N)) {
    if (!is.null(terrain)) {
      warning(.no_N_message(funct, "'terrain' is ignored", "raster"),
              call. = FALSE)
      terrain <- NULL
    } else if (N != 1) {
      warning(.no_N_message(funct, paste0("N = ", N, " is ignored"), "N"),
              call. = FALSE)
    }
  }

  # -- terrain raster: per-cell terrain factors (on the DTM grid) -------------
  terrain_rast <- NULL
  n_filled <- 0L
  if (!is.null(terrain)) {
    terrain <- .align_terrain(dtm, terrain)
    tv <- .terrain_values(dtm, terrain, terrain.table, N)
    n_filled <- tv$n.filled
    if (n_filled > 0) {
      message("Note: ", format(n_filled, big.mark = ","), " DTM cell(s) with ",
              "no terrain value take the scalar N = ", N, ".")
    }
    terrain_rast <- terra::rast(dtm)
    terra::values(terrain_rast) <- tv$values
    names(terrain_rast) <- "N"
    n_rng <- range(tv$values, na.rm = TRUE)
    params$terrain <- TRUE
    params$terrain.range <- n_rng
    params$terrain.filled <- n_filled
  }

  # -- adjacency: pairs of connected cells -----------------------------------
  n_cells <- terra::ncell(dtm)
  adj <- terra::adjacent(dtm, cells = seq_len(n_cells),
                         directions = move, pairs = TRUE)
  from <- adj[, 1]
  to <- adj[, 2]

  # elevations; drop edges touching NoData cells (sea, voids, ...)
  elev <- terra::values(dtm, mat = FALSE)
  keep <- !(is.na(elev[from]) | is.na(elev[to]))
  from <- from[keep]
  to <- to[keep]

  # -- distances and signed slope (rise over run, direction of travel) -------
  xy_from <- terra::xyFromCell(dtm, from)
  xy_to <- terra::xyFromCell(dtm, to)
  dist <- sqrt((xy_to[, 1] - xy_from[, 1])^2 + (xy_to[, 2] - xy_from[, 2])^2)
  slope <- (elev[to] - elev[from]) / dist

  # cognitive slope (Pingel 2013)
  if (isTRUE(cogn.slp)) {
    slope <- ifelse(slope > 0, slope * 1.99, slope * 2.31)
  }

  # -- apply the cost function -----------------------------------------------
  # with a terrain raster, N becomes a vector with one value per edge: the
  # mean of the two cells the edge joins. As half of each edge lies in each
  # cell, and every implemented function yields an edge cost that is linear
  # (or, for the Pandolf family, affine) in N, the mean reproduces exactly the
  # cost of crossing half of the edge on each terrain. All the cost-function
  # formulas are vectorised, so they accept the per-edge N unchanged.
  params_cf <- params
  if (!is.null(terrain_rast)) {
    n_cell <- terra::values(terrain_rast, mat = FALSE)
    params_cf$N <- (n_cell[from] + n_cell[to]) / 2
  }
  fun <- cf$make(params_cf)
  fvals <- fun(slope)

  # functions evaluated only down to a slope floor (Minetti's polynomial,
  # held at the lower end of its experimental range): report how many edges
  # were affected, so that the user can judge the extrapolation
  if (!is.null(cf$slope_floor)) {
    n_floor <- sum(slope < cf$slope_floor, na.rm = TRUE)
    if (n_floor > 0) {
      message("Note: ", format(n_floor, big.mark = ","), " of ",
              format(length(slope), big.mark = ","), " graph edges are steeper ",
              "downhill than ", cf$slope_floor, " (the lower end of the range ",
              "over which the '", funct, "' function was measured); their cost ",
              "is evaluated at ", cf$slope_floor, ". Consider 'hrz' for terrain ",
              "this steep.")
    }
  }

  # conductance: speed functions are km/h -> m/s (exact factor 1/3.6; movecost
  # <= 3.1.0 used the rounded 0.278, which shortened times by 0.08%),
  # then divided by distance; reciprocal functions divided by distance
  conductance <- if (cf$group == "speed") (fvals / 3.6) / dist else fvals / dist

  # topographic (surface) distance correction
  if (isTRUE(topo.dist)) {
    conductance <- conductance / sqrt(1 + slope^2)
  }

  # -- barrier ----------------------------------------------------------------
  if (!is.null(barrier)) {
    barrier_sf <- .align_crs(dtm, .as_sf(barrier, "barrier"), "barrier")
    bcells <- .cells_touched(dtm, barrier_sf)
    hit <- (from %in% bcells) | (to %in% bcells)
    conductance[hit] <- conductance[hit] * field
  }

  # -- per-cell summary raster (analogue of v2.x 'cost.surface') -------------
  ok <- is.finite(fvals)
  cell_sum <- rowsum(fvals[ok], group = from[ok])
  cell_n <- rowsum(rep(1, sum(ok)), group = from[ok])
  cost_vals <- rep(NA_real_, n_cells)
  cost_vals[as.integer(rownames(cell_sum))] <- cell_sum / cell_n
  cost.raster <- terra::rast(dtm)
  terra::values(cost.raster) <- cost_vals
  names(cost.raster) <- if (cf$group == "speed") "speed_kmh" else "conductance"

  # -- build the directed graph ----------------------------------------------
  # edges with zero/NaN conductance are dropped (impassable);
  # infinite conductance (e.g. the 'ree' function on flat ground, whose cost
  # tends to zero) yields zero-cost edges, which Dijkstra handles correctly
  pos <- !is.na(conductance) & conductance > 0
  w <- 1 / conductance[pos]          # Inf conductance -> 0 weight
  g <- igraph::make_empty_graph(n = n_cells, directed = TRUE)
  g <- igraph::add_edges(g, rbind(from[pos], to[pos]), weight = w)

  # label: with a terrain raster, report the range of N instead of a value
  params_lab <- params
  if (!is.null(terrain_rast)) {
    params_lab$N <- paste0("raster (", .fmt_num(n_rng[1]), " to ",
                           .fmt_num(n_rng[2]), ")")
  }

  out <- list(dtm = dtm,
              graph = g,
              cost.raster = cost.raster,
              funct = funct,
              funct.label = cf$label(params_lab),
              params = params,
              move = move,
              barrier = if (is.null(barrier)) NULL else barrier_sf,
              terrain = terrain_rast)
  class(out) <- "movecost_surface"
  out
}

#' @rdname movecost-methods
#' @export
print.movecost_surface <- function(x, ...) {
  cf <- .cf_get(x$funct)
  cat("movecost cost surface\n")
  cat("  cost function : ", cf$name, " ('", x$funct, "')\n", sep = "")
  cat("  cost unit     : ", cf$unit, "\n", sep = "")
  cat("  DTM           : ", terra::nrow(x$dtm), " x ", terra::ncol(x$dtm),
      " cells, resolution ", paste(round(terra::res(x$dtm), 2), collapse = " x "),
      "\n", sep = "")
  cat("  connectivity  : ", x$move, " directions; ",
      format(igraph::ecount(x$graph), big.mark = ","), " graph edges\n", sep = "")
  if (isTRUE(x$params$terrain)) {
    rng <- x$params$terrain.range
    cat("  terrain factor: raster, N from ", .fmt_num(rng[1]), " to ",
        .fmt_num(rng[2]), sep = "")
    if (isTRUE(x$params$terrain.filled > 0)) {
      cat("; ", format(x$params$terrain.filled, big.mark = ","),
          " cell(s) with scalar N = ", x$params$N, sep = "")
    }
    cat("\n")
  } else if (isTRUE(.cf_get(x$funct)$uses_N)) {
    cat("  terrain factor: N = ", x$params$N, "\n", sep = "")
  }
  if (!is.null(x$barrier)) cat("  barrier       : yes (field = ", x$params$field, ")\n", sep = "")
  if (isTRUE(x$params$cogn.slp)) cat("  cognitive slope: yes\n")
  if (isTRUE(x$params$topo.dist)) cat("  topographic distance correction: yes\n")
  cat("Reuse this object with mc_accum(), mc_paths(), mc_corridor(), mc_boundary(),\n")
  cat("mc_alloc(), mc_network(), mc_rank(); plot() it for a quick look.\n")
  invisible(x)
}

#' Internal: accumulated cost values (one row per origin cell)
#'
#' Single igraph::distances() call; unit conversion applied afterwards.
#' Returns a matrix [n_origins x n_cells] in the cost function's native unit
#' (seconds for time-output functions; J, J/kg, kJ, kcal x 1000 or abstract
#' units otherwise, before division by the registry 'divisor').
#'
#' @keywords internal
#' @noRd
.accum_matrix <- function(surface, origin_cells, mode = "out") {
  igraph::distances(surface$graph,
                    v = origin_cells,
                    mode = mode,
                    algorithm = "dijkstra")
}

#' Internal: convert native accumulated values to user-facing units
#' @keywords internal
#' @noRd
.convert_units <- function(vals, funct, time = "h") {
  cf <- .cf_get(funct)
  if (cf$time_output) vals <- .time_convert(vals, time)
  if (!is.null(cf$divisor) && cf$divisor != 1) vals <- vals / cf$divisor
  vals
}

#' Internal: cells of point locations, with validation
#' @keywords internal
#' @noRd
.points_to_cells <- function(surface, pts, arg_name = "origin") {
  pts_sf <- .align_crs(surface$dtm, .as_sf(pts, arg_name), arg_name)
  if (!all(sf::st_geometry_type(pts_sf) == "POINT")) {
    stop("'", arg_name, "' must contain POINT geometries only ",
         "(cast MULTIPOINTs first, e.g. sf::st_cast(x, 'POINT')).",
         call. = FALSE)
  }
  xy <- sf::st_coordinates(pts_sf)[, 1:2, drop = FALSE]
  cells <- terra::cellFromXY(surface$dtm, xy)
  if (anyNA(cells)) {
    stop("Some '", arg_name, "' locations fall outside the DTM extent.", call. = FALSE)
  }
  elev <- terra::values(surface$dtm, mat = FALSE)[cells]
  if (anyNA(elev)) {
    stop("Some '", arg_name, "' locations fall on NoData cells of the DTM.", call. = FALSE)
  }
  list(sf = pts_sf, cells = cells)
}
