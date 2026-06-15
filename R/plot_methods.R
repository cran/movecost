#' ggplot2-based visualisation of movecost results
#'
#' Every analysis function in movecost 3.0 returns an object that stores its
#' results; visualisation is performed afterwards, on demand, by the
#' \code{plot()} methods documented here. Each method returns a
#' \strong{ggplot object}: printed at the console it draws the map, assigned to
#' a variable it can be customised with standard ggplot2 syntax (additional
#' layers, themes, scales) and drawn any number of times - all without
#' re-running any computation. \code{ggplot2::autoplot()} methods are provided
#' as aliases.
#'
#' @name movecost-plots
#' @keywords internal
NULL

#' Internal S3 methods for movecost objects
#'
#' \code{print} and \code{autoplot} methods for the objects returned by the
#' package's analysis functions. The \code{print} methods give a concise
#' console summary of a result; the \code{autoplot} methods are thin aliases of
#' the corresponding \code{plot} methods (for example
#' \code{\link{plot.movecost_accum}}), provided so that the package works with
#' the \code{ggplot2::autoplot()} generic. These methods are grouped here and
#' kept out of the main help index to reduce clutter; they remain fully
#' documented and accessible, for example via \code{?print.movecost_accum} or
#' \code{?autoplot.movecost_accum}.
#'
#' @param x a movecost result object (for the \code{print} methods).
#' @param object a movecost result object (for the \code{autoplot} methods).
#' @param ... passed on to the corresponding \code{plot} method (for
#'   \code{autoplot}), or unused (for \code{print}).
#' @return The \code{print} methods return their argument invisibly; the
#'   \code{autoplot} methods return a ggplot object.
#' @name movecost-methods
#' @keywords internal
NULL

# tell R CMD check that the aes() column names are data variables
utils::globalVariables(c("x", "y", "value", "shade", "level", "rank",
                         "funct", "cost", "metric", "origin_id", "label",
                         "destination", "direction", "length", "lab"))

#' Internal: raster to data frame for ggplot
#' @keywords internal
#' @noRd
.r_df <- function(r, name = "value") {
  df <- terra::as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(df) <- c("x", "y", name)
  df
}

#' Internal: slope-shade data frame (used as a translucent overlay)
#' @keywords internal
#' @noRd
.hs_df <- function(dtm) {
  sl <- terra::terrain(dtm, v = "slope", unit = "radians")
  df <- terra::as.data.frame(sl, xy = TRUE, na.rm = TRUE)
  names(df) <- c("x", "y", "shade")
  df
}

#' Internal: common map scaffolding (coord, theme, titles)
#' @keywords internal
#' @noRd
.gg_finish <- function(p, dtm, title, subtitle) {
  crs <- sf::st_crs(terra::crs(dtm))
  p +
    ggplot2::coord_sf(crs = crs, datum = crs, expand = FALSE) +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 11, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 8, colour = "grey30"),
      axis.text = ggplot2::element_text(size = 6)
    )
}

#' Internal: translucent slope-shade layer
#' @keywords internal
#' @noRd
.gg_shade <- function(dtm, transp = 0.5) {
  list(
    ggplot2::geom_raster(data = .hs_df(dtm),
                         ggplot2::aes(x = x, y = y, alpha = shade),
                         fill = "black"),
    ggplot2::scale_alpha_continuous(range = c(0, transp), guide = "none")
  )
}

#' Internal: elevation base layer (terrain colours)
#' @keywords internal
#' @noRd
.gg_elev <- function(dtm) {
  list(
    ggplot2::geom_raster(data = .r_df(dtm, "value"),
                         ggplot2::aes(x = x, y = y, fill = value)),
    ggplot2::scale_fill_gradientn(colours = grDevices::terrain.colors(255),
                                  name = "Elevation\n(masl)")
  )
}

# ---------------------------------------------------------------------------
# mc_surface
# ---------------------------------------------------------------------------

#' Plot a movecost cost surface
#'
#' Returns a ggplot object showing the per-cell cost raster (walking speed in
#' km/h for speed-based functions, conductance otherwise) with a translucent
#' slope-shade overlay.
#'
#' @param x a \code{movecost_surface} object.
#' @param transp transparency of the slope-shade overlay (0.5 by default).
#' @param ... unused.
#' @seealso \code{\link{mc_surface}}
#' @return A ggplot object (draw it by printing; customise it with \code{+}).
#' @export
plot.movecost_surface <- function(x, transp = 0.5, ...) {
  cf <- .cf_get(x$funct)
  legend <- if (cf$group == "speed") "speed (km/h)" else "conductance"
  p <- ggplot2::ggplot() +
    ggplot2::geom_raster(data = .r_df(x$cost.raster, "value"),
                         ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::scale_fill_viridis_c(name = legend, option = "C")
  p <- p + .gg_shade(x$dtm, transp)
  .gg_finish(p, x$dtm, "Cost of movement across the terrain", x$funct.label)
}

#' @rdname movecost-methods
#' @export
autoplot.movecost_surface <- function(object, ...) plot(object, ...)

# ---------------------------------------------------------------------------
# mc_accum
# ---------------------------------------------------------------------------

#' Plot an accumulated cost surface
#'
#' Returns a ggplot object showing the accumulated cost surface with isolines
#' and origin location(s), or the isolines alone
#' (\code{type = "contour"}), reproducing the two output types of
#' movecost <= 2.x's \code{outp} parameter - but on demand, without re-running
#' the analysis.
#'
#' @param x a \code{movecost_accum} object.
#' @param type "raster" (default; accumulated surface, slope-shade, isolines)
#'   or "contour" (isolines only).
#' @param cont.lab TRUE (default) or FALSE: label the isolines.
#' @param transp transparency of the slope-shade overlay (0.5 by default).
#' @param cex.breaks size of the isoline labels (2.5 by default, in ggplot2
#'   text-size units).
#' @param iso.col colour of the isolines and their labels (default
#'   \code{"white"}, which reads clearly against the dark low-cost core of the
#'   viridis "C" palette; over the \code{"contour"} output, where there is no
#'   coloured background, \code{"grey20"} is used automatically unless the user
#'   overrides \code{iso.col}).
#' @param iso.lwd line width of the isolines (0.3 by default).
#' @param ... unused.
#' @seealso \code{\link{mc_accum}}
#' @return A ggplot object (draw it by printing; customise it with \code{+}).
#' @export
plot.movecost_accum <- function(x, type = c("raster", "contour"),
                                cont.lab = TRUE, transp = 0.5,
                                cex.breaks = 2.5, iso.col = "white",
                                iso.lwd = 0.3, ...) {
  type <- match.arg(type)
  cf <- .cf_get(x$funct)
  legend <- if (cf$time_output) paste0("cost (", x$time, ")")
            else paste0("cost (", cf$unit, ")")
  title <- if (cf$time_output) {
    paste0("Walking-time isochrones (in ", x$time, ") around origin")
  } else {
    "Accumulated cost isolines around origin"
  }
  # on the plain "contour" output there is no dark raster to read white
  # lines against, so fall back to a dark colour unless the user chose one
  if (type == "contour" && missing(iso.col)) iso.col <- "grey20"

  p <- ggplot2::ggplot()
  if (type == "raster") {
    p <- p +
      ggplot2::geom_raster(data = .r_df(x$accum, "value"),
                           ggplot2::aes(x = x, y = y, fill = value)) +
      ggplot2::scale_fill_viridis_c(name = legend, option = "C")
    p <- p + .gg_shade(x$dtm, transp)
  }
  if (!is.null(x$isolines)) {
    p <- p + ggplot2::geom_sf(data = x$isolines, colour = iso.col,
                              linewidth = iso.lwd)
    if (isTRUE(cont.lab)) {
      labpts <- suppressWarnings(sf::st_point_on_surface(x$isolines))
      coords <- sf::st_coordinates(labpts)
      labdf <- data.frame(x = coords[, 1], y = coords[, 2],
                          label = labpts$level)
      p <- p + ggplot2::geom_text(data = labdf,
                                  ggplot2::aes(x = x, y = y, label = label),
                                  size = cex.breaks, colour = iso.col)
    }
  }
  p <- p + ggplot2::geom_sf(data = x$origin, shape = 20, size = 2,
                            colour = "black")
  .gg_finish(p, x$dtm, title, x$funct.label)
}

#' @rdname movecost-methods
#' @export
autoplot.movecost_accum <- function(object, ...) plot(object, ...)

# ---------------------------------------------------------------------------
# mc_paths
# ---------------------------------------------------------------------------

#' Plot least-cost paths
#'
#' Returns a ggplot object showing the DTM (terrain colours plus slope-shade),
#' the least-cost path(s), the origin (black dot) and destination(s) (red
#' dots), optionally labelled with the accumulated cost, and the return paths
#' (dashed) if they were computed.
#'
#' @param x a \code{movecost_paths} object.
#' @param destin.lab TRUE (default) or FALSE: label each destination with its
#'   accumulated cost.
#' @param cex.lcp.lab size of the destination labels (2.5 by default).
#' @param transp transparency of the slope-shade overlay (0.5 by default).
#' @param lcp.lwd line width of the least-cost path(s) (0.3 by default; the
#'   return paths, if any, are drawn at 80\% of this width).
#' @param plot.barrier TRUE or FALSE (default): draw the barrier (blue) if one
#'   was used when building the surface.
#' @param ... unused.
#' @seealso \code{\link{mc_paths}}
#' @return A ggplot object (draw it by printing; customise it with \code{+}).
#' @export
plot.movecost_paths <- function(x, destin.lab = TRUE, cex.lcp.lab = 2.5,
                                transp = 0.5, lcp.lwd = 0.3,
                                plot.barrier = FALSE, ...) {
  p <- ggplot2::ggplot()
  p <- p + .gg_elev(x$dtm) + .gg_shade(x$dtm, transp)
  if (isTRUE(plot.barrier) && !is.null(x$barrier)) {
    p <- p + ggplot2::geom_sf(data = x$barrier, colour = "blue",
                              linewidth = 0.6)
  }
  p <- p + ggplot2::geom_sf(data = x$paths, colour = "black",
                            linewidth = lcp.lwd)
  if (!is.null(x$paths.back)) {
    p <- p + ggplot2::geom_sf(data = x$paths.back, colour = "black",
                              linewidth = lcp.lwd * 0.8, linetype = "dashed")
  }
  p <- p +
    ggplot2::geom_sf(data = x$origin, shape = 20, size = 2.5, colour = "black") +
    ggplot2::geom_sf(data = x$destin, shape = 20, size = 2.5, colour = "red")
  if (isTRUE(destin.lab)) {
    coords <- sf::st_coordinates(x$destin)
    labs <- if (!is.null(x$destin$cost_hms)) x$destin$cost_hms
            else round(x$destin$cost, 2)
    labdf <- data.frame(x = coords[, 1], y = coords[, 2], label = labs)
    p <- p + ggplot2::geom_text(data = labdf,
                                ggplot2::aes(x = x, y = y, label = label),
                                size = cex.lcp.lab, hjust = -0.15,
                                colour = "grey10")
  }
  .gg_finish(p, x$dtm, "Digital Terrain Model with Least-Cost Path(s)",
             paste0("LCP(s) based on the ", x$funct.label,
                    "\nblack dot = origin; red dot(s) = destination(s)"))
}

#' @rdname movecost-methods
#' @export
autoplot.movecost_paths <- function(object, ...) plot(object, ...)

# ---------------------------------------------------------------------------
# mc_corridor
# ---------------------------------------------------------------------------

#' Plot a least-cost corridor
#'
#' Returns a ggplot object showing the corridor surface (low values = cheap
#' movement) with the input locations. For the two-location case computed with
#' \code{lcp = TRUE}, the least-cost paths A-to-B and B-to-A are overlaid (they
#' differ under anisotropic cost functions).
#'
#' @param x a \code{movecost_corridor} object.
#' @param transp transparency of the slope-shade overlay (0.5 by default).
#' @param lcp.lwd line width of the overlaid least-cost paths (0.4 by default).
#' @param ... unused.
#' @seealso \code{\link{mc_corridor}}
#' @return A ggplot object (draw it by printing; customise it with \code{+}).
#' @export
plot.movecost_corridor <- function(x, transp = 0.5, lcp.lwd = 0.4, ...) {
  p <- ggplot2::ggplot() +
    ggplot2::geom_raster(data = .r_df(x$corridor, "value"),
                         ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::scale_fill_viridis_c(name = paste0("cost (", x$cost.unit, ")"),
                                  option = "C", direction = -1)
  p <- p + .gg_shade(x$dtm, transp)
  # overlay the two directional least-cost paths, if computed (two-location case)
  if (!is.null(x$lcp.AtoB)) {
    p <- p + ggplot2::geom_sf(data = x$lcp.AtoB, colour = "white",
                              linewidth = lcp.lwd)
  }
  if (!is.null(x$lcp.BtoA)) {
    p <- p + ggplot2::geom_sf(data = x$lcp.BtoA, colour = "white",
                              linewidth = lcp.lwd, linetype = "dashed")
  }
  p <- p + ggplot2::geom_sf(data = x$locations, shape = 20, size = 2.5,
                            colour = "red")
  subtitle <- paste0(x$funct.label,
                     if (x$method == "through") "\n(method: through; directional A-to-B)"
                     else "\n(method: reach; symmetric)")
  .gg_finish(p, x$dtm, "Least-cost corridor", subtitle)
}

#' @rdname movecost-methods
#' @export
autoplot.movecost_corridor <- function(object, ...) plot(object, ...)

# ---------------------------------------------------------------------------
# mc_boundary
# ---------------------------------------------------------------------------

#' Plot cost-limit boundaries
#'
#' Returns a ggplot object showing the DTM with the cost-limit boundary
#' polygon of each origin.
#'
#' @param x a \code{movecost_boundary} object.
#' @param transp transparency of the slope-shade overlay (0.5 by default).
#' @param ... unused.
#' @seealso \code{\link{mc_boundary}}
#' @return A ggplot object (draw it by printing; customise it with \code{+}).
#' @export
plot.movecost_boundary <- function(x, transp = 0.5, ...) {
  bnd <- x$boundaries
  bnd$origin_id <- factor(bnd$origin_id)
  p <- ggplot2::ggplot() + .gg_elev(x$dtm) + .gg_shade(x$dtm, transp) +
    ggplot2::geom_sf(data = bnd,
                     ggplot2::aes(colour = origin_id),
                     fill = NA, linewidth = 0.6) +
    ggplot2::scale_colour_viridis_d(name = "origin", option = "H") +
    ggplot2::geom_sf(data = x$origin, shape = 20, size = 2.5, colour = "black")
  .gg_finish(p, x$dtm,
             paste0("Cost-limit boundaries (limit: ", x$limit, " ",
                    x$cost.unit, ")"),
             x$funct.label)
}

#' @rdname movecost-methods
#' @export
autoplot.movecost_boundary <- function(object, ...) plot(object, ...)

# ---------------------------------------------------------------------------
# mc_alloc
# ---------------------------------------------------------------------------

#' Plot cost-allocation zones
#'
#' Returns a ggplot object showing the cost-allocation zones (one colour per
#' origin), with the origins and, optionally, the isolines of the minimum
#' accumulated cost surface.
#'
#' @param x a \code{movecost_alloc} object.
#' @param isolines TRUE or FALSE (default): overlay the cost isolines.
#' @param transp transparency of the slope-shade overlay (0.5 by default).
#' @param iso.col colour of the isolines (default \code{"white"}, which reads
#'   clearly against the coloured allocation zones).
#' @param iso.lwd line width of the isolines (0.3 by default).
#' @param origin.lab TRUE (default) or FALSE: label each origin with its index.
#'   The indices match the zone identifiers (the cell values of the allocation
#'   raster and the \code{origin_id} of \code{zones}), so each origin can be
#'   related to the zone allocated to it.
#' @param cex.origin.lab size of the origin labels (3 by default, in ggplot2
#'   text-size units).
#' @param ... unused.
#' @seealso \code{\link{mc_alloc}}
#' @return A ggplot object (draw it by printing; customise it with \code{+}).
#' @export
plot.movecost_alloc <- function(x, isolines = FALSE, transp = 0.5,
                                iso.col = "white", iso.lwd = 0.3,
                                origin.lab = TRUE, cex.origin.lab = 3, ...) {
  zdf <- .r_df(x$alloc, "value")
  zdf$value <- factor(zdf$value)
  p <- ggplot2::ggplot() +
    ggplot2::geom_raster(data = zdf,
                         ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::scale_fill_viridis_d(name = "allocated to\norigin", option = "H",
                                  alpha = 0.8)
  p <- p + .gg_shade(x$dtm, transp)
  if (isTRUE(isolines) && !is.null(x$isolines)) {
    p <- p + ggplot2::geom_sf(data = x$isolines, colour = iso.col,
                              linewidth = iso.lwd)
  }
  p <- p + ggplot2::geom_sf(data = x$origin, shape = 20, size = 2.5,
                            colour = "black")
  if (isTRUE(origin.lab)) {
    oxy <- sf::st_coordinates(x$origin)
    odf <- data.frame(x = oxy[, 1], y = oxy[, 2],
                      label = seq_len(nrow(x$origin)))
    p <- p + ggplot2::geom_text(data = odf,
                                ggplot2::aes(x = x, y = y, label = label),
                                size = cex.origin.lab, colour = "black",
                                fontface = "bold", hjust = -0.3, vjust = -0.3)
  }
  .gg_finish(p, x$dtm, "Cost allocation zones", x$funct.label)
}

#' @rdname movecost-methods
#' @export
autoplot.movecost_alloc <- function(object, ...) plot(object, ...)

# ---------------------------------------------------------------------------
# mc_network
# ---------------------------------------------------------------------------

#' Plot a least-cost network
#'
#' Returns a ggplot object showing either the network of least-cost paths over
#' the DTM (\code{type = "paths"}, default) or the path-density raster
#' (\code{type = "density"}, available if the network was computed with
#' \code{density = TRUE}).
#'
#' @param x a \code{movecost_network} object.
#' @param type "paths" (default) or "density".
#' @param transp transparency of the slope-shade overlay (0.5 by default).
#' @param lcp.lwd line width of the network paths (0.3 by default).
#' @param node.lab TRUE (default) or FALSE: label each node with its index.
#'   The indices match the row and column names of the returned
#'   \code{cost.matrix}, so a node on the map can be related to its entries in
#'   the matrix (reproducing the behaviour of movecost <= 2.x).
#' @param cex.node.lab size of the node labels (3 by default, in ggplot2
#'   text-size units).
#' @param ... unused.
#' @seealso \code{\link{mc_network}}
#' @return A ggplot object (draw it by printing; customise it with \code{+}).
#' @export
plot.movecost_network <- function(x, type = c("paths", "density"),
                                  transp = 0.5, lcp.lwd = 0.3,
                                  node.lab = TRUE, cex.node.lab = 3, ...) {
  type <- match.arg(type)
  # node coordinates and indices, in the order used by the cost matrix
  node_xy <- sf::st_coordinates(x$nodes)
  node_df <- data.frame(x = node_xy[, 1], y = node_xy[, 2],
                        label = seq_len(nrow(x$nodes)))

  if (type == "density") {
    if (is.null(x$density.perc)) {
      stop("No density raster stored; re-run mc_network(..., density = TRUE).",
           call. = FALSE)
    }
    # density raster first, then the translucent slope-shade on top (so the
    # terrain shows through and gives context), as in movecost <= 2.x
    p <- ggplot2::ggplot() +
      ggplot2::geom_raster(data = .r_df(x$density.perc, "value"),
                           ggplot2::aes(x = x, y = y, fill = value)) +
      ggplot2::scale_fill_viridis_c(name = "% of paths", option = "C") +
      .gg_shade(x$dtm, transp) +
      ggplot2::geom_sf(data = x$nodes, shape = 20, size = 2, colour = "red")
    if (isTRUE(node.lab)) {
      p <- p + ggplot2::geom_text(data = node_df,
                                  ggplot2::aes(x = x, y = y, label = label),
                                  size = cex.node.lab, colour = "white",
                                  hjust = -0.3, vjust = -0.3)
    }
    return(.gg_finish(p, x$dtm, "Least-cost paths density", x$funct.label))
  }
  p <- ggplot2::ggplot() + .gg_elev(x$dtm) + .gg_shade(x$dtm, transp) +
    ggplot2::geom_sf(data = x$paths, colour = "black", linewidth = lcp.lwd) +
    ggplot2::geom_sf(data = x$nodes, shape = 20, size = 2.5, colour = "red")
  if (isTRUE(node.lab)) {
    p <- p + ggplot2::geom_text(data = node_df,
                                ggplot2::aes(x = x, y = y, label = label),
                                size = cex.node.lab, colour = "black",
                                fontface = "bold", hjust = -0.3, vjust = -0.3)
  }
  .gg_finish(p, x$dtm,
             paste0("Least-cost network (",
                    if (x$type == "allpairs") "all pairs" else "nearest neighbours",
                    ")"),
             x$funct.label)
}

#' @rdname movecost-methods
#' @export
autoplot.movecost_network <- function(object, ...) plot(object, ...)

# ---------------------------------------------------------------------------
# mc_rank
# ---------------------------------------------------------------------------

#' Plot ranked least-cost paths
#'
#' Returns a ggplot object showing, depending on \code{type}, the ranked paths
#' over the terrain (\code{"map"}, default), the same paths drawn over the
#' least-cost corridor between the two locations (\code{"corridor"}, the
#' equivalent of movecost <= 2.x's \code{use.corr}), or a bubble chart of path
#' length against rank with bubble size proportional to cost (\code{"chart"},
#' the equivalent of movecost <= 2.x's \code{add.chart}). No manual data
#' assembly is required for any of them.
#'
#' @param x a \code{movecost_rank} object.
#' @param type \code{"map"} (default), \code{"corridor"}, or \code{"chart"}
#'   (see Details).
#' @param transp transparency of the slope-shade overlay (0.5 by default).
#' @param lcp.lwd line width of the ranked paths (0.6 by default).
#' @param ... unused.
#' @seealso \code{\link{mc_rank}}
#' @return A ggplot object (draw it by printing; customise it with \code{+}).
#' @export
plot.movecost_rank <- function(x, type = c("map", "corridor", "chart"),
                               transp = 0.5, lcp.lwd = 0.6, ...) {
  type <- match.arg(type)

  if (type == "chart") {
    df <- sf::st_drop_geometry(x$paths)
    df$length <- as.numeric(df$length)
    # bubble size encodes the cost, but instead of a size legend the actual
    # cost value is printed next to each bubble; for time-based functions use
    # the sexagesimal (hh:mm:ss) form, matching the returned object
    df$lab <- if (!is.null(df$cost_hms)) {
      df$cost_hms
    } else {
      paste0(signif(df$cost, 3), " ", x$cost.unit)
    }
    p <- ggplot2::ggplot(df, ggplot2::aes(x = rank, y = length, size = cost)) +
      ggplot2::geom_point(colour = "steelblue", alpha = 0.85,
                          show.legend = FALSE) +
      ggplot2::geom_text(ggplot2::aes(label = lab), size = 3, hjust = -0.25,
                         vjust = 0.5, show.legend = FALSE) +
      ggplot2::scale_size_continuous(range = c(3, 10), guide = "none") +
      ggplot2::scale_x_continuous(breaks = df$rank,
                                  expand = ggplot2::expansion(mult = c(0.1, 0.25))) +
      ggplot2::labs(title = "Ranked least-cost paths: length and cost by rank",
                    subtitle = paste0(x$funct.label,
                                      "\nbubble size and label = cost"),
                    x = "rank (1 = optimal)", y = "length (m)") +
      ggplot2::theme_minimal(base_size = 10)
    return(p)
  }

  pp <- x$paths
  pp$rank <- factor(pp$rank)

  if (type == "corridor") {
    # build the least-cost corridor between the two locations on the fly, so
    # the user does not have to compute and overlay it manually
    co <- x$corridor
    p <- ggplot2::ggplot() +
      ggplot2::geom_raster(data = .r_df(co, "value"),
                           ggplot2::aes(x = x, y = y, fill = value)) +
      ggplot2::scale_fill_viridis_c(name = paste0("corridor cost (", x$cost.unit, ")"),
                                    option = "C", direction = -1) +
      .gg_shade(x$dtm, transp)
    title <- "Ranked least-cost paths over the least-cost corridor"
  } else {
    p <- ggplot2::ggplot() + .gg_elev(x$dtm) + .gg_shade(x$dtm, transp)
    title <- "Ranked least-cost paths"
  }

  p <- p +
    ggplot2::geom_sf(data = pp, ggplot2::aes(colour = rank), linewidth = lcp.lwd) +
    ggplot2::scale_colour_viridis_d(name = "rank\n(1 = optimal)", option = "H") +
    ggplot2::geom_sf(data = x$origin, shape = 20, size = 2.5, colour = "black") +
    ggplot2::geom_sf(data = x$destin, shape = 20, size = 2.5, colour = "red")
  .gg_finish(p, x$dtm, title, x$funct.label)
}

#' @rdname movecost-methods
#' @export
autoplot.movecost_rank <- function(object, ...) plot(object, ...)

# ---------------------------------------------------------------------------
# mc_comp
# ---------------------------------------------------------------------------

#' Plot a cost-function comparison
#'
#' Returns a ggplot object showing either the least-cost paths produced by the
#' compared cost functions, colour-coded by function (\code{type = "paths"},
#' default), or a chart comparing the path length and the accumulated cost
#' across functions (\code{type = "chart"}).\cr
#'
#' For the chart, when several destinations are compared (three or more) the
#' panels are \strong{boxplots} showing, for each cost function, the
#' distribution across destinations of the path length and of the cost,
#' reproducing the \code{add.chart} output of movecost <= 2.x's
#' \code{movecomp()}. With fewer than three destinations a boxplot would be
#' degenerate, so bars are drawn instead. The cost panel is split by cost unit,
#' so cost is only ever compared between functions that share a unit (length,
#' in metres, is always directly comparable). When the comparison was computed
#' with \code{return.base = TRUE}, the outward and the return paths are shown in
#' separate rows ("outward" and "return"), as movecost <= 2.x did with a
#' separate set of boxplots for the paths back to the origin.
#'
#' @param x a \code{movecost_comp} object.
#' @param type "paths" (default) or "chart".
#' @param transp transparency of the slope-shade overlay (0.5 by default).
#' @param ... unused.
#' @seealso \code{\link{mc_comp}}
#' @return A ggplot object (draw it by printing; customise it with \code{+}).
#' @export
plot.movecost_comp <- function(x, type = c("paths", "chart"),
                               transp = 0.5, ...) {
  type <- match.arg(type)
  if (type == "chart") {
    n_dest <- nrow(x$destin)
    # assemble a long data frame of length and cost values, by function, and
    # by movement direction (outward, and return if return.base was used)
    .long <- function(tab, direction) {
      rbind(
        data.frame(funct = tab$funct,
                   metric = paste0("cost (", tab$cost_unit, ")"),
                   value = tab$cost, direction = direction),
        data.frame(funct = tab$funct,
                   metric = "length (m)",
                   value = as.numeric(tab$length), direction = direction)
      )
    }
    df <- .long(x$cost.table, "outward")
    has_back <- !is.null(x$paths.back)
    if (has_back) {
      df <- rbind(df, .long(sf::st_drop_geometry(x$paths.back), "return"))
    }
    df$funct <- factor(df$funct, levels = x$functs)
    df$direction <- factor(df$direction, levels = c("outward", "return"))

    p <- ggplot2::ggplot(df, ggplot2::aes(x = funct, y = value, fill = funct))
    # boxplots when there are enough destinations to form a distribution
    # (the movecost <= 2.x 'add.chart' behaviour); bars otherwise
    p <- p + if (n_dest >= 3) {
      ggplot2::geom_boxplot(show.legend = FALSE)
    } else {
      ggplot2::geom_col(show.legend = FALSE)
    }
    p <- p + ggplot2::scale_fill_viridis_d(option = "H")
    # outward and return paths in separate rows when both are present
    p <- p + if (has_back) {
      ggplot2::facet_grid(direction ~ metric, scales = "free_y")
    } else {
      ggplot2::facet_wrap(~ metric, scales = "free_y")
    }
    p <- p +
      ggplot2::labs(
        title = "Cost-function comparison",
        subtitle = if (n_dest >= 3)
          paste0("distribution across ", n_dest, " destinations")
        else NULL,
        x = "cost function", y = NULL) +
      ggplot2::theme_minimal(base_size = 10)
    return(p)
  }
  pp <- x$paths
  pp$funct <- factor(pp$funct, levels = x$functs)
  p <- ggplot2::ggplot() + .gg_elev(x$dtm) + .gg_shade(x$dtm, transp) +
    ggplot2::geom_sf(data = pp, ggplot2::aes(colour = funct),
                     linewidth = 0.6) +
    ggplot2::scale_colour_viridis_d(name = "cost\nfunction", option = "H")
  if (!is.null(x$paths.back)) {
    ppb <- x$paths.back
    ppb$funct <- factor(ppb$funct, levels = x$functs)
    p <- p + ggplot2::geom_sf(data = ppb, ggplot2::aes(colour = funct),
                              linewidth = 0.45, linetype = "dashed",
                              show.legend = FALSE)
  }
  p <- p +
    ggplot2::geom_sf(data = x$origin, shape = 20, size = 2.5, colour = "black") +
    ggplot2::geom_sf(data = x$destin, shape = 20, size = 2.5, colour = "red")
  .gg_finish(p, x$dtm, "LCPs compared across cost functions",
             paste(x$functs, collapse = ", "))
}

#' @rdname movecost-methods
#' @export
autoplot.movecost_comp <- function(object, ...) plot(object, ...)
