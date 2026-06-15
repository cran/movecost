#' Least-cost corridor between locations
#'
#' Computes the least-cost corridor between point locations, that is, the band
#' of terrain within which movement between the locations is comparatively
#' cheap. The corridor is built from accumulated cost surfaces produced by the
#' reused \code{\link{mc_surface}} graph. Two formulations are offered through
#' the \code{method} argument; they coincide for isotropic cost functions and
#' differ only when the cost function is anisotropic (direction-dependent, for
#' instance the Tobler hiking function).\cr
#'
#' \strong{method = "reach"} (default; the classic least-cost corridor):\cr
#' for every cell x the value is \eqn{cost(A \to x) + cost(B \to x)}, the sum of
#' the accumulated cost spreading outward from each location. This is the
#' definition used in the GIS literature and in movecost <= 2.x: "the cost of a
#' cell is the total cost to reach it from all the analysed locations" (Mitchell
#' A. 2012, The ESRI Guide to GIS Analysis, Vol. 3, Esri Press, 257-259). The
#' surface is symmetric in A and B, so it describes a direction-free corridor (a
#' band of ground that is cheap to reach from both endpoints), which matches the
#' usual conception of a corridor as bidirectional. Use this for the standard
#' corridor analysis, for reproducibility against movecost <= 2.x and other GIS
#' tools, and whenever the corridor is conceived as a two-way swath between
#' places.\cr
#'
#' \strong{method = "through"} (directional, near-optimal-route band):\cr
#' for every cell x the value is \eqn{cost(A \to x) + cost(x \to B)}, the cost of
#' a journey from A to B that is constrained to pass through x (the second term
#' is computed on the reversed graph). The minimum of this surface equals
#' exactly the cost of the least-cost path from A to B, and cells whose value is
#' within a chosen tolerance of that minimum form the band of near-optimal A-to-B
#' routes. Unlike "reach", this surface is direction-specific: the (A, B)
#' corridor differs from the (B, A) corridor under anisotropic costs. Use this
#' when the interest is explicitly in one-way journeys from A to B and in
#' identifying routes close to the optimal one in cost terms.\cr
#'
#' \strong{Why the two can differ}: with an anisotropic cost function the cost of
#' moving from B to x is not the cost of moving from x to B. The "reach" surface
#' adds two outward ("from") costs and is symmetric; the "through" surface adds
#' an outward and an inward cost and is directional. For symmetric (isotropic)
#' cost functions the two are identical. As a concrete illustration, on sloped
#' terrain with A low and B high the minimum of the "reach" surface corresponds
#' to the cheaper of the two one-way trips (typically the downhill direction),
#' whereas the minimum of the "through" surface is exactly the A-to-B (uphill)
#' least-cost-path cost.\cr
#'
#' \strong{Least-cost path overlay}: when exactly two locations are supplied
#' (\code{a} and \code{b}) and \code{lcp = TRUE} (default), the least-cost paths
#' from A to B and from B to A are also computed and returned, and are drawn over
#' the corridor by the plot method, as in movecost <= 2.x. The two paths differ
#' under anisotropic costs.\cr
#'
#' \strong{More than two locations}: if a set of locations is supplied via
#' \code{a} (with \code{b} left NULL), corridors are computed between all unique
#' pairs and summed, as in movecost <= 2.x; the individual pairwise corridors are
#' also returned. The accumulated surfaces are obtained from one multi-source
#' graph query (\code{"reach"}) or two (\code{"through"}), rather than two per
#' pair. No least-cost paths are overlaid in this case.
#'
#' @param surface a \code{movecost_surface} object created by
#'   \code{\link{mc_surface}}.
#' @param a first location, or (if \code{b} is NULL) a set of locations between
#'   whose unique pairs corridors are computed and summed (sf points,
#'   SpatVector, or legacy SpatialPointsDataFrame).
#' @param b second location (same accepted classes); NULL by default.
#' @param method \code{"reach"} (default) for the classic symmetric corridor
#'   (sum of two outward accumulated surfaces), or \code{"through"} for the
#'   directional A-to-B near-optimal-route band (see Details).
#' @param lcp TRUE (default) or FALSE: when exactly two locations are supplied,
#'   also compute the least-cost paths A-to-B and B-to-A and overlay them on the
#'   plot. Ignored when more than two locations are supplied.
#' @param rescale TRUE or FALSE (default): rescale the corridor values to the
#'   0-1 range (useful for comparing corridors computed with different cost
#'   functions).
#' @param time unit for time-based cost functions: "h" (hours, default) or
#'   "m" (minutes).
#'
#' @return An object of class \code{movecost_corridor}, a list with components
#' \itemize{
#'   \item \code{corridor}: the corridor surface (\code{SpatRaster});
#'   \item \code{corridors}: for more than two locations, a list of the
#'     individual pairwise corridor surfaces (NULL for the two-location case);
#'   \item \code{locations}: the input locations (\code{sf});
#'   \item \code{lcp.AtoB}, \code{lcp.BtoA}: the two least-cost paths (\code{sf}
#'     lines with \code{cost} and \code{length}) for the two-location case with
#'     \code{lcp = TRUE}, else NULL;
#'   \item \code{lcp.cost}: for \code{method = "through"} on two locations, the
#'     A-to-B least-cost-path cost (equal to the corridor minimum); NULL
#'     otherwise;
#'   \item \code{method}: the formulation used;
#'   \item metadata used by the plot method.
#' }
#' Plot it at any time with \code{plot()} (see
#' \code{\link{plot.movecost_corridor}}).
#'
#' @examples
#' dtm <- mc_volc()
#' destin <- mc_destin_loc()
#'
#' surf <- mc_surface(dtm, funct = "t", move = 8)
#'
#' # classic (symmetric) least-cost corridor, with the two LCPs overlaid
#' corr <- mc_corridor(surf, a = destin[1, ], b = destin[4, ])
#' plot(corr)
#'
#' # directional A-to-B near-optimal-route band
#' corr2 <- mc_corridor(surf, a = destin[1, ], b = destin[4, ], method = "through")
#' plot(corr2)
#'
#' @references
#' Mitchell, A. (2012). The ESRI Guide to GIS Analysis. Vol. 3: Modelling
#' Suitability, Movement, and Interaction. New York: Esri Press, 257-259.
#'
#' @seealso \code{\link{mc_surface}}, \code{\link{plot.movecost_corridor}},
#'   \code{\link{mc_paths}}
#' @export
mc_corridor <- function(surface, a, b = NULL, method = c("reach", "through"),
                        lcp = TRUE, rescale = FALSE, time = "h") {
  stopifnot(inherits(surface, "movecost_surface"))
  method <- match.arg(method)

  lcp.AtoB <- lcp.BtoA <- NULL
  lcp.cost <- NULL
  corridors <- NULL

  if (!is.null(b)) {
    # ---- two-location corridor ----
    A <- .points_to_cells(surface, a, "a")
    B <- .points_to_cells(surface, b, "b")
    if (length(A$cells) > 1L || length(B$cells) > 1L) {
      warning("Multiple features supplied in 'a' and/or 'b'; ",
              "using the first feature of each.", call. = FALSE)
    }
    a_cell <- A$cells[1]
    b_cell <- B$cells[1]

    from_a <- .accum_matrix(surface, a_cell, mode = "out")[1, ]
    second <- if (method == "reach") {
      .accum_matrix(surface, b_cell, mode = "out")[1, ]   # cost B -> x
    } else {
      .accum_matrix(surface, b_cell, mode = "in")[1, ]    # cost x -> B
    }
    vals <- from_a + second

    locations <- rbind(
      sf::st_sf(geometry = sf::st_geometry(A$sf[1, ])),
      sf::st_sf(geometry = sf::st_geometry(B$sf[1, ]))
    )

    if (method == "through") {
      lcp.cost <- .convert_units(min(vals[is.finite(vals)]), surface$funct, time)
    }
    if (isTRUE(lcp)) {
      lcp.AtoB <- .lcp_sf(surface, a_cell, b_cell, time)
      lcp.BtoA <- .lcp_sf(surface, b_cell, a_cell, time)
    }

  } else {
    # ---- multi-location: sum of all unique pairwise corridors ----
    A <- .points_to_cells(surface, a, "a")
    if (length(A$cells) < 2L) {
      stop("Provide either both 'a' and 'b', or at least two locations in 'a'.",
           call. = FALSE)
    }
    from_all <- .accum_matrix(surface, A$cells, mode = "out")
    to_all   <- if (method == "through") .accum_matrix(surface, A$cells, mode = "in") else NULL
    prs <- utils::combn(length(A$cells), 2)
    corridors <- vector("list", ncol(prs))
    vals <- 0
    for (k in seq_len(ncol(prs))) {
      i <- prs[1, k]; j <- prs[2, k]
      pij <- if (method == "reach") from_all[i, ] + from_all[j, ]
             else from_all[i, ] + to_all[j, ]
      cr <- terra::rast(surface$dtm)
      terra::values(cr) <- .convert_units(pij, surface$funct, time)
      cr <- terra::mask(cr, surface$dtm)
      names(cr) <- paste0("corridor_", i, "_", j)
      corridors[[k]] <- cr
      vals <- vals + pij
    }
    locations <- sf::st_sf(geometry = sf::st_geometry(A$sf))
  }

  vals <- .convert_units(vals, surface$funct, time)
  vals[!is.finite(vals)] <- NA_real_
  if (isTRUE(rescale)) {
    rng <- range(vals, na.rm = TRUE)
    vals <- (vals - rng[1]) / (rng[2] - rng[1])
    if (!is.null(lcp.cost)) lcp.cost <- (lcp.cost - rng[1]) / (rng[2] - rng[1])
  }

  corridor <- terra::rast(surface$dtm)
  terra::values(corridor) <- vals
  corridor <- terra::mask(corridor, surface$dtm)
  names(corridor) <- "corridor_cost"

  cf <- .cf_get(surface$funct)
  out <- list(corridor = corridor,
              corridors = corridors,
              locations = locations,
              lcp.AtoB = lcp.AtoB,
              lcp.BtoA = lcp.BtoA,
              lcp.cost = lcp.cost,
              method = method,
              rescaled = isTRUE(rescale),
              time = time,
              funct = surface$funct,
              funct.label = surface$funct.label,
              cost.unit = if (isTRUE(rescale)) "0-1 rescaled"
                          else if (cf$time_output) time else cf$unit,
              dtm = surface$dtm)
  class(out) <- "movecost_corridor"
  out
}

#' @rdname movecost-methods
#' @export
print.movecost_corridor <- function(x, ...) {
  cat("movecost least-cost corridor (method: ", x$method, ")\n", sep = "")
  cat("  cost function : ", x$funct, "\n", sep = "")
  cat("  locations     : ", nrow(x$locations), "\n", sep = "")
  if (!is.null(x$lcp.cost)) {
    cat("  A-to-B optimal path cost (corridor minimum): ",
        signif(x$lcp.cost, 4), " (", x$cost.unit, ")\n", sep = "")
  }
  if (!is.null(x$lcp.AtoB)) {
    cat("  least-cost paths A-to-B and B-to-A computed (overlaid by plot()).\n")
  }
  cat("Use plot() to visualise; mc_export() to write to disk.\n")
  invisible(x)
}
