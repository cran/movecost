#' Export movecost results to disk
#'
#' Writes the spatial components of any movecost result object to disk:
#' rasters as GeoTIFF, vector layers as GeoPackage (which, unlike the ESRI
#' shapefiles produced by movecost <= 2.x, preserves full attribute names and
#' supports mixed geometry without sidecar files). File names carry the cost
#' function code as a suffix, as in movecost <= 2.x.
#'
#' @param x a movecost result object (any of \code{movecost_surface},
#'   \code{movecost_accum}, \code{movecost_paths}, \code{movecost_corridor},
#'   \code{movecost_boundary}, \code{movecost_alloc},
#'   \code{movecost_network}, \code{movecost_rank}, \code{movecost_comp}).
#' @param dir output directory (default: current working directory).
#' @param prefix optional file-name prefix (default: the class-specific name).
#' @param overwrite TRUE or FALSE (default): overwrite existing files.
#'
#' @return Invisibly, a character vector with the paths of the written files.
#'
#' @examples
#' \dontrun{
#' dtm <- mc_volc()
#' start <- mc_volc_loc(); destin <- mc_destin_loc()
#' surf <- mc_surface(dtm, funct = "t", move = 8)
#' acc <- mc_accum(surf, start, breaks = 0.05)
#' mc_export(acc, dir = tempdir())
#' }
#'
#' @export
mc_export <- function(x, dir = ".", prefix = NULL, overwrite = FALSE) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  suffix <- if (!is.null(x$funct)) paste0("_", x$funct) else ""
  written <- character(0)

  wr <- function(r, name) {
    f <- file.path(dir, paste0(name, suffix, ".tif"))
    terra::writeRaster(r, f, overwrite = overwrite)
    written <<- c(written, f)
  }
  wv <- function(v, name) {
    f <- file.path(dir, paste0(name, suffix, ".gpkg"))
    sf::st_write(v, f, delete_dsn = overwrite, quiet = TRUE)
    written <<- c(written, f)
  }

  if (inherits(x, "movecost_surface")) {
    wr(x$cost.raster, prefix %||% "cost_surface")
    wr(x$dtm, "dtm")
  } else if (inherits(x, "movecost_accum")) {
    wr(x$accum, prefix %||% "accum_cost_surface")
    if (!is.null(x$isolines)) wv(x$isolines, "isolines")
  } else if (inherits(x, "movecost_paths")) {
    wv(x$paths, prefix %||% "LCPs")
    if (!is.null(x$paths.back)) wv(x$paths.back, "LCPs_back")
    wv(x$destin, "dest_loc_w_cost")
  } else if (inherits(x, "movecost_corridor")) {
    wr(x$corridor, prefix %||% "corridor")
    if (!is.null(x$lcp.AtoB)) wv(x$lcp.AtoB, "lcp_a_to_b")
    if (!is.null(x$lcp.BtoA)) wv(x$lcp.BtoA, "lcp_b_to_a")
  } else if (inherits(x, "movecost_boundary")) {
    wv(x$boundaries, prefix %||% "boundaries")
    wr(x$accum, "accum_cost_surface")
  } else if (inherits(x, "movecost_alloc")) {
    wr(x$alloc, prefix %||% "allocation")
    wv(x$zones, "allocation_zones")
    if (!is.null(x$isolines)) wv(x$isolines, "isolines")
  } else if (inherits(x, "movecost_network")) {
    wv(x$paths, prefix %||% "network_LCPs")
    if (!is.null(x$density.count)) {
      wr(x$density.count, "LCPs_density_count")
      wr(x$density.perc, "LCPs_density_perc")
    }
    f <- file.path(dir, paste0("cost_matrix", suffix, ".csv"))
    utils::write.csv(x$cost.matrix, f)
    written <- c(written, f)
  } else if (inherits(x, "movecost_rank")) {
    wv(x$paths, prefix %||% "ranked_LCPs")
  } else if (inherits(x, "movecost_comp")) {
    wv(x$paths, prefix %||% "compared_LCPs")
    if (!is.null(x$paths.back)) wv(x$paths.back, "compared_LCPs_back")
    f <- file.path(dir, "cost_table.csv")
    utils::write.csv(x$cost.table, f, row.names = FALSE)
    written <- c(written, f)
  } else {
    stop("Unsupported object of class: ", paste(class(x), collapse = "/"),
         call. = FALSE)
  }
  message("Written:\n", paste(" -", written, collapse = "\n"))
  invisible(written)
}

#' Save and load movecost objects across R sessions
#'
#' movecost result objects contain terra \code{SpatRaster} components, which
#' hold pointers to external memory and therefore cannot be stored with plain
#' \code{saveRDS()}/\code{readRDS()}. \code{mc_save()} wraps all raster
#' components (via \code{terra::wrap()}) before serialising; \code{mc_load()}
#' restores them. igraph components are serialised natively.
#'
#' @param x a movecost result object.
#' @param file path of the file to write/read (conventionally with extension
#'   \code{.rds}).
#'
#' @return \code{mc_save()} returns \code{file} invisibly; \code{mc_load()}
#'   returns the restored object.
#'
#' @examples
#' \dontrun{
#' surf <- mc_surface(mc_volc(), funct = "t", move = 8)
#' f <- tempfile(fileext = ".rds")
#' mc_save(surf, f)
#' surf2 <- mc_load(f)
#' }
#'
#' @export
mc_save <- function(x, file) {
  # recurse so that SpatRaster components are wrapped at any depth (e.g. the
  # list of pairwise corridors in a multi-location movecost_corridor object)
  y <- .mc_wrap(unclass(x))
  attr(y, "movecost_class") <- class(x)
  saveRDS(y, file)
  invisible(file)
}

#' @rdname mc_save
#' @export
mc_load <- function(file) {
  y <- readRDS(file)
  cl <- attr(y, "movecost_class")
  out <- .mc_unwrap(y)
  attr(out, "movecost_class") <- NULL
  class(out) <- cl
  out
}

#' Internal: recursively (un)wrap SpatRaster components for serialisation
#' @keywords internal
#' @noRd
.mc_wrap <- function(el) {
  if (inherits(el, "SpatRaster")) return(terra::wrap(el))
  # recurse only into plain (unclassed) lists, never into sf/data.frame, which
  # are also lists but must keep their structure and class
  if (is.list(el) && !is.object(el)) return(lapply(el, .mc_wrap))
  el
}

#' @keywords internal
#' @noRd
.mc_unwrap <- function(el) {
  if (inherits(el, "PackedSpatRaster")) return(terra::unwrap(el))
  if (is.list(el) && !is.object(el)) return(lapply(el, .mc_unwrap))
  el
}
