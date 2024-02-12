#' Create barrier cost surface
#'
#' Creates a cost surface that incorporates barriers that inhibit movement in the landscape.\cr
#'
#' The resultant Barrier Cost Surface is produced by assessing which areas of the raster
#' coincide with the Spatial* or RasterLayer object as specified in the barrier argument.
#' The areas of the raster that coincide with the barrier are given a conductance value of 0
#' (default value, with all other areas given a Conductance value of 1 (default value).
#' The conductance value of 0 ensures that movement is inhibited within these areas.
#' Examples of use include rivers, altitudes, and taboo areas.
#' If a RasterLayer object is supplied in the barrier argument then all cells with a value NOT NA will be
#' used as the barrier. \cr
#'
#' @param raster RasterLayer (raster package). The Resolution, Extent, and Spatial Reference System of the provided RasterLayer is used when creating the resultant Barrier Cost Surface.
#' @param barrier Spatial* (sp package) or RasterLayer (raster package). Area within the landscape that movement is inhibited.
#' @param neighbours Number of directions used in the Least Cost Path calculation.
#' @param field Value assigned to cells that coincide with the barrier Spatial* or RasterLayer object. Default is numeric value 0.
#' @param background Value assigned to cells that do not coincide with the Spatial* or RasterLayer object. Default is numeric value 1.
#'
#' @return
#' TransitionLayer (gdistance package) numerically expressing the barriers to movement in the landscape.
#' The resultant TransitionLayer can be incorporated with other TransitionLayer through Raster calculations.
#'
#' @keywords internal
#' @importFrom raster projection extent
#' @importFrom Matrix Matrix
#' @importFrom methods new

create_barrier_cs  <- function (raster, barrier, neighbours = 16, field = 0, background = 1)
{
  if (!inherits(raster, "RasterLayer")) {
    stop("raster argument is invalid. Expecting a RasterLayer object")
  }
  if ((!inherits(barrier, "Spatial")) & (!inherits(barrier,
                                                   "RasterLayer"))) {
    stop("barrier argument is invalid. Expecting a Spatial* or RasterLayer object")
  }
  if (any(!neighbours %in% c(4, 8, 16, 32, 48)) & (!inherits(neighbours,
                                                             "matrix"))) {
    stop("neighbours argument is invalid. Expecting 4, 8, 16, 32, 48, or matrix object")
  }
  if ((!inherits(barrier, "RasterLayer")) & (field == "mask")) {
    stop("field agument is invalid. Expecting numeric value when Spatial* object supplied in barrier argument. See details for more.")
  }
  if (inherits(neighbours, "numeric")) {
    if (neighbours == 32) {
      neighbours_32 <-  matrix(c(NA, 1, 1, NA, 1, 1, NA, 1, NA, 1, NA, 1, NA, 1, 1, 1, 1, 1, 1, 1, 1, NA, NA, 1, 0, 1, NA, NA, 1, 1, 1, 1, 1,
               1, 1, 1, NA, 1, NA, 1, NA, 1, NA, 1, 1, NA, 1, 1, NA), nrow = 7, ncol = 7, byrow = TRUE)
      neighbours <- neighbours_32
    }
    else if (neighbours == 48) {
      neighbours_48 <- matrix(c(NA, 1, NA, 1, NA, 1, NA, 1, NA, 1, NA, 1, 1, NA, 1, 1, NA, 1, NA, 1, NA, 1, NA, 1, NA, 1, NA, 1, 1, 1, 1, 1,
                                1, 1, 1, 1, NA, NA, NA, 1, 0, 1, NA, NA, NA, 1, 1, 1, 1, 1, 1, 1, 1, 1, NA, 1, NA, 1, NA, 1, NA, 1, NA, 1, NA, 1, 1, NA, 1, 1, NA,
                                1, NA, 1, NA, 1, NA, 1, NA, 1, NA), nrow = 9, ncol = 9, byrow = TRUE)
      neighbours <- neighbours_48
    }
  }
  barrier_cs <- new("TransitionLayer", nrows = as.integer(nrow(raster)),
                    ncols = as.integer(ncol(raster)), extent = extent(raster),
                    crs = projection(raster, asText = FALSE), transitionMatrix = Matrix(0,
                                                                                        ncell(raster), ncell(raster)), transitionCells = 1:ncell(raster))
  adj <- raster::adjacent(raster, cells = 1:raster::ncell(raster),
                          pairs = TRUE, directions = neighbours)
  if (inherits(barrier, "SpatialPoints")) {
    barrier_cells <- raster::cellFromXY(object = raster,
                                        xy = barrier)
  }
  else if (inherits(barrier, "SpatialLines")) {
    barrier_cells <- unlist(raster::cellFromLine(object = raster,
                                                 lns = barrier))
  }
  else if (inherits(barrier, "SpatialPolygons")) {
    barrier_cells <- unlist(raster::cellFromPolygon(object = raster,
                                                    p = barrier))
  }
  else if (inherits(barrier, "RasterLayer")) {
    barrier_cells <- which(!is.na(raster::getValues(barrier)))
    if (field == "mask") {
      field <- raster::extract(x = barrier, adj[adj[, 2] %in%
                                                  barrier_cells, ][, 2])
    }
  }
  barrier_cs[adj[adj[, 2] %in% barrier_cells, ]] <- field
  barrier_cs[adj[!adj[, 2] %in% barrier_cells, ]] <- background
  barrier_cs@transitionMatrix <- Matrix::drop0(barrier_cs@transitionMatrix)
  return(barrier_cs)
}
