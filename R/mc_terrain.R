#' Reference terrain factors (N) for use with a terrain raster
#'
#' Returns, as a data frame, the reference list of terrain factors reported in
#' the help page of \code{\link{mc_cost_functions}}. It is meant as a starting
#' point for building the lookup table (\code{terrain.table}) that
#' \code{\link{mc_surface}} uses to translate a categorical terrain (land-cover)
#' raster into terrain factors: pick the rows that describe the terrain types
#' of the study area, and pair them with the class codes (or labels) of the
#' raster (see Examples).\cr
#'
#' The values are taken from the data collected in Herzog, I. (2020), Spatial
#' Analysis Based on Cost Functions, in Gillings M., Haciguzeller P., Lock G.
#' (eds), "Archaeological Spatial Analysis. A Methodological Guide", Routledge:
#' New York, 340 (with previous references). They are split into two sets:
#' \code{type = "time"} for cost functions measuring walking time, and
#' \code{type = "other"} for functions measuring cost other than time (energy,
#' abstract cost). For ploughed fields, for which the source gives two values
#' (1.30 or 1.50), both are returned as separate rows.\cr
#'
#' \strong{A note of caution}: these factors were obtained in different
#' studies, under specific conditions, and for specific cost models. They are
#' convenient, documented defaults, not calibrated constants; when better
#' local evidence is available, it should be preferred.
#'
#' @param type character: \code{"time"} (default) for time-based cost
#'   functions, \code{"other"} for functions measuring cost other than time, or
#'   \code{"all"} for both sets.
#'
#' @return A data frame with columns \code{terrain} (description of the
#'   terrain type), \code{N} (terrain factor), and \code{type} (\code{"time"}
#'   or \code{"other"}).
#'
#' @examples
#' # reference values for time-based cost functions
#' tf <- mc_terrain_factors("time")
#' tf
#'
#' # build a lookup table for a raster whose cells are coded 1, 2, 3
#' # (1 = lawn grass, 2 = bog, 3 = rock)
#' tab <- data.frame(class = 1:3,
#'                   N = tf$N[match(c("Lawn grass", "Bog", "Rock"), tf$terrain)])
#' tab
#'
#' @seealso \code{\link{mc_surface}} (arguments \code{terrain} and
#'   \code{terrain.table}), \code{\link{mc_cost_functions}}
#' @export
mc_terrain_factors <- function(type = c("time", "other", "all")) {
  type <- match.arg(type)
  time_set <- data.frame(
    terrain = c("Blacktop roads, improved dirt paths, cement",
                "Lawn grass",
                "Loose beach sand",
                "Disturbed ground (former stone quarry)",
                "Horse riding path, flat trails and meadows",
                "Tall grassland (with thistle and nettles)",
                "Open space above the treeline (i.e., 2000 m asl)",
                "Bad trails, stony outcrops and river beds",
                "Off-paths",
                "Bog",
                "Off-path areas below the treeline (pastures, forests, heathland)",
                "Rock",
                "Swamp, water course"),
    N = c(1.00, 1.03, 1.19, 1.24, 1.25, 1.35, 1.50, 1.67, 1.67, 1.79, 2.00,
          2.50, 5.00),
    type = "time",
    stringsAsFactors = FALSE)
  other_set <- data.frame(
    terrain = c("Asphalt/blacktop",
                "Dirt road or grass",
                "Hard-surface road",
                "Light brush",
                "Ploughed field (lower value)",
                "Ploughed field (upper value)",
                "Heavy brush",
                "Hard-packed snow",
                "Swampy bog",
                "Sand dunes",
                "Loose sand"),
    N = c(1.00, 1.10, 1.20, 1.20, 1.30, 1.50, 1.50, 1.60, 1.80, 1.80, 2.10),
    type = "other",
    stringsAsFactors = FALSE)
  switch(type,
         time = time_set,
         other = other_set,
         all = rbind(time_set, other_set))
}

#' Internal: align a terrain raster to the DTM grid
#'
#' Coerces the input to a single-layer SpatRaster and, if its coordinate
#' system or grid (extent, resolution, origin) differs from the DTM's, aligns
#' it with \code{terra::project(method = "near")}. Nearest-neighbour is used
#' because it never invents values: it is the only safe choice for class codes
#' and it preserves the exact terrain factors of a numeric raster. A warning
#' is issued, because alignment changes the data and the user should know.
#'
#' @param dtm the DTM (SpatRaster).
#' @param terrain terrain raster (SpatRaster, RasterLayer, or file path).
#' @return a SpatRaster on the DTM grid.
#' @keywords internal
#' @noRd
.align_terrain <- function(dtm, terrain) {
  terrain <- .as_spatraster(terrain, "terrain")
  if (terra::nlyr(terrain) > 1) {
    warning("The terrain raster has ", terra::nlyr(terrain),
            " layers; only the first is used.", call. = FALSE)
    terrain <- terrain[[1]]
  }
  # a raster with no CRS is assumed to share the DTM's CRS
  if (!nzchar(terra::crs(terrain))) {
    warning("'terrain' has no coordinate reference system; ",
            "assuming the DTM's CRS.", call. = FALSE)
    terra::crs(terrain) <- terra::crs(dtm)
  }
  same_crs <- sf::st_crs(terra::crs(terrain)) == sf::st_crs(terra::crs(dtm))
  same_grid <- same_crs && isTRUE(terra::compareGeom(dtm, terrain, crs = FALSE,
                                                     stopOnError = FALSE,
                                                     messages = FALSE))
  if (!same_grid) {
    terrain <- terra::project(terrain, dtm, method = "near")
    warning("'terrain' did not match the DTM's ",
            if (same_crs) "grid" else "coordinate system and grid",
            "; it was aligned to the DTM with nearest-neighbour resampling.\n",
            "Align it beforehand (e.g. terra::project(terrain, dtm, method = 'near')) ",
            "to control this step yourself.", call. = FALSE)
  }
  terrain
}

#' Internal: per-cell terrain factors from a terrain raster
#'
#' Translates the (aligned) terrain raster into one terrain factor per DTM
#' cell: directly, if no lookup table is given, or through
#' \code{terrain.table} (class -> N) otherwise. Class matching uses the raster
#' codes when the table's classes are numeric, and the category labels of a
#' categorical (factor) raster when they are character. Cells without a
#' terrain value (NA in the raster, or a class missing from the table) take
#' the scalar \code{N}.
#'
#' @param dtm the DTM (SpatRaster).
#' @param terrain the terrain raster, already aligned to the DTM grid.
#' @param terrain.table NULL, or a data frame with columns \code{class} and
#'   \code{N} (or exactly two columns, read in that order).
#' @param N the scalar terrain factor used for cells with no terrain value.
#' @return a list with \code{values} (numeric, one per DTM cell; NA where the
#'   DTM is NA) and \code{n.filled} (number of valid DTM cells that took the
#'   scalar N).
#' @keywords internal
#' @noRd
.terrain_values <- function(dtm, terrain, terrain.table = NULL, N = 1) {
  codes <- terra::values(terrain, mat = FALSE)
  is_cat <- isTRUE(terra::is.factor(terrain)[1])

  if (is.null(terrain.table)) {
    if (is_cat) {
      stop("'terrain' is a categorical raster: supply 'terrain.table' to ",
           "translate its classes into terrain factors (see ",
           "mc_terrain_factors() for reference values).", call. = FALSE)
    }
    vals <- as.numeric(codes)
  } else {
    tab <- .check_terrain_table(terrain.table)
    if (is.character(tab$class)) {
      if (!is_cat) {
        stop("The classes in 'terrain.table' are labels (character), but ",
             "'terrain' is not a categorical raster; use the numeric cell ",
             "values as classes instead.", call. = FALSE)
      }
      lev <- terra::levels(terrain)[[1]]
      keys <- as.character(lev[[2]])[match(codes, lev[[1]])]
    } else {
      keys <- codes
    }
    idx <- match(keys, tab$class)
    unmatched <- unique(keys[!is.na(keys) & is.na(idx)])
    if (length(unmatched) > 0) {
      warning("Terrain class(es) not found in 'terrain.table': ",
              paste(utils::head(sort(unmatched), 20), collapse = ", "),
              if (length(unmatched) > 20) ", ..." else "",
              "; those cells take the scalar N = ", N, ".", call. = FALSE)
    }
    vals <- tab$N[idx]
  }

  bad <- !is.na(vals) & !(vals > 0)
  if (any(bad)) {
    stop("Terrain factors must be positive numbers (Inf is allowed and makes ",
         "a cell impassable); found ", sum(bad), " cell(s) with zero, ",
         "negative, or NaN values.", call. = FALSE)
  }

  elev <- terra::values(dtm, mat = FALSE)
  valid <- !is.na(elev)
  fill <- valid & is.na(vals)
  vals[fill] <- N
  vals[!valid] <- NA_real_
  list(values = vals, n.filled = sum(fill))
}

#' Internal: validate and normalise the terrain lookup table
#' @keywords internal
#' @noRd
.check_terrain_table <- function(tab) {
  if (!is.data.frame(tab)) {
    stop("'terrain.table' must be a data frame with columns 'class' and 'N'.",
         call. = FALSE)
  }
  if (all(c("class", "N") %in% names(tab))) {
    tab <- data.frame(class = tab$class, N = tab$N, stringsAsFactors = FALSE)
  } else if (ncol(tab) == 2L) {
    tab <- data.frame(class = tab[[1]], N = tab[[2]], stringsAsFactors = FALSE)
  } else {
    stop("'terrain.table' must have columns named 'class' and 'N' ",
         "(or exactly two columns: class, then N).", call. = FALSE)
  }
  if (is.factor(tab$class)) tab$class <- as.character(tab$class)
  if (!is.numeric(tab$N)) {
    stop("The 'N' column of 'terrain.table' must be numeric.", call. = FALSE)
  }
  if (anyNA(tab$N) || anyNA(tab$class)) {
    stop("'terrain.table' must not contain missing values.", call. = FALSE)
  }
  if (anyDuplicated(tab$class)) {
    stop("'terrain.table' lists the same class more than once: ",
         paste(unique(tab$class[duplicated(tab$class)]), collapse = ", "),
         call. = FALSE)
  }
  tab
}
