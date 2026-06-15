#' Download elevation data for a study area
#'
#' Acquires a Digital Terrain Model from online sources (via the \emph{elevatr}
#' package) for the area enclosed by the supplied polygon, and crops it to the
#' polygon's extent. This replaces the implicit download mechanism of
#' movecost <= 2.x's \code{studyplot} parameter (which remains available in
#' \code{\link{mc_surface}} and simply calls this function).\cr
#'
#' The zoom level \code{z} controls the resolution of the downloaded data
#' (default 9, a trade-off between resolution and download time). For details
#' on the data sources and the ground resolution per zoom level, see the
#' \emph{elevatr} documentation:
#' \url{https://cran.r-project.org/package=elevatr}.\cr
#'
#' \strong{Note}: the downloaded raster is returned in the CRS provided by the
#' source. Before using it with \code{\link{mc_surface}}, make sure it is in a
#' projected CRS (reproject with \code{terra::project()} if needed); a
#' geographic (lon/lat) raster is rejected by \code{mc_surface()} with an
#' informative error.
#'
#' @param studyplot polygon (sf, SpatVector, or legacy sp class) representing
#'   the study area for which elevation data are acquired.
#' @param z zoom level of the downloaded elevation data (0 to 15; 9 by
#'   default).
#'
#' @return A \code{SpatRaster}.
#'
#' @examples
#' \dontrun{
#' # elevation data for the area around Mt Etna (sample polygon in the package)
#' boundary <- mc_etna_boundary()
#' dtm <- mc_dtm(boundary, z = 9)
#' }
#'
#' @seealso \code{\link{mc_surface}}
#' @export
mc_dtm <- function(studyplot, z = 9) {
  if (!requireNamespace("elevatr", quietly = TRUE)) {
    stop("Downloading elevation data requires the 'elevatr' package.\n",
         "Install it with: install.packages('elevatr')", call. = FALSE)
  }
  sp_sf <- .as_sf(studyplot, "studyplot")
  elev <- elevatr::get_elev_raster(sp_sf, z = z, verbose = FALSE,
                                   override_size_check = TRUE)
  r <- .as_spatraster(elev)
  terra::crop(r, terra::vect(sp_sf))
}
