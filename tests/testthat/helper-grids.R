# Synthetic terrains used across tests (projected CRS, 10 m resolution).

.flat_dtm <- function(n = 41, res = 10) {
  r <- terra::rast(nrows = n, ncols = n,
                   xmin = 0, xmax = n * res, ymin = 0, ymax = n * res,
                   crs = "EPSG:32633")
  terra::values(r) <- 0
  r
}

# rows rise by 'rise' metres per cell going northwards (constant gradient)
.ramp_dtm <- function(n = 41, res = 10, rise = 2) {
  r <- .flat_dtm(n, res)
  rowidx <- terra::rowFromCell(r, seq_len(terra::ncell(r)))
  terra::values(r) <- (n - rowidx) * rise   # row 1 = top = highest
  r
}

# sf point at the centre of cell (row, col)
.pt <- function(r, row, col) {
  cell <- terra::cellFromRowCol(r, row, col)
  xy <- terra::xyFromCell(r, cell)
  sf::st_as_sf(data.frame(id = 1, x = xy[1], y = xy[2]),
               coords = c("x", "y"), crs = sf::st_crs(terra::crs(r)))
}

# Tobler on-path walking speed in m/s, with the same km/h -> m/s constant
# (exact 1/3.6) used by the package since 3.1.1 (movecost <= 3.1.0: 0.278)
.tobler_ms <- function(slope) (6 * exp(-3.5 * abs(slope + 0.05))) / 3.6
