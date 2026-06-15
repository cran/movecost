# Internal consistency laws on rugged synthetic terrain.

.rugged_dtm <- function(n = 35, res = 10, seed = 42) {
  set.seed(seed)
  r <- terra::rast(nrows = n, ncols = n,
                   xmin = 0, xmax = n * res, ymin = 0, ymax = n * res,
                   crs = "EPSG:32633")
  terra::values(r) <- stats::rnorm(n * n, 100, 25)
  terra::focal(r, w = 5, fun = "mean", na.policy = "omit", fillvalue = 100)
}

test_that("LCP cost equals accumulated cost at destination", {
  dtm <- .rugged_dtm()
  surf <- mc_surface(dtm, funct = "t", move = 16)
  o <- .pt(dtm, 4, 4); d <- .pt(dtm, 32, 32)

  acc <- mc_accum(surf, o, time = "h")
  lcp <- mc_paths(surf, o, d, time = "h")
  dcell <- terra::cellFromXY(dtm, sf::st_coordinates(d))
  expect_equal(lcp$paths$cost[1], terra::values(acc$accum)[dcell],
               tolerance = 1e-9)
  # endpoints of the path coincide with origin and destination cell centres
  coords <- sf::st_coordinates(lcp$paths)[, 1:2]
  expect_equal(unname(coords[1, ]), unname(sf::st_coordinates(o)[1, ]),
               tolerance = 1e-6)
  expect_equal(unname(coords[nrow(coords), ]), unname(sf::st_coordinates(d)[1, ]),
               tolerance = 1e-6)
})

test_that("corridor 'through' minimum equals the A-to-B optimal path cost", {
  dtm <- .rugged_dtm()
  surf <- mc_surface(dtm, funct = "t", move = 16)
  o <- .pt(dtm, 4, 4); d <- .pt(dtm, 32, 32)

  corr <- mc_corridor(surf, a = o, b = d, method = "through", time = "h")
  lcp <- mc_paths(surf, o, d, time = "h")
  expect_equal(corr$lcp.cost, lcp$paths$cost[1], tolerance = 1e-9)
  expect_equal(min(terra::values(corr$corridor), na.rm = TRUE),
               lcp$paths$cost[1], tolerance = 1e-9)
})

test_that("corridor 'reach' (default) is the symmetric sum of two outward surfaces", {
  dtm <- .ramp_dtm(n = 21, rise = 2)            # anisotropic terrain
  surf <- mc_surface(dtm, funct = "t", move = 8)
  o <- .pt(dtm, 18, 10); d <- .pt(dtm, 4, 10)   # A low, B high

  corr_ab <- mc_corridor(surf, a = o, b = d)              # default reach
  corr_ba <- mc_corridor(surf, a = d, b = o)
  # symmetry: swapping a and b leaves the reach surface unchanged
  expect_equal(terra::values(corr_ab$corridor),
               terra::values(corr_ba$corridor), tolerance = 1e-9)
  # the reach surface equals cost(A->x) + cost(B->x) by construction
  fa <- movecost:::.accum_matrix(surf, terra::cellFromXY(dtm, sf::st_coordinates(o)), "out")[1, ]
  fb <- movecost:::.accum_matrix(surf, terra::cellFromXY(dtm, sf::st_coordinates(d)), "out")[1, ]
  expected <- (fa + fb) / 3600  # to hours
  expect_equal(terra::values(corr_ab$corridor, mat = FALSE),
               ifelse(is.finite(expected), expected, NA_real_), tolerance = 1e-9)
  # reach default carries no lcp.cost but does carry the two overlaid LCPs
  expect_null(corr_ab$lcp.cost)
  expect_s3_class(corr_ab$lcp.AtoB, "sf")
  expect_s3_class(corr_ab$lcp.BtoA, "sf")
})

test_that("allocation assigns each origin its own cell and is consistent", {
  dtm <- .rugged_dtm()
  surf <- mc_surface(dtm, funct = "t", move = 8)
  origins <- rbind(.pt(dtm, 3, 3), .pt(dtm, 3, 33), .pt(dtm, 33, 17))

  al <- mc_alloc(surf, origins)
  cells <- terra::cellFromXY(dtm, sf::st_coordinates(origins))
  expect_equal(terra::values(al$alloc)[cells], c(1, 2, 3))
  expect_s4_class(al$alloc, "SpatRaster")
  expect_s3_class(al$zones, "sf")
})

test_that("barriers block or penalise movement as requested", {
  dtm <- .rugged_dtm()
  o <- .pt(dtm, 4, 4); d <- .pt(dtm, 32, 32)
  # vertical wall through the middle, leaving a gap at the bottom
  wall <- sf::st_sf(geometry = sf::st_sfc(
    sf::st_linestring(rbind(c(175, 350), c(175, 60))),
    crs = sf::st_crs(terra::crs(dtm))))

  free <- mc_paths(mc_surface(dtm, funct = "t", move = 8), o, d)
  blocked <- mc_paths(mc_surface(dtm, funct = "t", move = 8,
                                 barrier = wall, field = 0), o, d)
  penalised <- mc_paths(mc_surface(dtm, funct = "t", move = 8,
                                   barrier = wall, field = 0.01), o, d)

  expect_gt(blocked$paths$cost[1], free$paths$cost[1])
  expect_gte(blocked$paths$cost[1], penalised$paths$cost[1] - 1e-9)
  expect_gt(penalised$paths$cost[1], free$paths$cost[1] - 1e-9)

  # the blocked path must not touch barrier cells
  bcells <- movecost:::.cells_touched(dtm, wall)
  pcells <- terra::cellFromXY(dtm, sf::st_coordinates(blocked$paths)[, 1:2])
  expect_false(any(pcells %in% bcells))
})

test_that("NoData cells are never traversed (irregular DTMs need no workaround)", {
  dtm <- .rugged_dtm()
  # carve a NoData 'sea inlet'
  v <- terra::values(dtm, mat = FALSE)
  bad <- which(terra::rowFromCell(dtm, seq_along(v)) %in% 16:20 &
                 terra::colFromCell(dtm, seq_along(v)) <= 20)
  v[bad] <- NA
  terra::values(dtm) <- v

  surf <- mc_surface(dtm, funct = "t", move = 16)
  o <- .pt(dtm, 4, 4); d <- .pt(dtm, 32, 32)
  lcp <- mc_paths(surf, o, d)
  pcells <- terra::cellFromXY(dtm, sf::st_coordinates(lcp$paths)[, 1:2])
  expect_false(any(pcells %in% bad))

  acc <- mc_accum(surf, o)
  expect_true(all(is.na(terra::values(acc$accum)[bad])))
})

test_that("ranked paths have non-decreasing costs and start optimal", {
  dtm <- .rugged_dtm()
  surf <- mc_surface(dtm, funct = "t", move = 16)
  o <- .pt(dtm, 4, 4); d <- .pt(dtm, 32, 32)

  rk <- mc_rank(surf, o, d, k = 3)
  lcp <- mc_paths(surf, o, d)
  expect_equal(rk$paths$cost[1], lcp$paths$cost[1], tolerance = 1e-9)
  expect_true(all(diff(rk$paths$cost) >= -1e-9))
  expect_equal(rk$paths$rank, 1:3)
})

test_that("network cost matrix is asymmetric on sloped terrain", {
  dtm <- .ramp_dtm(n = 25, rise = 2)
  surf <- mc_surface(dtm, funct = "t", move = 8)
  nodes <- rbind(.pt(dtm, 5, 13), .pt(dtm, 21, 13))   # one high, one low

  nw <- mc_network(surf, nodes, type = "allpairs")
  cm <- nw$cost.matrix
  expect_gt(cm[2, 1], cm[1, 2])   # uphill > downhill (Tobler)
  expect_equal(nrow(nw$paths), 2L)

  nwd <- mc_network(surf, nodes, type = "neigh", density = TRUE)
  expect_s4_class(nwd$density.count, "SpatRaster")
  expect_lte(max(terra::values(nwd$density.perc), na.rm = TRUE), 100 + 1e-9)
})

test_that("multi-origin accumulation takes the per-cell minimum", {
  dtm <- .rugged_dtm()
  surf <- mc_surface(dtm, funct = "t", move = 8)
  o1 <- .pt(dtm, 3, 3); o2 <- .pt(dtm, 33, 33)

  a1 <- mc_accum(surf, o1)
  a2 <- mc_accum(surf, o2)
  both <- mc_accum(surf, rbind(o1, o2))
  expected <- pmin(terra::values(a1$accum), terra::values(a2$accum))
  expect_equal(terra::values(both$accum), expected, tolerance = 1e-9)
})
