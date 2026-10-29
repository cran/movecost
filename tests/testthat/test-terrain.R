# Terrain factor supplied as a raster (mc_surface(terrain = , terrain.table = )).

# raster on the grid of 'dtm' filled with 'value'
.const_rast <- function(dtm, value) {
  r <- terra::rast(dtm)
  terra::values(r) <- value
  r
}

# codes of the cost functions that respond to N
.n_functs <- function() {
  cf <- mc_cost_functions()
  cf$code[cf$uses_N]
}

test_that("a constant terrain raster reproduces the scalar N exactly, for every N-aware function", {
  dtm <- .ramp_dtm(n = 11, rise = 2)
  k <- 1.79
  for (f in .n_functs()) {
    for (V in c(1.2, 0)) {           # V = 0: Tobler-derived speed (Pandolf, Van Leusen, Ardigo)
      s_scalar <- mc_surface(dtm, funct = f, move = 16, N = k, V = V)
      s_raster <- mc_surface(dtm, funct = f, move = 16, V = V,
                             terrain = .const_rast(dtm, k))
      expect_equal(igraph::E(s_raster$graph)$weight,
                   igraph::E(s_scalar$graph)$weight,
                   tolerance = 1e-12, info = paste0(f, ", V = ", V))
    }
  }
})

test_that("two terrain zones: accumulated time matches the half-edge analytic value", {
  dtm <- .flat_dtm()                         # 41 x 41, 10 m cells, flat
  cols <- terra::colFromCell(dtm, seq_len(terra::ncell(dtm)))
  terr <- .const_rast(dtm, ifelse(cols <= 20, 1, 2))
  surf <- mc_surface(dtm, funct = "t", move = 4, terrain = terr)
  acc <- mc_accum(surf, .pt(dtm, 21, 11), time = "h")

  # from column 11 to column 31 along a row: 9 edges on N = 1, the boundary
  # edge (columns 20-21) with mean N = 1.5, and 10 edges on N = 2
  target <- terra::cellFromRowCol(dtm, 21, 31)
  v0 <- .tobler_ms(0)
  expected <- (10 / v0) * (9 * 1 + 1.5 + 10 * 2) / 3600
  expect_equal(terra::values(acc$accum)[target], expected, tolerance = 1e-9)
})

test_that("cells with no terrain value take the scalar N (with a message)", {
  dtm <- .flat_dtm(n = 11)
  cols <- terra::colFromCell(dtm, seq_len(terra::ncell(dtm)))
  partial <- .const_rast(dtm, ifelse(cols <= 5, 2.5, NA))
  full <- .const_rast(dtm, ifelse(cols <= 5, 2.5, 1.67))

  expect_message(s_partial <- mc_surface(dtm, move = 8, N = 1.67,
                                         terrain = partial),
                 "take the scalar N = 1.67")
  s_full <- mc_surface(dtm, move = 8, terrain = full)
  expect_equal(igraph::E(s_partial$graph)$weight,
               igraph::E(s_full$graph)$weight, tolerance = 1e-12)
  expect_equal(s_partial$params$terrain.filled, 66L)
  expect_equal(range(terra::values(s_partial$terrain)), c(1.67, 2.5))
})

test_that("a categorical raster with a lookup table equals the numeric raster", {
  dtm <- .ramp_dtm(n = 11)
  rows <- terra::rowFromCell(dtm, seq_len(terra::ncell(dtm)))
  codes <- .const_rast(dtm, ifelse(rows <= 4, 1, ifelse(rows <= 8, 2, 3)))
  tab <- data.frame(class = 1:3, N = c(1.03, 1.79, 2.50))
  numeric_N <- .const_rast(dtm, tab$N[terra::values(codes, mat = FALSE)])

  s_num <- mc_surface(dtm, funct = "hrz", move = 16, terrain = numeric_N)
  s_code <- mc_surface(dtm, funct = "hrz", move = 16, terrain = codes,
                       terrain.table = tab)
  expect_equal(igraph::E(s_code$graph)$weight, igraph::E(s_num$graph)$weight)

  # unnamed two-column table: read as class, then N
  s_unnamed <- mc_surface(dtm, funct = "hrz", move = 16, terrain = codes,
                          terrain.table = data.frame(a = 1:3, b = tab$N))
  expect_equal(igraph::E(s_unnamed$graph)$weight, igraph::E(s_num$graph)$weight)

  # terra factor raster matched by its category labels
  fac <- codes
  levels(fac) <- data.frame(id = 1:3, cover = c("grass", "bog", "rock"))
  s_lab <- mc_surface(dtm, funct = "hrz", move = 16, terrain = fac,
                      terrain.table = data.frame(class = c("rock", "grass", "bog"),
                                                 N = c(2.50, 1.03, 1.79)))
  expect_equal(igraph::E(s_lab$graph)$weight, igraph::E(s_num$graph)$weight)
})

test_that("classes missing from the table warn and take the scalar N", {
  dtm <- .flat_dtm(n = 11)
  cols <- terra::colFromCell(dtm, seq_len(terra::ncell(dtm)))
  codes <- .const_rast(dtm, ifelse(cols <= 5, 1, 9))
  expect_warning(
    expect_message(s <- mc_surface(dtm, move = 8, N = 1.25, terrain = codes,
                                   terrain.table = data.frame(class = 1, N = 2)),
                   "scalar N"),
    "not found in 'terrain.table': 9")
  expect_setequal(unique(terra::values(s$terrain, mat = FALSE)), c(2, 1.25))
})

test_that("an infinite terrain factor makes a cell impassable", {
  dtm <- .flat_dtm(n = 11)
  blocked <- terra::cellFromRowCol(dtm, 6, 6)
  vals <- rep(1, terra::ncell(dtm))
  vals[blocked] <- Inf
  surf <- mc_surface(dtm, move = 16, terrain = .const_rast(dtm, vals))
  expect_equal(unname(igraph::degree(surf$graph, v = blocked)), 0)
  expect_s3_class(plot(surf, type = "terrain"), "ggplot")
  expect_match(surf$funct.label, "raster \\(1 to Inf\\)")
})

test_that("a misaligned or reprojected raster is aligned to the DTM, with a warning", {
  dtm <- .flat_dtm(n = 21)
  k <- 2
  ref <- mc_surface(dtm, move = 8, N = k)

  # coarser grid covering a larger extent
  coarse <- terra::rast(xmin = -15, xmax = 225, ymin = -15, ymax = 225,
                        resolution = 15, crs = "EPSG:32633")
  terra::values(coarse) <- k
  expect_warning(s_c <- mc_surface(dtm, move = 8, terrain = coarse),
                 "nearest-neighbour")
  expect_equal(igraph::E(s_c$graph)$weight, igraph::E(ref$graph)$weight)

  # different coordinate system (cells left uncovered fall back to N = k)
  utm34 <- terra::project(.const_rast(dtm, k), "EPSG:32634", method = "near")
  expect_warning(s_p <- suppressMessages(mc_surface(dtm, move = 8, N = k,
                                                    terrain = utm34)),
                 "coordinate system and grid")
  expect_equal(igraph::E(s_p$graph)$weight, igraph::E(ref$graph)$weight)

  # a raster with no CRS is assumed to be in the DTM's
  nocrs <- .const_rast(dtm, k)
  terra::crs(nocrs) <- ""
  expect_warning(s_n <- mc_surface(dtm, move = 8, terrain = nocrs),
                 "no coordinate reference system")
  expect_equal(igraph::E(s_n$graph)$weight, igraph::E(ref$graph)$weight)
})

test_that("functions that do not take a terrain factor ignore N and 'terrain', with a warning", {
  dtm <- .ramp_dtm(n = 11)
  ref <- mc_surface(dtm, funct = "tofp", move = 8)
  expect_warning(s <- mc_surface(dtm, funct = "tofp", move = 8,
                                 terrain = .const_rast(dtm, 2)),
                 "'terrain' is ignored")
  expect_null(s$terrain)
  expect_equal(igraph::E(s$graph)$weight, igraph::E(ref$graph)$weight)
  expect_warning(mc_surface(dtm, funct = "alb", move = 8, N = 1.5),
                 "N = 1.5 is ignored")
  expect_silent(mc_surface(dtm, funct = "alb", move = 8))

  # the warning gives the reason, and points off-path functions to their
  # on-path version (which takes N); Alberti's is not pointed anywhere
  expect_warning(mc_surface(dtm, funct = "tofp", move = 8, N = 1.5),
                 "already embeds off-path walking.*on-path version \\('t'\\)")
  expect_warning(mc_surface(dtm, funct = "icmoffp", move = 8,
                            terrain = .const_rast(dtm, 2)),
                 "To model terrain with a raster, use the on-path version \\('icmonp'\\)")
  expect_warning(mc_surface(dtm, funct = "icfoffp", move = 8, N = 2),
                 "To apply a terrain factor, use the on-path version \\('icfonp'\\)")
  w <- tryCatch(mc_surface(dtm, funct = "alb", move = 8,
                           terrain = .const_rast(dtm, 2)),
                warning = conditionMessage)
  expect_match(w, "walking pace of a flock")
  expect_false(grepl("on-path", w))
})

test_that("the functions that ignore N are exactly the off-path ones and Alberti's, and the on-path versions take N", {
  cf <- mc_cost_functions()
  expect_setequal(cf$code[!cf$uses_N], c("tofp", "icmoffp", "icfoffp", "alb"))
  expect_true(all(cf$uses_N[match(c("t", "icmonp", "icfonp"), cf$code)]))
})

test_that("invalid terrain inputs are rejected with informative errors", {
  dtm <- .flat_dtm(n = 11)
  codes <- .const_rast(dtm, 1)
  fac <- codes
  levels(fac) <- data.frame(id = 1, cover = "grass")

  expect_error(mc_surface(dtm, N = 0), "single positive number")
  expect_error(mc_surface(dtm, N = -1), "single positive number")
  expect_error(mc_surface(dtm, N = c(1, 2)), "single positive number")
  expect_error(mc_surface(dtm, N = NA), "single positive number")
  expect_error(mc_surface(dtm, terrain = .const_rast(dtm, 0)), "must be positive")
  expect_error(mc_surface(dtm, terrain = .const_rast(dtm, -1)), "must be positive")
  expect_error(mc_surface(dtm, terrain = fac), "supply 'terrain.table'")
  expect_error(mc_surface(dtm, terrain.table = data.frame(class = 1, N = 1)),
               "without a 'terrain' raster")
  expect_error(mc_surface(dtm, terrain = codes,
                          terrain.table = data.frame(class = "grass", N = 1)),
               "not a categorical raster")
  expect_error(mc_surface(dtm, terrain = codes,
                          terrain.table = data.frame(class = c(1, 1), N = 1:2)),
               "more than once")
  expect_error(mc_surface(dtm, terrain = codes,
                          terrain.table = data.frame(class = 1, N = NA_real_)),
               "missing values")
  expect_error(mc_surface(dtm, terrain = codes,
                          terrain.table = data.frame(a = 1, b = 1, c = 1)),
               "columns named 'class' and 'N'")
  expect_error(mc_surface(dtm, terrain = codes, terrain.table = list(1, 2)),
               "must be a data frame")
  expect_error(mc_surface(dtm, terrain = 5), "'terrain' must be a SpatRaster")
  expect_error(plot(mc_surface(dtm)), NA)
  expect_error(plot(mc_surface(dtm), type = "terrain"), "without a terrain raster")
})

test_that("terrain surfaces work across all analysis, I/O and display functions", {
  dtm <- mc_volc()
  start <- mc_volc_loc()
  destin <- mc_destin_loc()
  terr <- terra::classify(dtm, rbind(c(-Inf, 120, 1), c(120, Inf, 2)))
  tab <- data.frame(class = 1:2, N = c(1.03, 2.50))
  surf <- mc_surface(dtm, funct = "t", move = 8, terrain = terr,
                     terrain.table = tab)

  # the terrain makes travel slower than on the uniform N = 1 surface
  base <- mc_surface(dtm, funct = "t", move = 8)
  expect_gt(mc_paths(surf, start, destin[2, ])$paths$cost,
            mc_paths(base, start, destin[2, ])$paths$cost)

  expect_s3_class(mc_accum(surf, start), "movecost_accum")
  expect_s3_class(mc_paths(surf, start, destin, return.base = TRUE), "movecost_paths")
  expect_s3_class(mc_corridor(surf, a = destin[1, ], b = destin[4, ]), "movecost_corridor")
  expect_s3_class(mc_boundary(surf, start, limit = 0.1), "movecost_boundary")
  expect_s3_class(mc_alloc(surf, destin[1:3, ]), "movecost_alloc")
  expect_s3_class(mc_network(surf, destin[1:4, ]), "movecost_network")
  expect_s3_class(mc_rank(surf, start, destin[2, ], k = 2), "movecost_rank")

  # the label carried by every result names the terrain raster
  expect_match(mc_accum(surf, start)$funct.label, "N=raster \\(1.03 to 2.5\\)")
  expect_output(print(surf), "terrain factor: raster, N from 1.03 to 2.5")
  expect_output(print(base), "terrain factor: N = 1")
  expect_s3_class(plot(surf, type = "terrain"), "ggplot")
  expect_s3_class(plot(surf), "ggplot")

  # save / load keeps the terrain raster usable
  f <- tempfile(fileext = ".rds")
  mc_save(surf, f)
  back <- mc_load(f)
  expect_equal(terra::values(back$terrain), terra::values(surf$terrain))
  expect_s3_class(plot(back, type = "terrain"), "ggplot")

  # export writes the terrain factors
  d <- file.path(tempdir(), "mc_terrain_export")
  written <- suppressMessages(mc_export(surf, dir = d, overwrite = TRUE))
  expect_true(any(grepl("terrain_N_t\\.tif$", written)))
})

test_that("mc_comp aligns a terrain raster once and passes it to every function", {
  dtm <- mc_volc()
  start <- mc_volc_loc()
  destin <- mc_destin_loc()
  coarse <- terra::aggregate(.const_rast(dtm, 2), 2, fun = "mean", na.rm = TRUE)

  warns <- character(0)
  cmp <- withCallingHandlers(
    suppressMessages(mc_comp(dtm, start, destin[2, ], functs = c("t", "ks", "tofp"),
                             move = 8, N = 2, terrain = coarse)),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_equal(sum(grepl("nearest-neighbour", warns)), 1L)
  expect_equal(sum(grepl("'terrain' is ignored", warns)), 1L)

  ref <- mc_comp(dtm, start, destin[2, ], functs = c("t", "ks"), move = 8, N = 2)
  expect_equal(cmp$cost.table$cost[1:2], ref$cost.table$cost, tolerance = 1e-9)
})

test_that("mc_terrain_factors returns the documented reference values", {
  tt <- mc_terrain_factors("time")
  to <- mc_terrain_factors("other")
  ta <- mc_terrain_factors("all")
  expect_equal(nrow(tt), 13L)
  expect_equal(nrow(to), 11L)
  expect_equal(nrow(ta), 24L)
  expect_named(tt, c("terrain", "N", "type"))
  expect_true(all(tt$N >= 1 & tt$N <= 5))
  expect_equal(tt$N[tt$terrain == "Rock"], 2.5)
  expect_equal(to$N[to$terrain == "Loose sand"], 2.1)
  expect_error(mc_terrain_factors("speed"))
})
