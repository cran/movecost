# Interface behaviour: classes, plots, I/O, utilities, sample data.

test_that("plot methods return ggplot objects without re-computation", {
  dtm <- .flat_dtm(21)
  surf <- mc_surface(dtm, funct = "t", move = 8)
  o <- .pt(dtm, 11, 11); d <- .pt(dtm, 5, 17)

  acc <- mc_accum(surf, o)
  lcp <- mc_paths(surf, o, d)
  expect_s3_class(plot(surf), "ggplot")
  expect_s3_class(plot(acc), "ggplot")
  expect_s3_class(plot(acc, type = "contour"), "ggplot")
  expect_s3_class(plot(lcp), "ggplot")
  expect_s3_class(ggplot2::autoplot(acc), "ggplot")

  # a stored plot can be customised afterwards, ggplot2-style
  p <- plot(acc) + ggplot2::ggtitle("custom title")
  expect_s3_class(p, "ggplot")
})

test_that("hms formatting and unit conversion are correct", {
  expect_equal(movecost:::.fmt_hms(1.5, "h"), "01:30:00")
  expect_equal(movecost:::.fmt_hms(90, "m"), "01:30:00")
  expect_equal(movecost:::.fmt_hms(c(NA, 2), "h")[1], NA_character_)
  expect_equal(movecost:::.time_convert(7200, "h"), 2)
  expect_equal(movecost:::.time_convert(120, "m"), 2)
})

test_that("time-based results carry cost_hms; minutes unit works", {
  dtm <- .flat_dtm(21)
  surf <- mc_surface(dtm, funct = "t", move = 8)
  o <- .pt(dtm, 11, 11); d <- .pt(dtm, 11, 21)

  lcp_h <- mc_paths(surf, o, d, time = "h")
  lcp_m <- mc_paths(surf, o, d, time = "m")
  expect_true("cost_hms" %in% names(lcp_h$paths))
  expect_equal(lcp_m$paths$cost[1], lcp_h$paths$cost[1] * 60, tolerance = 1e-9)
  expect_equal(lcp_h$paths$cost_hms[1], lcp_m$paths$cost_hms[1])
})

test_that("geographic (lon/lat) DTMs are rejected with a clear error", {
  r <- terra::rast(nrows = 10, ncols = 10, xmin = 14, xmax = 15,
                   ymin = 35, ymax = 36, crs = "EPSG:4326")
  terra::values(r) <- 0
  expect_error(mc_surface(r), "projected")
})

test_that("locations outside the DTM or on NoData cells are rejected", {
  dtm <- .flat_dtm(21)
  surf <- mc_surface(dtm, funct = "t", move = 8)
  outside <- sf::st_as_sf(data.frame(x = 1e6, y = 1e6),
                          coords = c("x", "y"),
                          crs = sf::st_crs(terra::crs(dtm)))
  expect_error(mc_accum(surf, outside), "outside")
})

test_that("SpatVector input is accepted (input coercion)", {
  dtm <- .flat_dtm(21)
  surf <- mc_surface(dtm, funct = "t", move = 8)
  o <- terra::vect(.pt(dtm, 11, 11))
  expect_s3_class(mc_accum(surf, o), "movecost_accum")
})

test_that("mc_save / mc_load round-trip preserves results", {
  dtm <- .flat_dtm(21)
  surf <- mc_surface(dtm, funct = "t", move = 8)
  acc <- mc_accum(surf, .pt(dtm, 11, 11))

  f <- tempfile(fileext = ".rds")
  mc_save(acc, f)
  acc2 <- mc_load(f)
  expect_s3_class(acc2, "movecost_accum")
  expect_equal(terra::values(acc2$accum), terra::values(acc$accum))
  unlink(f)
})

test_that("mc_export writes GeoTIFF and GeoPackage files", {
  dtm <- .flat_dtm(21)
  surf <- mc_surface(dtm, funct = "t", move = 8)
  lcp <- mc_paths(surf, .pt(dtm, 11, 11), .pt(dtm, 5, 17))

  dir <- file.path(tempdir(), "mc_export_test")
  files <- suppressMessages(mc_export(lcp, dir = dir, overwrite = TRUE))
  expect_true(all(file.exists(files)))
  expect_true(any(grepl("LCPs_t\\.gpkg$", files)))
  unlink(dir, recursive = TRUE)
})

test_that("mc_comp compares functions and returns a usable cost table", {
  dtm <- .ramp_dtm(n = 25, rise = 2)
  o <- .pt(dtm, 21, 5); d <- .pt(dtm, 5, 21)
  cmp <- mc_comp(dtm, o, d, functs = c("t", "ks"), move = 8)
  expect_equal(sort(unique(cmp$cost.table$funct)), c("ks", "t"))
  expect_s3_class(plot(cmp), "ggplot")
  expect_s3_class(plot(cmp, type = "chart"), "ggplot")
})

test_that("mc_cost_functions() documents all functions coherently", {
  tab <- mc_cost_functions()
  expect_equal(nrow(tab), 26)
  expect_true(all(c("t", "tofp", "mp", "icmonp", "icmoffp", "icfonp",
                    "icfoffp", "ug", "ma", "alb", "gkrs", "r", "ks", "trp",
                    "wcs", "ree", "b", "e", "p", "pcf", "m", "hrz", "vl",
                    "ls", "a", "h") %in% tab$code))
  expect_error(mc_surface(.flat_dtm(11), funct = "nope"), "Unknown cost function")
})

test_that("bundled sample data load and interoperate", {
  skip_if_not(file.exists(system.file("extdata", "volc.tif",
                                      package = "movecost")))
  dtm <- mc_volc()
  expect_s4_class(dtm, "SpatRaster")
  o <- mc_volc_loc(); d <- mc_destin_loc()
  expect_s3_class(o, "sf")
  expect_equal(nrow(d), 9)

  surf <- mc_surface(dtm, funct = "t", move = 8)
  lcp <- mc_paths(surf, o, d[2, ])
  expect_true(is.finite(lcp$paths$cost[1]))
  expect_gt(as.numeric(lcp$paths$length[1]), 0)
})

test_that("boundary polygons carry exact area and perimeter", {
  dtm <- .flat_dtm(41)
  surf <- mc_surface(dtm, funct = "t", move = 8)
  o <- .pt(dtm, 21, 21)
  bnd <- mc_boundary(surf, o, limit = 2, time = "m")
  expect_s3_class(bnd$boundaries, "sf")
  expect_true(all(c("area", "perimeter") %in% names(bnd$boundaries)))
  expect_gt(as.numeric(bnd$boundaries$area[1]), 0)
})
