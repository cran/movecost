# Analytic validation of the computational core on synthetic terrains.
# These tests mirror the Python/SciPy validation performed during development
# (see inst/validation/), where the same checks all passed against an
# independent implementation.

test_that("flat terrain: accumulated time equals distance / Tobler speed", {
  dtm <- .flat_dtm()
  surf <- mc_surface(dtm, funct = "t", move = 16)
  origin <- .pt(dtm, 21, 21)

  acc <- mc_accum(surf, origin, time = "h")
  v0 <- .tobler_ms(0)

  # 10 cells due east = 100 m
  east <- terra::cellFromRowCol(dtm, 21, 31)
  expect_equal(terra::values(acc$accum)[east],
               (100 / v0) / 3600, tolerance = 1e-9)

  # 10 cells diagonal = sqrt(2) * 100 m
  diag <- terra::cellFromRowCol(dtm, 31, 31)
  expect_equal(terra::values(acc$accum)[diag],
               (sqrt(2) * 100 / v0) / 3600, tolerance = 1e-9)

  # knight-move neighbour = sqrt(5) * 10 m
  kn <- terra::cellFromRowCol(dtm, 23, 22)
  expect_equal(terra::values(acc$accum)[kn],
               (sqrt(5) * 10 / v0) / 3600, tolerance = 1e-9)
})

test_that("constant slope: up/down asymmetry matches Tobler analytically", {
  dtm <- .ramp_dtm(rise = 2)          # gradient 0.2 northwards
  surf <- mc_surface(dtm, funct = "t", move = 4)
  origin <- .pt(dtm, 21, 21)
  acc <- mc_accum(surf, origin, time = "h")

  up <- terra::cellFromRowCol(dtm, 11, 21)    # 10 cells uphill
  down <- terra::cellFromRowCol(dtm, 31, 21)  # 10 cells downhill
  t_up <- terra::values(acc$accum)[up]
  t_down <- terra::values(acc$accum)[down]

  expect_equal(t_up,   (100 / .tobler_ms(+0.2)) / 3600, tolerance = 1e-9)
  expect_equal(t_down, (100 / .tobler_ms(-0.2)) / 3600, tolerance = 1e-9)
  expect_gt(t_up, t_down)   # Tobler: climbing 20% is slower than descending it
})

test_that("topographic distance multiplies time by sqrt(1 + slope^2)", {
  dtm <- .ramp_dtm(rise = 2)
  origin <- .pt(dtm, 21, 21)
  up <- terra::cellFromRowCol(dtm, 11, 21)

  s1 <- mc_surface(dtm, funct = "t", move = 4)
  s2 <- mc_surface(dtm, funct = "t", move = 4, topo.dist = TRUE)
  a1 <- terra::values(mc_accum(s1, origin)$accum)[up]
  a2 <- terra::values(mc_accum(s2, origin)$accum)[up]
  expect_equal(a2, a1 * sqrt(1 + 0.2^2), tolerance = 1e-9)
})

test_that("cognitive slope rescales gradients as per Pingel (2013)", {
  dtm <- .ramp_dtm(rise = 2)
  origin <- .pt(dtm, 21, 21)
  up <- terra::cellFromRowCol(dtm, 11, 21)

  s <- mc_surface(dtm, funct = "t", move = 4, cogn.slp = TRUE)
  a <- terra::values(mc_accum(s, origin)$accum)[up]
  expect_equal(a, (100 / .tobler_ms(0.2 * 1.99)) / 3600, tolerance = 1e-9)
})

test_that("energy function (Herzog) accumulates exactly on flat terrain", {
  dtm <- .flat_dtm()
  origin <- .pt(dtm, 21, 21)
  east <- terra::cellFromRowCol(dtm, 21, 31)

  s <- mc_surface(dtm, funct = "hrz", move = 4)
  a <- terra::values(mc_accum(s, origin)$accum)[east]
  expect_equal(a, 1.64 * 100, tolerance = 1e-9)   # 1.64 J/(kg*m) * 100 m

  # terrain factor N = 2 doubles the cost
  s2 <- mc_surface(dtm, funct = "hrz", move = 4, N = 2)
  a2 <- terra::values(mc_accum(s2, origin)$accum)[east]
  expect_equal(a2, 2 * a, tolerance = 1e-9)
})

test_that("'ree' zero-cost flat edges are traversable, not impassable", {
  dtm <- .flat_dtm(21)
  surf <- mc_surface(dtm, funct = "ree", move = 8)
  acc <- mc_accum(surf, .pt(dtm, 11, 11))
  east <- terra::cellFromRowCol(dtm, 11, 21)
  expect_equal(terra::values(acc$accum)[east], 0, tolerance = 1e-12)
})

test_that("Pandolf-with-correction applies the downhill branch (3.0 fix)", {
  reg <- movecost:::.cf_registry()
  # check with a non-default terrain factor too: since 3.0, N is applied
  # once (Pandolf's eta), so pcf and p must coincide uphill for ANY N
  for (Nval in c(1, 2)) {
    p <- list(N = Nval, sl.crit = 10, W = 70, L = 0, V = 1.2)
    pcf <- reg$pcf$make(p)
    pp <- reg$p$make(p)
    # uphill: pcf coincides with plain Pandolf
    expect_equal(pcf(0.2), pp(0.2), tolerance = 1e-9)
    # downhill: the Yokota correction must change the cost
    expect_false(isTRUE(all.equal(pcf(-0.2), pp(-0.2))))
  }
})

test_that("Irmischer-Clarke variants are mutually coherent", {
  reg <- movecost:::.cf_registry()
  p <- list(N = 1, sl.crit = 10, W = 70, L = 0, V = 1.2)
  m_off <- reg$icmoffp$make(p)
  f_off <- reg$icfoffp$make(p)
  m_on <- reg$icmonp$make(p)
  # published relation: female variants are 0.95 * male variants
  s <- c(-0.3, -0.1, 0, 0.1, 0.3)
  expect_equal(f_off(s), 0.95 * m_off(s), tolerance = 1e-12)
  # off-path speed cannot exceed on-path speed at the same slope
  expect_true(all(m_off(s) <= m_on(s) + 1e-9))
})

test_that("Kondo-Seino and Marin Arroyo use the signed slope (3.0 fix)", {
  reg <- movecost:::.cf_registry()
  p <- list(N = 1, sl.crit = 10, W = 70, L = 0, V = 1.2)

  ks <- reg$ks$make(p)
  # at slope -0.2 (steep downhill) the gentle 1.5 decay must apply
  expect_equal(ks(-0.2), 5.1 * exp(-1.5 * abs(-0.2 + 0.07)), tolerance = 1e-12)
  # at slope +0.2 the steeper 2.25 decay applies
  expect_equal(ks(+0.2), 5.1 * exp(-2.25 * abs(0.2 + 0.07)), tolerance = 1e-12)

  ma <- reg$ma$make(p)
  # downhill: divisor 23; uphill: divisor 11
  expect_equal(1 / ma(-0.2), 0.6 * (20 / 23 + 1), tolerance = 1e-12)
  expect_equal(1 / ma(+0.2), 0.6 * (20 / 11 + 1), tolerance = 1e-12)
})

test_that("all 26 registered cost functions return positive values", {
  reg <- movecost:::.cf_registry()
  expect_length(reg, 26)
  p <- list(N = 1, sl.crit = 10, W = 70, L = 0, V = 1.2)
  slopes <- c(-0.4, -0.1, 0, 0.1, 0.4)
  for (nm in names(reg)) {
    f <- reg[[nm]]$make(p)
    out <- f(slopes)
    expect_true(all(out > 0), label = nm)
    # 'ree' legitimately returns Inf conductance on flat ground (zero cost);
    # all other functions must be finite everywhere
    if (nm != "ree") {
      expect_true(all(is.finite(out)), label = nm)
    } else {
      expect_true(all(is.finite(out[slopes != 0])), label = "ree (non-flat)")
      expect_true(is.infinite(f(0)), label = "ree (flat = zero cost)")
    }
  }
  # V = 0 variants (Tobler-derived speed)
  p0 <- p; p0$V <- 0
  for (nm in c("p", "pcf", "vl", "a")) {
    f <- reg[[nm]]$make(p0)
    out <- f(slopes)
    expect_true(all(is.finite(out)) && all(out > 0), label = paste0(nm, " V=0"))
  }
})
