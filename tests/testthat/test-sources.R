# Cost functions checked against their published equations (movecost 3.1.1).
# Expected values are written out independently from the source equations
# (not from the package registry), with the slope signed in the direction of
# travel (negative downhill), as in every source cited below.

.reg <- function() movecost:::.cf_registry()
.par <- function(N = 1, W = 70, L = 0, V = 1.2) {
  list(N = N, sl.crit = 10, W = W, L = L, V = V)
}

test_that("Irmischer-Clarke (2018, eqs. 2-5) use the signed percent slope", {
  reg <- .reg(); p <- .par()
  ic <- function(s, k, a, f = 1) f * (0.11 + a * exp(-(100 * s + k)^2 / 1800)) * 3.6
  s <- c(-0.30, -0.15, -0.05, -0.02, 0, 0.05, 0.30)
  expect_equal(reg$icmonp$make(p)(s), ic(s, 5, 1), tolerance = 1e-12)
  expect_equal(reg$icfonp$make(p)(s), ic(s, 5, 1, 0.95), tolerance = 1e-12)
  expect_equal(reg$icmoffp$make(p)(s), ic(s, 2, 0.67), tolerance = 1e-12)
  expect_equal(reg$icfoffp$make(p)(s), ic(s, 2, 0.67, 0.95), tolerance = 1e-12)
  # peak speed at -5% (on-path) and -2% (off-path), not on level ground
  grid <- seq(-0.3, 0.3, by = 0.001)
  expect_equal(grid[which.max(reg$icmonp$make(p)(grid))], -0.05, tolerance = 1e-9)
  expect_equal(grid[which.max(reg$icmoffp$make(p)(grid))], -0.02, tolerance = 1e-9)
})

test_that("Rees (2004, eq. 1) uses the signed slope", {
  f <- .reg()$r$make(.par())
  s <- c(-0.2, 0.2)
  expect_equal(f(s), 3.6 / (0.75 + 0.09 * s + 14.6 * s^2), tolerance = 1e-12)
  expect_gt(f(-0.2), f(0.2))   # descending is faster than ascending
})

test_that("Tripcevich (2008) peaks on a -4.127 degree descent", {
  f <- .reg()$trp$make(.par())
  deg <- c(-20, -4.127, 0, 20)
  s <- tan(deg * pi / 180)
  expect_equal(f(s), 4.028 * 46^2 / ((deg + 4.127)^2 + 46^2), tolerance = 1e-12)
  expect_equal(f(tan(-4.127 * pi / 180)), 4.028, tolerance = 1e-12)
})

test_that("Minetti et al. (2002) use the signed gradient, held at -0.45 below it", {
  f <- .reg()$m$make(.par())
  cw <- function(i) 280.5 * i^5 - 58.7 * i^4 - 76.8 * i^3 + 51.9 * i^2 + 19.6 * i + 2.5
  i <- c(-0.45, -0.10, 0, 0.10, 0.45)
  expect_equal(1 / f(i), cw(i), tolerance = 1e-12)
  expect_equal(1 / f(-0.10), 1.1271, tolerance = 1e-4)   # descent cheaper than level
  expect_equal(f(-0.60), f(-0.45), tolerance = 1e-12)    # floor
  expect_true(all(f(seq(-1, 1, 0.01)) > 0))
})

test_that("mc_surface reports edges beyond Minetti's measured range", {
  steep <- .ramp_dtm(n = 11, res = 10, rise = 6)          # 60% gradient
  expect_message(mc_surface(steep, funct = "m", move = 4), "steeper")
  expect_no_message(mc_surface(steep, funct = "hrz", move = 4))
})

test_that("Herzog (2016) polynomial uses the signed slope", {
  f <- .reg()$hrz$make(.par())
  h <- function(s) 1337.8 * s^6 + 278.19 * s^5 - 517.39 * s^4 - 78.199 * s^3 +
    93.419 * s^2 + 19.825 * s + 1.64
  s <- c(-0.3, -0.11, 0, 0.11, 0.3)
  expect_equal(1 / f(s), h(s), tolerance = 1e-12)
  expect_lt(1 / f(-0.11), 1 / f(0))
  expect_true(all(f(seq(-1, 1, 0.01)) > 0))
})

test_that("Llobera-Sluckin (2007, eq. 16) use the signed slope", {
  f <- .reg()$ls$make(.par())
  q <- function(s) 2.635 + 17.37 * s + 42.37 * s^2 - 21.43 * s^3 + 14.93 * s^4
  s <- c(-0.3, -0.1, 0, 0.1, 0.3)
  expect_equal(1 / f(s), q(s), tolerance = 1e-12)
  expect_true(all(f(seq(-1, 1, 0.01)) > 0))
})

test_that("Ardigo et al. (2003, eq. 1) use the signed gradient and 45.72", {
  f <- .reg()$a$make(.par(V = 1.2))
  cw <- function(i, v) 1.866 * exp(4.911 * i) * v^2 - 3.773 * exp(3.416 * i) * v +
    45.72 * i^2 + 18.90 * i + 4.456
  i <- c(-0.3, -0.1, 0, 0.1, 0.3)
  expect_equal(1 / f(i), cw(i, 1.2), tolerance = 1e-12)
  f0 <- .reg()$a$make(.par(V = 0))
  expect_true(all(f(seq(-1, 1, 0.01)) > 0) && all(f0(seq(-1, 1, 0.01)) > 0))
})

test_that("Pandolf with downhill correction follows Yokota et al. (2004)", {
  W <- 70; L <- 0; V <- 1.2
  M <- function(G, N) 1.5 * W + 2 * (W + L) * (L / W)^2 + N * (W + L) * (1.5 * V^2 + 0.35 * V * G)
  CF <- function(G, N) N * ((G * (W + L) * V) / 3.5 - ((W + L) * (G + 6)^2) / W + (25 - V^2))
  for (N in c(1, 1.5)) {
    f <- .reg()$pcf$make(.par(N = N))
    # cost per metre (J/m) = metabolic rate (W) / speed (m/s)
    expect_equal(1 / f(-0.10), (M(-10, N) - CF(-10, N)) / V, tolerance = 1e-12)
    expect_equal(1 / f(0.10), M(10, N) / V, tolerance = 1e-12)
  }
  # with W = 70, L = 0, V = 1.2: 195 W at -10% against 256 W on the level
  f <- .reg()$pcf$make(.par())
  expect_equal(V / f(-0.10), 194.6, tolerance = 1e-3)
  expect_equal(V / f(0), 256.2, tolerance = 1e-9)
  # positive over the whole slope range, for several loads and speeds
  for (L2 in c(0, 30)) for (V2 in c(0, 0.8, 1.5)) {
    g <- .reg()$pcf$make(.par(L = L2, V = V2))
    expect_true(all(g(seq(-1, 1, 0.01)) > 0), label = paste("pcf L", L2, "V", V2))
  }
})

test_that("Van Leusen (2002) is symmetric about -10% and Pandolf about 0", {
  vl <- .reg()$vl$make(.par())
  expect_equal(vl(-0.05), vl(-0.15), tolerance = 1e-12)
  expect_false(isTRUE(all.equal(vl(0.05), vl(-0.05))))
  p <- .reg()$p$make(.par())
  expect_equal(p(-0.1), p(0.1), tolerance = 1e-12)   # |G|, as in Herzog (2020)
})

test_that("anisotropic functions distinguish ascent from descent", {
  reg <- .reg(); p <- .par()
  aniso <- c("t", "tofp", "mp", "icmonp", "icmoffp", "icfonp", "icfoffp", "ma",
             "alb", "r", "ks", "trp", "pcf", "m", "hrz", "vl", "ls", "a", "h")
  for (nm in aniso) {
    f <- reg[[nm]]$make(p)
    expect_false(isTRUE(all.equal(f(0.15), f(-0.15))), label = nm)
  }
  symm <- c("ug", "gkrs", "wcs", "ree", "b", "e", "p")
  for (nm in symm) {
    f <- reg[[nm]]$make(p)
    expect_equal(f(0.15), f(-0.15), tolerance = 1e-12, label = nm)
  }
})

test_that("energy functions accumulate energy in the declared unit", {
  dtm <- .flat_dtm()
  origin <- .pt(dtm, 21, 21)
  east <- terra::cellFromRowCol(dtm, 21, 31)          # 100 m on the level
  acc <- function(f, ...) {
    terra::values(mc_accum(mc_surface(dtm, funct = f, move = 4, ...), origin)$accum)[east]
  }
  # Pandolf: (M / V) J/m * 100 m, in kJ
  M0 <- 1.5 * 70 + 70 * 1.5 * 1.2^2
  expect_equal(acc("p"), M0 / 1.2 * 100 / 1000, tolerance = 1e-9)
  expect_equal(acc("pcf"), M0 / 1.2 * 100 / 1000, tolerance = 1e-9)
  # V = 0: speed from Tobler on level ground
  v0 <- .tobler_ms(0)
  M0t <- 1.5 * 70 + 70 * 1.5 * v0^2
  expect_equal(acc("p", V = 0), M0t / v0 * 100 / 1000, tolerance = 1e-9)
  # Hare: kcal per km over 0.1 km
  expect_equal(acc("h"), (48 + 30 / (6 * exp(-3.5 * 0.05))) * 0.1, tolerance = 1e-9)
  # Llobera-Sluckin: kJ/m * 100 m
  expect_equal(acc("ls"), 2.635 * 100, tolerance = 1e-9)
  tab <- mc_cost_functions()
  u <- setNames(tab$cost_unit, tab$code)
  expect_equal(unname(u[c("p", "pcf", "vl", "ls", "m", "hrz", "a", "h")]),
               c("kJ", "kJ", "kJ", "kJ", "J/kg", "J/kg", "J/kg", "kcal"))
})
