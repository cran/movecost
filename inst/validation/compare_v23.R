# ---------------------------------------------------------------------------
# Cross-version regression check: movecost 3.0 vs movecost 2.x
#
# For each cost function this script builds the accumulated cost surface with
#   (a) the movecost 2.x computational core, replicated verbatim through
#       gdistance (transition -> geoCorrection -> apply cost function ->
#       geoCorrection -> accCost), exactly as movecost() did in 2.x, and
#   (b) the movecost 3.0 core (mc_surface() + mc_accum()),
# then reports the maximum relative difference between the two surfaces on the
# bundled volc DTM (16-directions), and classifies each function as a MATCH or
# an expected DIFFERENCE.
#
# Only gdistance + raster + sp are required (not the old movecost itself): the
# pre-3.0 maths is reproduced here directly. The cost-function expressions
# below are transcribed verbatim from the pre-3.0 codebase, INCLUDING the
# errors that 3.0 corrects, so that the unchanged functions agree to machine
# precision and the corrected ones differ in a known, explainable way.
#
# EXPECTED OUTCOME:
#   * MATCH (max rel diff ~1e-15, i.e. floating-point noise) for every function
#     except the four below;
#   * DIFFERENCE for:
#       - ks   (Kondo-Seino): 2.x signed-slope test written on the absolute
#               slope, so the steep-downhill branch never fired;
#       - ma   (Marin Arroyo): same kind of dead downhill branch;
#       - pcf  (Pandolf with downhill correction): downhill correction never
#               applied in 2.x;
#       - p, vl with N != 1 (Pandolf, Van Leusen): 2.x applied the terrain
#               factor N twice; with N = 1 they MATCH, with N != 1 they differ.
# ---------------------------------------------------------------------------

stopifnot(requireNamespace("gdistance", quietly = TRUE),
          requireNamespace("raster", quietly = TRUE),
          requireNamespace("sp", quietly = TRUE))

suppressPackageStartupMessages(library(movecost))   # version 3.0

dtm3   <- mc_volc()
dtm2   <- raster::raster(dtm3)
origin <- mc_volc_loc()
origin_sp <- sf::as_Spatial(origin)

# ---- shared 2.x ingredients (computed once) --------------------------------
# the slope (rise/run, signed) transition layer, and the adjacency index, are
# the same for every cost function; gdistance warns that the altitude-
# difference transition has negative values (downhill), which is intended.
altDiff <- function(x) x[2] - x[1]
hd    <- suppressWarnings(gdistance::transition(dtm2, altDiff, directions = 16,
                                                symm = FALSE))
slope <- gdistance::geoCorrection(hd)
adj   <- raster::adjacent(dtm2, cells = seq_len(raster::ncell(dtm2)),
                          pairs = TRUE, directions = 16)

# 2.x default parameters
N <- 1; W <- 70; L <- 0; V <- 1.2; sl.crit <- 10

# ---- the function table ----------------------------------------------------
# group : "speed" (km/h, x0.278 then geoCorrection) or "recip" (1/cost)
# div   : unit divisor to match mc_accum(time = "h"): 3600 (time), 1e6 (MW), 1
# make  : the 2.x per-edge expression, verbatim, as a function of the slope
#         transition layer x (using x[adj], exactly as movecost 2.x did)
# pars  : list of mc_surface() parameters for the 3.0 side
tests <- list(
  list(code="t",    group="speed", div=3600, pars=list(),
       make=function(x) (6 * exp(-3.5 * abs(x[adj] + 0.05))) * (1/N)),
  list(code="tofp", group="speed", div=3600, pars=list(),
       make=function(x) (6 * exp(-3.5 * abs(x[adj] + 0.05))) * 0.6),
  list(code="mp",   group="speed", div=3600, pars=list(),
       make=function(x) (4.8 * exp(-5.3 * abs((x[adj] * 0.7) + 0.03))) * (1/N)),
  list(code="alb",  group="speed", div=3600, pars=list(),
       make=function(x) (6 * exp(-3.5 * abs(x[adj] + 0.05))) * 0.25),
  list(code="icmonp", group="speed", div=3600, pars=list(),
       make=function(x) ((0.11 + exp(-(abs(x[adj])*100 + 5)^2 / (2 * 30^2))) * 3.6) * (1/N)),
  list(code="icmoffp", group="speed", div=3600, pars=list(),
       make=function(x) (0.11 + 0.67 * exp(-(abs(x[adj])*100 + 2)^2 / (2 * 30^2))) * 3.6),
  list(code="icfonp", group="speed", div=3600, pars=list(),
       make=function(x) ((0.95 * (0.11 + exp(-(abs(x[adj]) * 100 + 5)^2/(2 * 30^2)))) * 3.6) * (1/N)),
  list(code="icfoffp", group="speed", div=3600, pars=list(),
       make=function(x) (0.95 * (0.11 + 0.67 * exp(-(abs(x[adj]) * 100 + 2)^2/(2 * 30^2)))) * 3.6),
  list(code="gkrs", group="speed", div=3600, pars=list(),
       make=function(x) (4 * exp(-0.008 * ((atan(abs(x[adj]))*180/pi)^2))) * (1/N)),
  list(code="r",    group="speed", div=3600, pars=list(),
       make=function(x) ((1 / (0.75 + 0.09 * abs(x[adj]) + 14.6 * (abs(x[adj]))^2)) * 3.6) * (1/N)),
  list(code="trp",  group="speed", div=3600, pars=list(),
       make=function(x) ((4.028*46^2)/(((atan(abs(x[adj]))*180/pi)+4.127)^2+46^2))*(1/N)),
  # --- amended: 2.x dead downhill branch (abs() in the test) ---
  list(code="ks",   group="speed", div=3600, pars=list(),
       make=function(x) ifelse(abs(x[adj]) >= -0.07, (5.1 * exp(-2.25 * abs(x[adj] + 0.07))) * (1/N), (5.1 * exp(-1.5 * abs(x[adj] + 0.07)))) * (1/N)),

  list(code="ug",   group="recip", div=3600, pars=list(),
       make=function(x) 1 / ((0.0277 * (abs(x[adj])*100) + 0.6115) * N)),
  # --- amended: 2.x dead downhill branch ---
  list(code="ma",   group="recip", div=3600, pars=list(),
       make=function(x) ifelse((abs(x[adj])*100) < 0, 1 / ((0.6 * ((abs(x[adj])*100)/23+1))*N), 1 / ((0.6 * ((abs(x[adj])*100)/11+1))*N))),

  list(code="wcs",  group="recip", div=1, pars=list(),
       make=function(x) 1 / ((1 + ((abs(x[adj])*100) / sl.crit)^2) * N)),
  list(code="ree",  group="recip", div=1, pars=list(),
       make=function(x) 1 / ((tan((atan(abs(x[adj]))*180/pi)*pi/180) / tan(1*pi/180)) * N)),
  list(code="b",    group="recip", div=1, pars=list(),
       make=function(x) 1 / (((atan(abs(x[adj]))*180/pi)+1) * N)),
  list(code="e",    group="recip", div=1, pars=list(),
       make=function(x) 1 / ((0.031 * (atan(abs(x[adj]))*180/pi)^2 - 0.025 * (atan(abs(x[adj]))*180/pi) + 1)*N)),
  list(code="m",    group="recip", div=1, pars=list(),
       make=function(x) 1 / (((280.5 * abs(x[adj])^5) - (58.7 * abs(x[adj])^4) - (76.8 * abs(x[adj])^3) + (51.9 * abs(x[adj])^2) + (19.6 * abs(x[adj])) + 2.5) * N)),
  list(code="hrz",  group="recip", div=1, pars=list(),
       make=function(x) 1 / (((1337.8 * abs(x[adj])^6) + (278.19 * abs(x[adj])^5) - (517.39 * abs(x[adj])^4) - (78.199 * abs(x[adj])^3) + (93.419 * abs(x[adj])^2) + (19.825 * abs(x[adj])) + 1.64) * N)),
  list(code="ls",   group="recip", div=1, pars=list(),
       make=function(x) 1 / ((2.635 + (17.37 * abs(x[adj])) + (42.37 * abs(x[adj])^2) - (21.43 * abs(x[adj])^3) + (14.93 * abs(x[adj])^4)) * N)),
  list(code="a",    group="recip", div=1, pars=list(),
       make=function(x) 1 / ((1.866 * exp(4.911*abs(x[adj])) * V^2 - 3.773 * exp(3.416*abs(x[adj])) * V + (45.71*abs(x[adj]^2) + 18.90*abs(x[adj])) + 4.456) * N)),
  list(code="h",    group="recip", div=1, pars=list(),
       make=function(x) 1 / ((48 + 30 / (6 * exp(-3.5 * abs(x[adj] + 0.05)))) * N)),
  list(code="p",    group="recip", div=1e6, pars=list(),
       make=function(x) 1 / ((1.5 * W + 2.0 * (W + L) * (L / W)^2 + N * (W + L) * (1.5 * (V^2) + 0.35 * V * (abs(x[adj])*100)))*N)),
  list(code="vl",   group="recip", div=1e6, pars=list(),
       make=function(x) 1 / ((1.5 * W + 2.0 * (W + L) * (L / W)^2 + N * (W + L) * (1.5 * (V^2) + 0.35 * V * ((abs(x[adj])*100) + 10)))*N)),
  # --- amended: 2.x dead downhill branch ---
  list(code="pcf",  group="recip", div=1e6, pars=list(),
       make=function(x) ifelse(abs(x[adj])*100 > 0 ,
              1 / (1.5 * W + 2.0 * (W+L) * (L/W)^2 + N * (W+L) * (1.5 * V^2 + 0.35 * V * (abs(x[adj])*100))),
              1 / ((1.5 * W + 2.0 * (W+L) * (L/W)^2 + N * (W+L) * (1.5 * V^2 + 0.35 * V * (abs(x[adj])*100))) - (N * ((abs(x[adj])*100) * (W+L) * V/3.5) - ((W+L) * ((abs(x[adj])*100)+6)^2/W) + (25-V^2)))))
)

# the functions for which a difference from 2.x is expected (default params)
expected_diff <- c("ks", "ma", "pcf")

# ---- run the comparison ----------------------------------------------------
compare_one <- function(tt) {
  # 2.x surface
  val <- slope
  val[adj] <- tt$make(slope)
  if (tt$group == "speed") val <- val * 0.278
  Cond <- suppressWarnings(gdistance::geoCorrection(val))
  a2 <- suppressWarnings(gdistance::accCost(Cond, sp::coordinates(origin_sp))) / tt$div
  # 3.0 surface
  surf <- do.call(mc_surface, c(list(dtm = dtm3, funct = tt$code, move = 16),
                                tt$pars))
  a3 <- terra::values(mc_accum(surf, origin, time = "h")$accum, mat = FALSE)
  v2 <- raster::values(a2)
  ok <- is.finite(v2) & is.finite(a3)
  rel <- abs(v2[ok] - a3[ok]) / pmax(abs(v2[ok]), 1e-12)
  max(rel)
}

cat(sprintf("%-9s %-14s %s\n", "function", "max rel diff", "verdict"))
cat(strrep("-", 45), "\n")
for (tt in tests) {
  md <- tryCatch(compare_one(tt), error = function(e) NA_real_)
  expect <- tt$code %in% expected_diff
  verdict <- if (is.na(md)) "ERROR (see above)"
             else if (md < 1e-6) "MATCH"
             else if (expect) "DIFFERS (expected)"
             else "DIFFERS (UNEXPECTED!)"
  cat(sprintf("%-9s %-14s %s\n", tt$code, format(md, digits = 3), verdict))
}

cat("\nNote: p and vl MATCH here because the default N = 1; the 2.x double-N\n")
cat("difference appears only for N != 1. To see it, re-run a single function\n")
cat("with N = 2 on both sides, e.g. compare mc_surface(dtm, 'p', N = 2) with\n")
cat("the 2.x expression evaluated at N = 2.\n")
