# movecost 3.1.1

Calculation of slope-dependent accumulated cost surfaces, least-cost paths,
least-cost corridors, least-cost networks, ranked alternative routes, cost
allocation, and cost boundaries related to human movement across the
landscape, with 26 selectable cost functions. See Alberti (2019)
<doi:10.1016/j.softx.2019.100331> for an earlier version of the package.

Version 3.0.0 is a ground-up architectural redesign of the 2.x series. The
analytical capabilities and the cost functions are retained; what changes is
how the package computes, scales, and visualises. See `NEWS.md` for the full
list of changes and the 2.x-to-3.0 function correspondence table.

Version 3.1.0 adds a long-requested feature: the terrain factor can now vary
across the landscape, supplied as a raster of terrain factors or of terrain
classes (with a lookup table), alongside the single-value terrain factor `N`.

Version 3.1.1 is a corrections release: twelve cost functions are brought into
line with their published equations (mainly the sign of the slope, which makes
descents cheaper than ascents where the sources say so), and the energy
functions now accumulate energy in stated units. Results change for those
functions; see `NEWS.md`.

## Design in one paragraph

Build the cost graph **once** with `mc_surface()`; reuse it with every
analysis function; visualise **on demand** with `plot()`. Computation and
visualisation are fully decoupled: every function returns an object that
stores its results, and the `plot()`/`autoplot()` methods (ggplot2) can be
invoked, customised, and re-invoked at any time without re-running anything.

## Quick start

```r
library(movecost)

# sample data shipped with the package
dtm    <- mc_volc()          # SpatRaster
start  <- mc_volc_loc()      # sf points
destin <- mc_destin_loc()    # sf points

# 1. build the cost surface ONCE (Tobler's on-path hiking function)
surf <- mc_surface(dtm, funct = "t", move = 16)

# 2. reuse it across analyses - nothing is recomputed
acc <- mc_accum(surf, origin = start, breaks = 0.05)   # isochrones
lcp <- mc_paths(surf, origin = start, destin = destin) # least-cost paths
nw  <- mc_network(surf, nodes = destin)                # LCP network + cost matrix
al  <- mc_alloc(surf, origin = destin)                 # cost allocation
rk  <- mc_rank(surf, start, destin[2, ], k = 4)        # ranked alternatives
bd  <- mc_boundary(surf, destin[1:3, ], limit = 5, time = "m")
co  <- mc_corridor(surf, a = destin[1, ], b = destin[4, ])

# 3. visualise on demand; plots are ggplot objects, customise at will
plot(acc)
p <- plot(lcp) + ggplot2::ggtitle("My custom title")
p
ggplot2::ggsave("lcp.pdf", p, width = 7, height = 7)

# a terrain factor varying across the landscape (here a fictitious
# two-class raster: lawn grass below 120 m, rock above), via a lookup table
terr <- terra::classify(dtm, rbind(c(-Inf, 120, 1), c(120, Inf, 2)))
tab  <- data.frame(class = 1:2, N = c(1.03, 2.50))  # see mc_terrain_factors()
surf_t <- mc_surface(dtm, funct = "t", terrain = terr, terrain.table = tab)
plot(surf_t, type = "terrain")

# compare cost functions
cmp <- mc_comp(dtm, start, destin[2, ], functs = c("t", "ks", "hrz"))
plot(cmp); plot(cmp, type = "chart")

# the 26 implemented cost functions, annotated
mc_cost_functions()

# persistence and GIS export
mc_save(surf, "surf.rds"); surf <- mc_load("surf.rds")
mc_export(lcp, dir = "outputs")
```

## Installation (from this source folder)

```r
install.packages(c("terra", "sf", "igraph", "ggplot2", "devtools"))
devtools::install()        # run from the package root
# optional extras:
install.packages("elevatr")   # online DTM download via mc_dtm()
```

To (re)generate the help files from the roxygen2 annotations and run the
full quality pipeline:

```r
devtools::document()   # builds man/*.Rd from the roxygen2 annotations
devtools::test()       # runs the analytic + consistency test suite
devtools::check()      # full R CMD check
```

## Why 3.0.0 is faster and more scalable

1. **Compute-once architecture**: the conductance graph is built a single
   time per cost function; in 2.x it was rebuilt inside every call (an
   n-location network implied O(n) redundant constructions).
2. **Batched multi-source queries**: allocation, boundaries, and cost
   matrices are served by single igraph Dijkstra calls over all origins.
3. **Modern stack**: terra (C++-backed raster I/O and algebra) replaces
   raster/sp; igraph replaces the gdistance accumulation machinery; chron is
   replaced by an internal formatter. Legacy sp/raster *inputs* are still
   accepted and converted automatically.
4. **NoData-aware graph**: cells with NoData never enter the graph, so
   irregular DTMs (coastlines) need no special handling - and the graph gets
   smaller, not larger, when the study area is irregular.

## Validation

The computational core was validated against analytic results on synthetic
terrains, against internal consistency laws, and against the 2.x
(gdistance-based) procedure: see `inst/validation/` and the testthat suite.
Two deliberate methodological corrections (corridor directionality;
signed-slope handling in the `ks` and `ma` cost functions) are documented in
`NEWS.md` and in `?mc_cost_functions`.
