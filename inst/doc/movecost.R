## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 6,
  dpi = 96,
  fig.align = "center"
)
# Build the figures only if the spatial stack is present, so the vignette
# never fails on a minimal check machine.
have_deps <- requireNamespace("terra", quietly = TRUE) &&
  requireNamespace("sf", quietly = TRUE) &&
  requireNamespace("igraph", quietly = TRUE) &&
  requireNamespace("ggplot2", quietly = TRUE)
knitr::opts_chunk$set(eval = have_deps)

## ----start--------------------------------------------------------------------
library(movecost)

dtm    <- mc_volc()        # sample DTM (a volcano), a terra SpatRaster
origin <- mc_volc_loc()    # a start location, an sf point
destin <- mc_destin_loc()  # nine destination locations, sf points

## ----task1--------------------------------------------------------------------
surf <- mc_surface(dtm, funct = "t", move = 16)
acc  <- mc_accum(surf, origin = origin, time = "m", breaks = 2)
plot(acc)

## ----task2--------------------------------------------------------------------
surf_N <- mc_surface(dtm, funct = "t", move = 16, N = 1.19)
acc_N  <- mc_accum(surf_N, origin = origin, time = "m", breaks = 2)
plot(acc_N)

## ----task3--------------------------------------------------------------------
surf_e <- mc_surface(dtm, funct = "pcf", move = 16, W = 60, L = 5)
acc_e  <- mc_accum(surf_e, origin = origin)
plot(acc_e)

## ----task4--------------------------------------------------------------------
surf_v <- mc_surface(dtm, funct = "pcf", move = 16, W = 60, L = 5, V = 0)
acc_v  <- mc_accum(surf_v, origin = origin)
plot(acc_v)

## ----task5--------------------------------------------------------------------
surf_c <- mc_surface(dtm, funct = "t", move = 16, cogn.slp = TRUE)
acc_c  <- mc_accum(surf_c, origin = origin, time = "m", breaks = 2)
plot(acc_c)

## ----task6--------------------------------------------------------------------
lcp <- mc_paths(surf, origin = origin, destin = destin, time = "m")
plot(lcp)

## ----task7--------------------------------------------------------------------
lcp_b <- mc_paths(surf, origin = origin, destin = destin, time = "m",
                  return.base = TRUE)
plot(lcp_b)

## ----task8--------------------------------------------------------------------
surf_wcs <- mc_surface(dtm, funct = "wcs", move = 16, sl.crit = 10)
corr <- mc_corridor(surf_wcs, a = origin, b = destin[2, ])
plot(corr)

## ----task9--------------------------------------------------------------------
cmp <- mc_comp(dtm, origin = origin, destin = destin,
               functs = c("t", "ma", "ug", "gkrs"), move = 16)
plot(cmp)

## ----task9b, fig.height = 4---------------------------------------------------
plot(cmp, type = "chart")

## ----task10, eval = FALSE-----------------------------------------------------
# boundary <- mc_etna_boundary()                 # sample study-area polygon
# dtm_etna <- mc_dtm(boundary, z = 9)            # download elevation data
# cmp_etna <- mc_comp(dtm_etna,
#                     origin = mc_etna_start(),
#                     destin = mc_etna_end(),
#                     functs = c("t", "wcs"), move = 16)
# plot(cmp_etna)

## ----task11-------------------------------------------------------------------
malta_b <- mc_malta_dtm()
springs_b <- mc_springs()
surf_b <- mc_surface(malta_b, funct = "tofp", move = 8)
bnd <- mc_boundary(surf_b, origin = springs_b[c(5, 15, 30), ],
                   limit = 30, time = "m")
plot(bnd)
bnd

## ----task12-------------------------------------------------------------------
surf_a <- mc_surface(dtm, funct = "a", move = 16)
al <- mc_alloc(surf_a, origin = destin[c(3, 7, 9), ])
plot(al)

## ----task13-------------------------------------------------------------------
plot(al, isolines = TRUE)

## ----task14a------------------------------------------------------------------
surf_t <- mc_surface(dtm, funct = "t", move = 16)
nw_all <- mc_network(surf_t, nodes = destin, type = "allpairs")
plot(nw_all)

## ----task14b------------------------------------------------------------------
nw_nei <- mc_network(surf_t, nodes = destin, type = "neigh")
plot(nw_nei)

## ----task14c------------------------------------------------------------------
round(nw_nei$cost.matrix, 2)

## ----task15-------------------------------------------------------------------
nw_d <- mc_network(surf_t, nodes = destin, type = "allpairs", density = TRUE)
plot(nw_d, type = "density")

## ----task16-------------------------------------------------------------------
malta <- mc_malta_dtm()
springs <- mc_springs()
surf_malta <- mc_surface(malta, funct = "t", move = 8)
lcp_coast  <- mc_paths(surf_malta, origin = springs[5, ], destin = springs[15, ])
plot(lcp_coast)

## ----task17-------------------------------------------------------------------
# a path to reuse as a barrier
surf8 <- mc_surface(dtm, funct = "t", move = 8)
bar   <- mc_paths(surf8, origin = destin[1, ], destin = destin[4, ])$paths

# without the barrier
free <- mc_paths(surf8, origin = destin[3, ], destin = destin[6, ])
plot(free)

# with the barrier (built into a new surface)
surf_bar <- mc_surface(dtm, funct = "t", move = 8, barrier = bar)
blocked  <- mc_paths(surf_bar, origin = destin[3, ], destin = destin[6, ])
plot(blocked, plot.barrier = TRUE)

## ----task18-------------------------------------------------------------------
rk <- mc_rank(surf8, origin = destin[3, ], destin = destin[6, ], k = 3)
plot(rk)
rk$paths                                   # rank, cost, length

## ----task19-------------------------------------------------------------------
rk2 <- mc_rank(surf8, origin = destin[1, ], destin = destin[4, ], k = 3)
plot(rk2, type = "corridor")

## ----task20-------------------------------------------------------------------
bar2 <- mc_paths(surf8, origin = destin[9, ], destin = destin[2, ])$paths
surf_bar2 <- mc_surface(dtm, funct = "t", move = 8, barrier = bar2)
rk3 <- mc_rank(surf_bar2, origin = destin[1, ], destin = destin[4, ], k = 3)
plot(rk3)

## ----task21, fig.height = 4---------------------------------------------------
rk_h <- mc_rank(surf8, origin = destin[1, ], destin = destin[4, ],
                k = 5, time = "m")
plot(rk_h, type = "chart")

## ----io, eval = FALSE---------------------------------------------------------
# # the 26 cost functions, with families, units, and parameter usage
# mc_cost_functions()
# 
# # persist a result across sessions (use mc_save, NOT saveRDS: SpatRasters
# # hold external pointers)
# mc_save(acc, "acc.rds")
# acc <- mc_load("acc.rds")
# 
# # export to GIS formats (GeoTIFF + GeoPackage)
# mc_export(lcp, dir = "outputs")

