# movecost 3.1.1

A corrections release. Twelve cost functions are brought into line with their
published equations; no feature is added or removed, and the interface is
unchanged. Results change for the functions listed below (correctly), so
analyses run with these functions in earlier versions should be re-run. Each
point is documented, with its source, in the help page of
`mc_cost_functions()`, which also states and justifies the few places where
movecost deliberately departs from a source.

## Corrected cost-function errors (results change)

The 3.0.0 review had checked the branch tests and coefficients of each cost
function against its source, but not the sign of the slope inside each
formula. Several functions applied the absolute value to a slope that the
source defines as signed (negative downhill), which made descents as costly as
ascents of the same gradient. Every function below has now been checked
against its primary source.

* Irmischer-Clarke functions (`icmonp`, `icfonp`, `icmoffp`, `icfoffp`): the
  percent slope is signed, as in equations 2-5 of Irmischer and Clarke (2018),
  so that speed peaks on a -5% (on-path) or -2% (off-path) descent.
* Rees (`r`): signed slope (Rees 2004, eq. 1). Small effect.
* Tripcevich (`trp`): signed slope in degrees, so that speed peaks on a
  -4.1 degree descent, as in the asymmetrical Cauchy function fitted by
  Tripcevich (2008); the reference is completed with the poster's details and
  URL.
* Minetti et al. (`m`): signed gradient (Minetti et al. 2002); at -10% the
  cost is now 1.13 instead of 4.90 J/(kg*m). Below -0.45, the lower end of the
  measured range, where the polynomial decreases and then turns negative, the
  gradient is held at -0.45, and `mc_surface()` reports how many edges are
  affected (documented departure from the source).
* Herzog (`hrz`): signed slope (Herzog 2016, 2020).
* Llobera-Sluckin (`ls`): signed slope (Llobera and Sluckin 2007, eq. 16).
* Ardigo et al. (`a`): signed gradient, and coefficient 45.72 instead of 45.71
  (Ardigo et al. 2003, eq. 1).
* Pandolf et al. with downhill correction (`pcf`): the grade is signed in both
  the Pandolf rate and the Yokota et al. correction, and the terrain factor
  multiplies the whole correction (Yokota et al. 2004, p. 4). Versions 3.0.0-3.1.0 used the absolute grade, so that a
  descent cost more than an ascent of the same gradient (e.g., 543 W instead
  of 195 W at -10%, with W = 70, L = 0, V = 1.2).
* Van Leusen (`vl`): the grade term is abs(G + 10), with G signed, placing the
  axis of cost symmetry at -10% as stated by Van Leusen (2002); 3.1.0 used
  abs(G) + 10, which equals it uphill but not downhill. The printed formula of
  the source reads abs(G + 6); the reasons for following the stated axis are
  given in the help page (documented departure from the source).

## Corrected units of the energy functions (results change)

* Pandolf (`p`), Pandolf with correction (`pcf`), Van Leusen (`vl`): these
  equations give a metabolic rate (W). Earlier versions accumulated
  rate x distance (W*m, labelled "Megawatts"); the cost per metre is now
  rate / speed (J/m), so that the accumulated value is energy, reported in kJ.
  With `V = 0` (speed derived from Tobler's function cell by cell) the
  least-cost paths of earlier versions did not minimise energy; they now do.
* Hare (`h`): the function gives kcal per km; the accumulated value is now
  divided by 1000 and reported in kcal (earlier versions reported the raw sum
  under the label "cal/km").
* The `cost_unit` column of `mc_cost_functions()`, the printed surfaces, and
  the plot legends now give the unit of the accumulated value (kJ, J/kg, kcal)
  rather than the per-metre unit (J/(kg*m), kJ/m).

## Exact speed conversion (results change slightly)

* Walking speeds are converted from km/h to m/s with the exact factor 1/3.6;
  earlier versions (including 2.x) used the rounded 0.278, which overstated
  speeds, and so understated walking times, by 0.08%.

## Documented, unchanged choices

* Pandolf et al. (`p`) keeps the absolute grade, the symmetric form reported by
  Herzog (2020, Table 18.2): the original equation covers level and uphill
  walking only (Yokota et al. 2004) and, with a negative grade, turns negative
  on steep descents. The help page says so and points to `pcf` for descents.
* Uriarte Gonzalez (`ug`), Garmy et al. (`gkrs`), Bellavia (`b`), Eastman
  (`e`), the wheeled-vehicle function (`wcs`), and the relative energetic
  expenditure (`ree`) are symmetric in their sources; the help page now states
  this where it was not stated, and that `ree` costs nothing on level ground.
* `pcf`: the help page notes the small step in cost between level and slightly
  downhill movement that comes from the source equation, and that Looney et al.
  (2019) print the last term of the correction differently from Yokota et al.
  (2004), whom movecost follows.

## Tests and validation

* New `tests/testthat/test-sources.R`: each corrected function is checked
  against its published equation, written out independently of the package
  code; anisotropic functions must distinguish ascent from descent and
  symmetric ones must not; energy functions must accumulate in the declared
  unit; every function must stay positive over slopes from -1 to 1.
* `inst/validation/` notes and the Python core replica updated accordingly.

# movecost 3.1.0

## New feature: a terrain factor that varies across the landscape

* `mc_surface()` gains two arguments, `terrain` and `terrain.table`, so that
  the terrain factor N no longer has to be a single value for the whole study
  area. `terrain` takes a raster (SpatRaster, RasterLayer, or file path) in
  either of two forms: a numeric raster whose cells hold terrain factors
  directly, or a categorical raster of terrain classes (integer codes, or a
  terra factor raster with category labels) translated into terrain factors
  by the lookup table `terrain.table` (columns `class` and `N`).
* The scalar `N` keeps its meaning and default, and now acts as the fallback
  value: cells with no terrain value (NA cells, cells outside the raster's
  extent, classes missing from the table) take `N`. A raster can therefore map
  only the terrain types of interest and leave the rest to a single general
  value. A message reports how many cells were filled; a warning lists any
  class missing from the table. `Inf` makes a cell impassable.
* Each edge takes the mean of the terrain factors of the two cells it joins.
  Since every implemented cost function yields an edge cost that is linear in
  N (affine for the Pandolf-type functions), the mean reproduces exactly the
  cost of crossing half of the edge on each terrain. N is applied inside the
  cost function, exactly as the scalar N is, rather than as a multiplier of
  the final cost: for the Pandolf (`p`), Pandolf with correction (`pcf`), and
  Van Leusen (`vl`) functions, where N multiplies only the terrain-dependent
  term of the metabolic equation, the two would differ.
* A terrain raster whose coordinate system or grid differs from the DTM's is
  aligned to the DTM with nearest-neighbour resampling (safe for class codes,
  and preserving the exact terrain factors of a numeric raster), with a
  warning.
* The returned surface stores the per-cell terrain factors actually used
  (component `terrain`); `print()` reports their range, the cost-function label
  carried by every downstream result reads e.g. "terrain factor
  N=raster (1.03 to 2.5)", `plot(surface, type = "terrain")` maps them, and
  `mc_export()` writes them to `terrain_N_<funct>.tif`. `mc_save()` /
  `mc_load()` preserve them.
* New `mc_terrain_factors()` returns, as a data frame, the reference list of
  terrain factors already documented in `mc_cost_functions()` (from Herzog
  2020, for time-based functions and for functions measuring other costs),
  as a starting point for building `terrain.table`.
* `mc_comp()` accepts `terrain` and `terrain.table` through `...`; the raster
  is aligned to the DTM once and shared by all the compared functions.
* All the analysis functions (`mc_accum()`, `mc_paths()`, `mc_corridor()`,
  `mc_boundary()`, `mc_alloc()`, `mc_network()`, `mc_rank()`) work unchanged on
  surfaces built with a terrain raster.

## Minor changes

* `mc_surface()` now validates `N`, which must be a single positive number
  (previously any value was accepted silently).
* The cost functions that do not take a terrain factor now warn when one is
  supplied (`N` other than 1, or `terrain`), which they ignore; previously `N`
  was ignored silently. The warning gives the reason: the off-path functions
  (`tofp`, `icmoffp`, `icfoffp`) already embed off-path walking, and the
  warning points to their on-path versions (`t`, `icmonp`, `icfonp`), which do
  take N; Alberti's function (`alb`) models the walking pace of a flock, to
  which a terrain factor does not apply.
* `print()` of a cost surface reports the terrain factor in use.

## Documentation

* `mc_surface()` help: new Details section on the terrain raster, and a worked
  example. `mc_cost_functions()`, `mc_comp()`, `mc_export()`,
  `plot.movecost_surface()` and the package help page updated accordingly.
* Vignette: new Task 3 ("Can the terrain factor vary across the study area?");
  the following tasks are renumbered. Task 10 now also shows `mc_comp()` with
  a terrain raster, the export section mentions the terrain raster written
  by `mc_export()`, and Task 2 explains what happens when a terrain factor is
  passed to a function that does not take one. The vignette title now reads
  "movecost 3.x: a task-based guide".
* `mc_cost_functions()` help: the exceptions to the terrain factor are now
  described accurately. Alberti's function (`alb`) was said to "natively
  embed" a terrain factor; its 0.25 coefficient is in fact the ratio between
  flock and human walking speed, and the function does not take N because it
  models a flock's pace. Each of the four functions concerned now states that
  it does not take N, and the off-path ones name their on-path version.
* `inst/CITATION`: the package version in the citation is now read from
  DESCRIPTION, so it can no longer fall out of date.
* Corrected the package-name capitalisation in the citation information
  (`inst/CITATION`) and in the vignette reference: the title of the 2019
  SoftwareX paper is "movecost" (the name carries no capital letters), not
  "MovEcost" or "Movecost".
* README: the release is now named with its full version string where the
  version itself is referenced.
* `mc_cost_functions()` help: completed and corrected the literature references
  for several cost functions. The Eastman function (`"e"`) now cites the
  original source, Eastman (1999: 61, the IDRISI32 guide), with Barge and
  Chataigner (2004) and Vaissie (2021) added as the secondary sources from
  which the equation was drawn; previously only Vaissie (2021) was given. The
  Bellavia function (`"b"`) now cites the primary Bellavia (2002) reference
  alongside Herzog (2020). The `gkrs` function now cites the primary Garmy,
  Kaddouri, Rozenblat and Schneider (2005) reference alongside Herzog (2020).
  The Tripcevich (`"trp"`) and Van Leusen (`"vl"`) references were completed
  with their previously missing publication venues. The Alberti function
  (`"alb"`) reference was corrected and completed to the full Alberti, Grima
  and Vella (2020) book-chapter citation (Temple Landscapes volume), replacing
  the incomplete and mis-dated "Alberti (2019)". No cost-function formula,
  result, or argument is affected.

# movecost 3.0.0

A ground-up redesign of the package's architecture, while retaining all the
analytical capabilities (and the 26 cost functions) of the previous releases.
The baseline for the comparisons below is the last version published on CRAN,
movecost 2.2.

## New architecture: compute once, reuse everywhere

* New `mc_surface()` builds the cost graph once; every analysis function
  (`mc_accum()`, `mc_paths()`, `mc_corridor()`, `mc_boundary()`, `mc_alloc()`,
  `mc_network()`, `mc_rank()`) reuses it. In the 2.x series each call to
  `movecost()`, `movecorr()`, `movealloc()`, `movebound()`, `movenetw()`,
  `moverank()` rebuilt the conductance matrix from scratch (an n-location
  network implied O(n) redundant constructions); that redundancy is gone.
* Multi-origin computations (allocation, boundaries, cost matrices) are
  served by single batched Dijkstra queries (igraph).

## Modern spatial stack

* The computational core now rests on terra + sf + igraph. The retired or
  maintenance-mode dependencies of the 2.x series (raster, sp, gdistance,
  chron, Matrix as a direct import, elevatr as a hard dependency) are gone;
  elevatr is now an optional (Suggested) dependency used only by `mc_dtm()`.
* Legacy sp/raster inputs are still accepted everywhere and coerced
  internally, so existing data files keep working.
* Sample datasets are now shipped as standard GIS files (GeoTIFF/GeoPackage,
  read via `mc_volc()`, `mc_malta_dtm()`, `mc_volc_loc()`, `mc_destin_loc()`,
  `mc_springs()`, `mc_etna_start()`, `mc_etna_end()`, `mc_etna_boundary()`),
  removing the last data-level dependency on raster/sp.

## Decoupled, ggplot2-based visualisation

* Computation and visualisation are fully separated: every result object is
  plotted on demand via `plot()` / `autoplot()` methods returning ggplot
  objects, which can be stored, customised with standard ggplot2 syntax, and
  re-drawn at will - without re-running any analysis (in 2.x, obtaining a new
  visualisation required re-running the computation).
* Styling is exposed through plot arguments: isoline colour and width
  (`iso.col`, default white, `iso.lwd`), least-cost-path width (`lcp.lwd`),
  and label sizes (`cex.breaks` for isoline labels, `cex.lcp.lab` for
  destination labels, `cex.node.lab` / `cex.origin.lab` for network and
  allocation labels).
* Network nodes and allocation origins are labelled with their indices by
  default (`node.lab`, `origin.lab`), matching `cost.matrix` and the
  allocation zones respectively, reproducing the 2.x numbering.

## Corrected cost-function errors (results change because 2.2 was wrong)

These are fixes to mistakes present in the last CRAN release (movecost 2.2),
verified against the source publications. Analyses using the affected
functions will give different (correct) results.

* Kondo-Seino (`ks`), Marin Arroyo (`ma`), and Pandolf-with-correction
  (`pcf`): the slope-sign tests are now applied to the signed slope, as the
  original publications intend. In 2.2 those tests were applied to the
  absolute slope, so the downhill branches (the gentler `ks` decay, the
  gentler `ma` divisor, the Yokota et al. downhill correction in `pcf`)
  could never be selected. Results on downhill gradients therefore differ -
  correctly - from 2.2.
* Pandolf-family functions (`p`, `vl`): the terrain factor N is now applied
  once, as the coefficient eta of Pandolf et al. (1977), multiplying the
  velocity-dependent term; 2.2 additionally multiplied the whole expression
  by N, counting terrain twice and inflating the basal-metabolism terms.
  With the default N = 1 results are unchanged; for N != 1 they differ
  (correctly) from 2.2.
* Hare metabolic cost function (`funct = "h"`): the 2.2 code contained a
  double-reciprocal that reduced algebraically to `1 / ((48 + 30 * v) * N)`
  rather than the intended `1 / ((48 + 30 / v) * N)`, where `v` is the
  cell-by-cell on-path Tobler walking speed. The corrected formula reproduces
  the body-mass-anchored approximation of Ralston (1958) cited by Hare (2004),
  agreeing with it to within ~5% across 2-6 km/h (the 2.2 code diverged by a
  factor of three or more). The plotted cost-unit label is corrected from
  "cal/km" to "kcal/km" to match.
* Metabolic functions with dynamic speed (`funct = "vl"`, `"p"`, `"pcf"`,
  `"a"` with `V = 0`): in 2.2 the value substituted for the velocity was the
  reciprocal of Tobler's on-path speed rather than the speed itself,
  mis-stating the metabolic rate (about half the canonical value on flat
  ground, a multiple on steep terrain). Version 3.0 substitutes Tobler's speed
  in m/s directly. Results obtained with `V = 0` therefore differ from 2.2;
  analyses using a non-zero `V` (including the default 1.2) are unaffected.

## Design changes and enhancements (method improved or generalised)

These are deliberate refinements of how a result is produced. Where 2.x
followed a specific published method, the 2.x behaviour is preserved as the
default and the refinement is offered as an option.

* Least-cost corridors (`mc_corridor()`): the classic definition of the 2.x
  series, the sum of two accumulated cost surfaces spread outward from the two
  locations after Mitchell (2012), is kept as the default
  (`method = "reach"`), is symmetric, and reproduces 2.x. A new
  `method = "through"` option adds a directional alternative,
  `cost(A->x) + cost(x->B)` (the second term on the reversed graph), whose
  minimum equals the A-to-B least-cost-path cost; it is intended for one-way
  journey analysis under anisotropic cost functions. The two coincide for
  isotropic functions. For the two-location case the least-cost paths A-to-B
  and B-to-A are computed, returned, and overlaid on the plot, as in 2.x.
* Cost-limit boundaries (`mc_boundary()`): the boundary is now obtained by
  polygonising the reachable area directly (cells at or below the limit),
  yielding closed polygons with exact area and perimeter; 2.x derived it from
  contour lines of the accumulated surface. The result is a more rigorous,
  measurable catchment polygon.
* Ranked sub-optimal paths (`mc_rank()`): same iterative-penalisation idea as
  2.x's `moverank()`, but operating directly on the graph edge weights, so any
  number of ranks can be requested (2.x was capped at 6) and the avoidance
  strength is exposed through a tunable `penalty` parameter.
* Irregular DTMs: NoData cells are excluded from the cost graph by
  construction, so accumulated costs and paths never cross them and the
  `irregular.dtm` workaround (and parameter) of 2.x is no longer needed.
* Least-cost networks (`mc_network()`): the path-density rasters are computed
  exactly from the graph vertices of each path (no rasterisation
  approximation), and the network nodes are numbered on the plot and in the
  returned `nodes` object, matching the row/column indices of `cost.matrix`,
  as in 2.x.
* Sample data: the bundled datasets carry explicit EPSG codes (e.g. Etna
  EPSG:23032, verified to fall on Mount Etna; volc EPSG:27200; Malta
  EPSG:23033), so their CRS is encoded identically across the GeoTIFF and
  GeoPackage files.

## Other improvements

* Exports now use GeoPackage instead of ESRI shapefile (no 10-character
  field-name truncation); `mc_export()` handles every result class.
* `mc_save()`/`mc_load()` provide reliable on-disk persistence of result
  objects (SpatRaster components cannot survive a bare `saveRDS()`).
* `mc_cost_functions()` returns the full annotated registry of the 26 cost
  functions (codes, families, units, parameter usage), which is also the
  single internal source of cost-function definitions (in 2.x the
  definitions were duplicated across functions).
* Walking-time formatting no longer depends on chron (internal HH:MM:SS
  formatter).
* Informative validation errors: geographic-CRS rasters, CRS mismatches,
  locations outside the DTM or on NoData cells, unknown cost-function codes.
* Comprehensive testthat suite with analytic and consistency checks, plus
  independent validation materials in `inst/validation/`.

## API correspondence with the 2.x series

| movecost 2.x                  | movecost 3.0                                |
|-------------------------------|---------------------------------------------|
| `movecost()` (accumulated)    | `mc_surface()` + `mc_accum()`               |
| `movecost()` (LCPs)           | `mc_surface()` + `mc_paths()`               |
| `movecorr()`                  | `mc_corridor()`                             |
| `movebound()`                 | `mc_boundary()`                             |
| `movealloc()`                 | `mc_alloc()`                                |
| `movenetw()`                  | `mc_network()`                              |
| `moverank()`                  | `mc_rank()`                                 |
| `movecomp()`                  | `mc_comp()`                                 |
| `studyplot=` (download DTM)   | `mc_dtm()` (or `studyplot=` in `mc_surface()`) |
| `export=TRUE`                 | `mc_export()`                               |

# movecost 2.2

* Added the `topo.dist` parameter to `movecost()`: when TRUE, accumulated cost
  uses topographic (surface) distance rather than planar distance, accounting
  for the extra distance over sloped terrain (surface = planar * sqrt(1 +
  slope^2); about +4.4% on a 30% slope, +11.8% on a 50% slope). Particularly
  recommended for time-based hiking functions; to be used with caution for
  metabolic functions, which were derived against planar distance.
* Added a documentation section explaining the methodology and its
  appropriateness per cost-function type.
* Propagated `topo.dist` to all wrapper functions (`movealloc()`,
  `movebound()`, `movecomp()`, `movecorr()`, `movenetw()`, `moverank()`); the
  default `FALSE` preserves backward compatibility.

# movecost 2.1

* Resolved lost-braces issues in `\itemize` sections across several Rd files.
* Corrected the `create_barrier_cs` documentation (missing `fields` argument
  in `\usage`).
* Minor help-documentation amendments.

# movecost 2.0

* Integrated `create_barrier_cs` (from an earlier version of the
  'leastcostpath' package) as an internal function; dependency on
  'leastcostpath' dropped; dependency on 'Matrix' added.
* Minor fixes and documentation amendments.

# movecost 1.9

* Added the `netw.type` parameter to `movenetw()`.
* Documentation amendments for the `origin`/`destination` parameters.
* Dropped the dependence on 'rgdal', 'maptools', 'rgeos', 'spatstat.geom'.

# movecost 1.8

* Added the `moverank()` function.
* Added the facility (via `irregular.dtm`) to prevent LCPs crossing NoData
  (sea) areas.
* Added length/cost distribution charts to `movecomp()` (via `add.chart`); the
  returned LCPs now store the `cost`.
* Fixed the closest-location calculation in `movenetw()`.
* Added the Tripcevich hiking function; added the `malta_dtm_40` and `springs`
  datasets.
* Replaced the hillshade with a slopeshade raster in all rendered plots.
* Documentation updates.

# movecost 1.7

* Added the Hare metabolic and Eastman abstract cost functions (25 functions
  total).
* Added barrier handling to all functions except `movealloc()`.
* `movecost()` returns the conductance transitional layer (used to optimise
  `movenetw()`).
* `movenetw()` also computes LCP networks between neighbouring locations and
  returns a pairwise cost matrix; internal optimisation (~3x faster than 1.6);
  `oneplot` dropped; export fix.
* Documentation fixes; vignette link added throughout.

# movecost 1.6

* Fixed a projection issue in the in-built `volc` dataset.
* Internal optimisations in `movenetw()`.
* Documentation improvements.

# movecost 1.5

* `movenetw()`: fixed exported-GeoTiff naming; documentation tweaks.
* Added the `volc` RasterLayer to the sample datasets.

# movecost 1.4

* Added `movecomp()` and `movenetw()`.
* Added cost functions: Minetti et al., Marin Arroyo, Ardigo et al., Pandolf
  et al. with correction factor (23 functions total).
* Added perimeter/area output to `movebound()`.

# movecost 1.3

* Added `movealloc()`.
* Added the Kondo-Seino modified Tobler function.
* Pandolf and Van Leusen functions can derive walking speed from the Tobler
  function by setting `V = 0`.
* `movecorr()` plots point locations on top of the multi-origin corridor.

# movecost 1.2

* Added `movebound()`.
* Overlaid a hillshade on the rendered surfaces (transparency via `transp`).
* LCP destination labels show time in sexagesimal form; added the `cost_hms`
  variable to the returned destination layer; that layer is now exportable.
* Minor plot and documentation amendments.

# movecost 1.1

* Fixed the Bellavia cost function (introduced in 1.0).
* Extended the terrain factor (N) to almost all cost functions.
* Plot and documentation improvements.

# movecost 1.0

* Added the option to compute and plot LCPs back to the origin.
* Added the Garmy et al., Rees, and Bellavia cost functions (18 total).
* Documentation updates with references.

# movecost 0.9

* `movecorr()` can produce a corridor between two locations or among multiple
  locations; added the 0-1 rescaling option.
* Added the 'cognitive slope' option to `movecost()` and `movecorr()`.

# movecost 0.8

* Fixed a bug preventing the downloaded DTM being returned by `movecorr()`.
* Added point-label customisation and optional export of DTM, LCPs, and
  accumulated surfaces to `movecorr()`; the accumulated rasters are now
  returned.
* `movecost()` now returns (without plotting) the cost surface.
* Minor optimisations and documentation amendments.

# movecost 0.7

* Added `movecorr()`; added the `graph.out` option to `movecost()`; renamed the
  Marquez-Perez function code from `mt` to `mp`.

# movecost 0.6

* Added online elevation-data download (via 'elevatr') when no DTM is supplied;
  added three sample datasets; documentation updated.

# movecost 0.5

* Added the Irmischer-Clarke female functions; renamed the parameters to
  `icmonp`, `icmoffp`, `icfonp`, `icfoffp`; documentation layout improvements.

# movecost 0.4

* Added the `move` parameter (number of directions in which cells connect).

# movecost 0.3

* Added the Llobera-Sluckin metabolic function and the Alberti modified Tobler
  function for animal foraging excursions.

# movecost 0.2

* Cost calculation now follows van Etten (2017), computing slope as the
  altitudinal difference between cells divided by the distance between cell
  centres (rise over run, direction-dependent), making the cost functions
  truly anisotropic. Output colour ramp changed.

# movecost 0.1

* First release to CRAN.
