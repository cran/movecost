# Validation materials

- `compare_v23.R` - cross-version regression check. For every cost function it
  builds the accumulated cost surface with the movecost 2.x core (replicated
  verbatim through gdistance) and with the movecost 3.0 core
  (`mc_surface()` + `mc_accum()`), on the bundled volc DTM, and prints the
  maximum relative difference with a MATCH / DIFFERS verdict. Unchanged
  functions agree to floating-point precision (~1e-15); the deliberately
  corrected functions (`ks`, `ma`, `pcf`, and `p`/`vl` when N != 1) show a
  known, explainable difference. Needs only gdistance + raster + sp installed.
  Run with:
  `source(system.file("validation", "compare_v23.R", package = "movecost"))`

- `validate_core.py` - independent Python/SciPy re-implementation of the
  movecost 3.0 computational core, validating it against analytic results,
  internal consistency laws, and a replica of the 2.x gdistance procedure on
  random rugged DTMs (developer material; excluded from the built package via
  .Rbuildignore; requires numpy + scipy).

- The package's testthat suite (`tests/testthat/`) mirrors the analytic and
  consistency checks in R; run it with `devtools::test()`.
