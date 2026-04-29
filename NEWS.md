# movecost 2.2:

* Added `topo.dist` parameter to `movecost()`: when set to `TRUE`, the function calculates accumulated cost using topographic (surface) distance rather than planar (horizontal) distance. This accounts for the additional distance incurred when traversing sloped terrain, providing more realistic walking-time estimates in mountainous areas. The correction is calculated as surface distance = planar distance times the square root of (1 + slope squared). For example, on a 30% slope, surface distance is approximately 4.4% greater than planar distance; on a 50% slope, approximately 11.8% greater.

* The `topo.dist` parameter works with all available cost functions, though it is particularly recommended for time-based hiking functions. Users should exercise caution when applying it to metabolic energy expenditure functions, as these were typically derived empirically using planar distance measurements.

* Added detailed documentation section explaining the topographic distance methodology, its mathematical basis, and guidance on appropriateness for different cost function types.

* Propagated the new `topo.dist` parameter to all wrapper functions that rely on `movecost()` as internal workhorse: `movealloc()`, `movebound()`, `movecomp()`, `movecorr()`, `movenetw()`, and `moverank()`. Users can now control whether topographic (surface) or planar distance is used directly from each wrapper. The default value (`topo.dist = FALSE`) preserves backward compatibility, so existing code continues to behave exactly as in version 2.1.

