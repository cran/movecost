# movecost 1.3:

* The 'movealloc()' function has been added.
* Kondo-Seino's modified Tobler hiking function added.
* The Pandolf et al. and Van Leusen's cost functions can now internally work out the walking speed on the basis of the Tobler function (on-path hiking); to achieve this, the user has to simply set the V parameter to 0.
* Help documentation updated accordingly.
* In the corridor raster produced by 'movecorr()', when multiple origins are used, the point locations are now plotted on top of the rendered plot.
