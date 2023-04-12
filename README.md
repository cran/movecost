# movecost
vers 2.0


## History

`version 1.9`:
* the 'netw.type' parameter added to 'movenetw()' to allow the users to select which type of network will be calculated;
* amendments to the help documentation of the 'origin' and 'destination' parameters in the 'movebound()', 'movecomp()', 'movenetw()', and 'moverank()' functions;
* the dependence on the following packages has been dropped: rgdal, maptools, rgeos, spatstat.geom

`version 1.8`:
* 'moverank()' function added;
* facility to prevent LCPs to cross NoData areas (i.e., areas corresponding to the sea) added (controlled by the 'irregular.dtm' parameter);
* facility to plot the distribution of the length and cost of LCPs produced by different cost function added
to  'movecomp()' (controlled by the 'add.chart' parameter);
* the 'LCPs' and 'LCPs.back' data returned by the 'movecomp()' function now store the cost associated to each path (the variable is labelled 'cost');
* fix to an error to the calculation of the closest location in the 'movenetw()' function; the error was affecting the rendered plot;
* Tripcevich's hiking function added;
* 'malta_dtm_40' and 'springs' datasets added;
* in all the rendered plots, for visulisation purposes, the hillshade raster has been replaced by a slopeshade raster;
* updates, amendements, and improvements to the help documentation.


`version 1.7`:
* Hare's metabolic cost function and Eastman's abstract cost function added (the number of implemented cost functions is now 25);
* facility to take into account barriers (i.e., areas where the movement is inhibited) added to all the functions ('movealloc()' excluded);
* 'movecost()' now returns the conductance TransitionalLayer; this is instrumental to the internal optimisation of the 'movenetw()' function;
* the 'movenetw()' now also calculates LCPs network between pairs of neighboring locations;
* the 'movenetw()' function now returns a matrix reporting the cost incurred when moving between all the locations pair-wisely;
* internal optimisation of the 'movenetw()' function, which is now about 3 times faster than in ver.1.6;
* 'oneplot' parameter dropped from the 'movenetw()' function;
* fix to an error to the export facility in the 'movenetw()' function;
* fixes, improvements, and updates to the help documentation;
* link to the package's vignette added to the help documentation of all the implemented functions.

`version 1.6`:
* fix to a projection issue in the in-built 'volc' raster dataset;
* internal optimizations in the 'movenetw()' function;
* improvements layout- and content-wise to the functions' help documentation.

`version 1.5`:
* 'movenetw()': fix to an error in the way the exported GeoTiff files were given a name.
* minor adjustements and amendements to the help documentation of 'movenetw()'.
* 'volc' RasterLayer added to the body of sample datasets.

`version 1.4`:
* The 'movecomp()' and 'movenetw() functions have been added.
* Cost functions added: Minetti et al.'s, Marin Arroyo's, Ardigo et al.'s, Pandolf et al.'s with correction factor (the total number of implemented functions is now 23).
* The option has been added to 'movebound()' to calculate the perimeter and area enclosed by each cost-limit boundary.
* Help documentions updated accordingly and improved layout-wise.


`version 1.3`:
* The 'movealloc()' function has been added.
* Kondo-Seino's modified Tobler hiking function added.
* The Pandolf et al. and Van Leusen's cost functions can now internally work out the walking speed on the basis of the Tobler function (on-path hiking); to achieve this, the user has to simply set the V parameter to 0.
* Help documentation updated accordingly.
* In the corridor raster produced by 'movecorr()', when multiple origins are used, the point locations are now plotted on top of the rendered plot.


`version 1.2`:
* 'movebound()' function added.
* In the plots rendered by 'movecost()' and 'movecorr()', the accumulated cost surface, the DTM, and the least-cost corridor are overlaid by a hillshade raster to better convey the features of the terrain; the transparency of the hillshade can be adjusted using the 'transp' parameter.
* In the plot showing the least-cost path(s) rendered by 'movecost()', destination locations' label(s) shows time in sexagesimal numbers (hours, minutes, seconds).
* The destination location layer ('dest.loc.w.cost') returned by 'movecost()' is now included in the layers than can be exported using the 'export' parameter.
* In the destination location layer ('dest.loc.w.cost') returned by 'movecost()', when the cost is expressed in terms of time, a new variable ('cost_hms') is added (besides 'cost') expressing time in sexagesimal numbers. The variable 'cost' still stores the cost at the destination location(s), but it is expressed in decimal numbers.
* Minor adjustments to the subtitle of the plot rendered by 'movecorr()' to indicate the used terrain factor.
* Minor adjustments to the annotations of the LCP plot rendered by 'movecost()'.
* Minor edits and amendements to the help documentation of 'movecost()' function.


`version 1.1`:
* Fix to an error in the Bellavia's cost function, which has been introduced in ver 1.0.
* The use of a terrain factor (N) has been extended to almost all the implemented cost functions; help documentation updated accordingly.
* Minor adjustements to the returned plots to indicate the value of the used terrain factor.
* Improvements to the content and layout of the 'movecost()'s help documentation.


`version 1.0`:
* The option to calculate and plot the LCPs back to the origin has been added to 'movecost()'.
* The following cost functions have been implemented in 'movecost()' (and, consequently, in 'movecorr()'): Garmy, Kaddouri, Rozenblat, and Schneider's; Rees'; Bellavia's. The total number of the implemented cost functions is now eighteen.
* Help documentation updated accordingly, with bibliographical references for the new cost fucntions added.
* Minor amendements and improvements to some of the labels attached to the returned plots.


`version 0.9`:
* The 'movecorr()' function now can produce either a least-cost corridor between two locations, or corridors between multiple locations.
* The option has been added to 'movecorr()' to rescale the least-cost values in the output raster to be constrained between 0 and 1.
* The ficility to use the so-called 'cognitive slope' in place of the true slope has been added to 'movecost()' and 'movecorr()'.


`version 0.8`:
* A bug has been fixed that prevented the downloaded DTM to be returned in the list provided by the 'movecorr()' function.
* The option to customize the labels attached to the point locations on the output map has been added to 'movecorr()'.
* 'movecorr()' now optionally exports the DTM (if not provided by the user but acquired online), the two LCPs, and the accumulated cost surfaces around a and b.
* The accumulated cost rasters (around point a and around point b) are now returned by 'movecorr()'.
* The cost-surface is now returned (but not plotted) by 'movecost()'.
* Minor internal optimizations.
* Minor amendements and improvements to the help documentation.

`version 0.7`:
Function 'movecorr()' added. In the 'movecost()' function, the option has been added to have the function producing or not a graphical output (via the 'graph.out'parameter).The Marquez-Perez et al.'s modified Tobler's function has been renamed for consistency with the names given to the other functions: from 'mt' (which was standing for "modified Tobler"") to 'mp' (which stands for "Marquez-Perez"").

`version 0.6`:
The facility to download online elevation data (when an input DTM is not provided) has been added; it relies on the 'elevatr' package. Three new datasets (to be used to showcase the added functionality) have been added. Help documentation updated accordingly.

`version 0.5`:
Irmischer-Clarke's modified Tobler hiking functions for female added. This entailed modifying the name of the parameters as follows: `icmonp` for on-path male, `icmoffp` for off-path male, `icfonp` for on-path female, `icfoffp` for off-path female. Minor improvements to the layout of the help documentation.

`version 0.4`:
Parameter `move` added, which provides the option to set the number of directions in which cells are connected in the cost calculation.

`version 0.3`:
added cost functions: Llobera-Sluckin's metabolic energy expenditure; Alberti's Tobler hiking function modified for animal foraging excursions.

`version 0.2`: 
unlike the previous version, the cost calculation is not based on a slope raster (in degrees) derived from the input dtm; instead, the procedure described in van Etten, "R Package gdistance: Distances and Routes on Geographical Grids" in Journal of Statistical Software 76(13), 2017, pp. 14-15 is followed. The altitudinal difference between dtm cells is first calculated and then divided by cell centres distance, so obtaining slope values (expressed as rise over run) changing according to the direction of movement. This has the important implication that the implemented cost functions are now truly anisotropic. The color ramp of the output accumulated cost raster has been changed.

`version 0.1`: 
first release to CRAN.
