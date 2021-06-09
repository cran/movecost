# movecost
vers 1.0


## History
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
