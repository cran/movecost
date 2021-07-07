# movecost 1.2:

* 'movebound()' function added.
* In the plots rendered by 'movecost()' and 'movecorr()', the accumulated cost surface, the DTM, and the least-cost corridor are overlaid by a hillshade raster to better convey the features of the terrain; the transparency of the hillshade can be adjusted using the 'transp' parameter.
* In the plot showing the least-cost path(s) rendered by 'movecost()', destination locations' label(s) shows time in sexagesimal numbers (hours, minutes, seconds).
* The destination location layer ('dest.loc.w.cost') returned by 'movecost()' is now included in the layers than can be exported using the 'export' parameter.
* In the destination location layer ('dest.loc.w.cost') returned by 'movecost()', when the cost is expressed in terms of time, a new variable ('cost_hms') is added (besides 'cost') expressing time in sexagesimal numbers. The variable 'cost' still stores the cost at the destination location(s), but it is expressed in decimal numbers.
* Minor adjustments to the subtitle of the plot rendered by 'movecorr()' to indicate the used terrain factor.
* Minor adjustments to the annotations of the LCP plot rendered by 'movecost()'.
* Minor edits and amendements to the help documentation of 'movecost()' function.
