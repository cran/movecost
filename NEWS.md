# movecost 1.8:

* 'moverank()' function added;
* facility to prevent LCPs to cross NoData areas (i.e., areas corresponding to the sea) added (controlled by the 'irregular.dtm' parameter);
* facility to plot the distribution of the length and cost of LCPs produced by different cost function added
to  'movecomp()' (controlled by the 'add.chart' parameter);
* the 'LCPs' and 'LCPs.back' data returned by the 'movecomp()' function now store the cost associated to each path (the variable is labelled 'cost');
* fix to an error to the calculation of the closest location in the 'movenetw()' function; the error was affecting the rendered plot;
* 'malta_dtm_40' and 'springs' datasets added;
* in all the rendered plots, for visulisation purposes, the hillshade raster has been replaced by a slopeshade raster;
* updates, amendements, and improvements to the help documentation.
