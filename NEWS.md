# movecost 1.7:

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

