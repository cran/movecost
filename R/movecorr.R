#' R function for calculating least-cost corridor between point locations
#'
#' The function provides the facility to calculate the least-cost corridor between point locations.
#' It just requires an input DTM and at least two point locations ('SpatialPointsDataFrame' class) representing the locations between which the corridor is calculated.
#' Under the hood, 'movecorr()' relies on the \code{\link{movecost}} function and, needless to say, implements the same
#' cost functions. See the help documentation of 'movecost()' for further details.\cr
#' Visit this \href{https://drive.google.com/file/d/1gLDrkZFh1b_glzCEqKdkPrer72JJ9Ffa/view?usp=sharing}{LINK} to access the package's vignette.\cr
#'
#' If only two locations are provided (one via parameter 'a', one via parameter 'b'),
#' the function renders a raster representing the least cost corridor (which can be optionally exported as GeoTiff) with least-cost paths superimposed.
#' If more than 2 locations are fed into the function via the 'a' parameter, the function calculates the least-cost corridor between pairs of locations.
#' All the pair-wise corridor rasters are returned (but not individually plotted) in a list.
#' All those rasters will be summed, and the resulting raster will be plotted (and can be, optionally, exported as GeoTiff).\cr
#'
#' The function returns a list containing a number of components (see 'Value' below). For more details about exporting the function's outputs,  see 'Arguments' below. \cr
#'
#' If the user wants to calculate the least-cost corridor between two locations only, (s)he may want to use parameter 'a' and 'b' to indicate
#' the two locations of interest respectively. For example, using the datasets provided by this package: \cr
#'
#' result <- movecorr(a=Etna_start_location, b=Etna_end_location[1,], studyplot=Etna_boundary, funct="tofp") \cr
#'
#' The above will produce the least-cost corridor between two locations close to Mt Etna (Sicily, Italy), using the
#' Tobler's cost function (for off-path hiking). Side note: the elevation data will be acquired online. \cr
#'
#' If the interest lies in using more than 2 locations, the user may want to feed the dataset storing all the locations
#' into parameter 'a' (disregarding 'b'). As explained above, in this case the function calculates the least-cost corridor between pairs of locations.
#' All the pair-wise corridor rasters are returned in a list. Those rasters will be summed, and the resulting raster will be plotted (and can be, optionally, exported as GeoTiff).
#' For example, to calculate the least-cost corridors between every individual unique pair of the 9 locations stored in the 'destin.loc' dataset:\cr
#'
#' volc <- raster::raster(system.file("external/maungawhau.grd", package="gdistance")) \cr
#'
#' result <- movecorr(dtm=volc, a=destin.loc, funct="ree", rescale=TRUE) \cr
#'
#' Note that only parameter 'a' has been used. The function returns and plots the sum of the 36 individual corridors; the latter are not plotted,
#' but are stored in a list. If the user wants to plot the least-cost corridor, say, n 4, and then add the two locations
#' between which the corridor has been calculated, (s)he can first plot the corridor raster n 4: \cr
#'
#' raster::plot(result$corridors[[4]]) \cr
#'
#' Then, identifying which locations are related to corridor n 4 can be easily accomplished by looking up the values stored in
#' the 4th column of the returned matrix: \cr
#'
#' result$locations.matrix \cr
#'
#' The locations are the n 1 and n 5, so the user can add them to the plot previosly produced using: \cr
#'
#' raster::plot(destin.loc[1,], pch=20, add=T)\cr
#' raster::plot(destin.loc[5,], pch=20, add=T)\cr
#'
#' Note that the resulting plot can be produced (with a nicer outlook) directly by 'movecorr()' by feeding those two locations in the
#' parameter 'a' and 'b' respectively: \cr
#'
#' result <- movecorr(dtm=volc, a=destin.loc[1,], b=destin.loc[5,], funct="ree") \cr
#'
#' Overall, what 'movecorr()' does is to calculate (via the \code{\link{movecost}} function) the accumulated cost surface around each location.
#' Those are eventually summed to produce the least-cost corridor between locations. On the produced corridor raster, the cost of a cell is the total cost to reach it
#' from all the analysed locations. About least-cost corridors between pairs of locations, see for instance: \cr
#' Mitchell A. (2012), "The ESRI Guide to GIS Analysis. Vol 3. Modelling Suitability, Movement, and Interaction", New York: Esri Press (257-259). \cr
#'
#'
#' @param dtm Digital Terrain Model (RasterLayer class); if not provided, elevation data will be acquired online for the area enclosed by the 'studyplot' parameter (see \code{\link{movecost}}).
#' @param a first location from which the least-cost corridor is calculated (SpatialPointsDataFrame class); if it contains more than two locations, see the 'Description' section above.
#' @param b second location from which the least-cost corridor is calculated (SpatialPointsDataFrame class); if parameter 'a' stores more than two locations, this parameter is disregarded; see the 'Description' section above.
#' @param lab.a string to be used to label point a on the outplut plot (A is the default)
#' @param lab.b string to be used to label point a on the outplut plot (B is the default).
#' @param cex.labs scaling factor for the size of the points' labels (0.8 by default)
#' @param studyplot polygon (SpatialPolygonDataFrame class) representing the study area for which online elevation data are aquired (see \code{\link{movecost}}); NULL is default.
#' @param barrier area where the movement is inhibited (SpatialLineDataFrame or SpatialPolygonDataFrame class) (see \code{\link{movecost}}.
#' @param funct cost function to be used (for details on each of the following, see \code{\link{movecost}}):\cr
#'
#' \strong{-functions expressing cost as walking time-}\cr
#' \strong{t} (default) uses the on-path Tobler's hiking function;\cr
#' \strong{tofp} uses the off-path Tobler's hiking function;\cr
#' \strong{mp} uses the Marquez-Perez et al.'s modified Tobler's function;\cr
#' \strong{icmonp} uses the Irmischer-Clarke's hiking function (male, on-path);\cr
#' \strong{icmoffp} uses the Irmischer-Clarke's hiking function (male, off-path);\cr
#' \strong{icfonp} uses the Irmischer-Clarke's hiking function (female, on-path);\cr
#' \strong{icfoffp} uses the Irmischer-Clarke's hiking function (female, off-path);\cr
#' \strong{ug} uses the Uriarte Gonzalez's walking-time cost function;\cr
#' \strong{ma} uses the Marin Arroyo's walking-time cost function;\cr
#' \strong{alb} uses the Alberti's Tobler hiking function modified for pastoral foraging excursions;\cr
#' \strong{gkrs} uses the Garmy, Kaddouri, Rozenblat, and Schneider's hiking function;\cr
#' \strong{r} uses the Rees' hiking function;\cr
#' \strong{ks} uses the Kondo-Seino's hiking function;\cr
#'
#' \strong{-functions for wheeled-vehicles-}\cr
#' \strong{wcs} uses the wheeled-vehicle critical slope cost function;\cr
#'
#' \strong{-functions expressing abstract cost-}\cr
#' \strong{ree} uses the relative energetic expenditure cost function;\cr
#' \strong{b} uses the Bellavia's cost function;\cr
#' \strong{e} uses the Eastman's cost function;\cr
#'
#' \strong{-functions expressing cost as metabolic energy expenditure-}\cr
#' \strong{p} uses the Pandolf et al.'s metabolic energy expenditure cost function;\cr
#' \strong{pcf} uses the Pandolf et al.'s cost function with correction factor for downhill movements;\cr
#' \strong{m} uses the Minetti et al.'s metabolic energy expenditure cost function;\cr
#' \strong{hrz} uses the Herzog's metabolic energy expenditure cost function;\cr
#' \strong{vl} uses the Van Leusen's metabolic energy expenditure cost function;\cr
#' \strong{ls} uses the Llobera-Sluckin's metabolic energy expenditure cost function;\cr
#' \strong{a} uses the Ardigo et al.'s metabolic energy expenditure cost function;\cr
#' \strong{h} uses the Hare's metabolic energy expenditure cost function (for all the mentioned cost functions, see \code{\link{movecost}});\cr
#' @param time time-unit expressed by the accumulated raster if Tobler's and other time-related cost functions are used;
#' 'h' for hour, 'm' for minutes.
#' @param move number of directions in which cells are connected: 4 (rook's case), 8 (queen's case), 16 (knight and one-cell queen moves; default).
#' @param field value assigned to the cells coincidinng with the barrier (0 by default) (see \code{\link{movecost}}.
#' @param cogn.slp  TRUE or FALSE (default) if the user wants or does not want the 'cognitive slope' to be used in place of the real slope (see \code{\link{movecost}}).
#' @param sl.crit critical slope (in percent), typically in the range 8-16 (10 by default) (used by the wheeled-vehicle cost function; see \code{\link{movecost}}).
#' @param W walker's body weight (in Kg; 70 by default; used by the Pandolf's and Van Leusen's cost function; see \code{\link{movecost}}).
#' @param L carried load weight (in Kg; 0 by default; used by the Pandolf's and Van Leusen's cost function; see \code{\link{movecost}}).
#' @param N coefficient representing ease of movement (1 by default) (see \code{\link{movecost}}).
#' @param V speed in m/s (1.2 by default) (used by the Pandolf et al.'s, Pandolf et al.s with correction factor, Van Leusen's, and Ardigo et al.'s cost function; if set to 0, it is internally worked out on the basis of Tobler on-path hiking function (see \code{\link{movecost}}).
#' @param z zoom level for the elevation data downloaded from online sources (from 0 to 15; 9 by default) (see \code{\link{movecost}} and \code{\link[elevatr]{get_elev_raster}}).
#' @param rescale TRUE or FALSE (default) if the user wants or does not want the output least-coast corridor raster to be rescaled between 0 and 1.
#' @param transp set the transparency of the hillshade raster that is plotted over the least-cost corridor raster (0.5 by default).
#' @param export TRUE or FALSE (default) if the user wants or does not want the output to be exported; if TRUE, the least-cost corridor, the dtm (if not provided by the user but acquired online),
#' and the accumulated cost surface around a and b are exported as a GeoTiff file, while the two LCPs (from a to b, and from b to a) as individual shapefiles. If multiple locations are analysed, only the
#' least-cost corridor (and the DTM if originally not provided) will be exported. All the exported files (excluding the DTM) will bear a suffix corresponding to the cost function selected by the user.
#'
#' @return The function returns a list storing the following components \itemize{
##'  \item{dtm: }{Digital Terrain Model ('RasterLayer' class)}
##'  \item{lc.corridor: }{raster of the least-cost corridor ('RasterLayer' class); if more than two locations are analysed, this raster is the sum of all the corridors between all the pairs of locations}
##'  \item{lcp_a_to_b: }{least-cost past from a to b ('SpatialLinesDataFrame' class); returned only when the corridor is calculated between two locations}
##'  \item{lcp_b_to_a: }{least-cost past from b to a ('SpatialLinesDataFrame' class); returned only when the corridor is calculated between two locations}
##'  \item{accum_cost_surf_a: }{accumulated cost-surface around a ('RasterLayer' class); returned only when the corridor is calculated between two locations}
##'  \item{accum_cost_surf_b: }{accumulated cost-surface around b ('RasterLayer' class); returned only when the corridor is calculated between two locations}
##'  \item{corridors: }{list of rasters ('RasterLayer' class) representing the least-cost corridor between all the unique pairs of locations; returned only when more than two locations are analysed}
##'  \item{locations.matrix: }{matrix whose columns indicate the identifiers for all the unique pairs of locations for which each corridor is calculated; returned only when more than two locations are analysed}
##' }
##'
#' @keywords movecorr
#' @export
#' @importFrom raster ncell mask crop stack cellStats raster hillShade terrain
#' @importFrom elevatr get_elev_raster
#' @importFrom graphics layout par
#' @importFrom utils combn setTxtProgressBar txtProgressBar
#' @importFrom grDevices terrain.colors topo.colors grey
#'
#' @examples
#' # load a sample Digital Terrain Model
#' data(volc)
#'
#'
#' # load the sample destination locations on the above DTM
#' data(destin.loc)
#'
#'
#' # calculate the least-cost corridor between two locations, using the
#' # relative energetic expenditure cost function, and store the results
#' # in the 'result' object
#'
#' result <- movecorr(dtm=volc, a=destin.loc[1,], b=destin.loc[3,], funct="ree", move=8)
#'
#'
#' #same as above, but using the 'cognitive slope'
#'
#' # result <- movecorr(dtm=volc, a=destin.loc[1,], b=destin.loc[3,],
#' # funct="ree", move=8, cogn.slp=TRUE)
#'
#'
#' @seealso \code{\link{movecost}}
#'
#'
movecorr <- function (dtm=NULL, a, b, lab.a="A", lab.b="B", cex.labs=0.8, studyplot=NULL, barrier=NULL, funct="t", time="h", move=16, field=0, cogn.slp=FALSE, sl.crit=10, W=70, L=0, N=1, V=1.2, z=9, rescale=FALSE, transp=0.5, export=FALSE){

  #if no dtm is provided
  if (is.null(dtm)==TRUE) {
    #get the elvation data using the elevatr's get_elev_raster() function, using the studyplot dataset (SpatialPolygonDataFrame)
    #to select the area whose elevation data are to be downloaded;
    #z sets the resolution of the elevation datataset
    elev.data <- elevatr::get_elev_raster(studyplot, z = z, verbose=FALSE, override_size_check = TRUE)
    #crop the elevation dataset to the exact boundary of the studyplot dataset
    dtm <- raster::crop(elev.data, studyplot)
  }

  #if the input dataset 'a' contains just 1 feature...
  if (length(a) < 2) {

  #calculate the accum cost surface and LCP around and from point a
  res.a <- movecost::movecost(dtm=dtm, origin=a, destin=b, studyplot=studyplot, barrier=barrier, funct=funct, time=time, move=move, field=field, cogn.slp=cogn.slp, sl.crit=sl.crit, W=W, L=L, N=N, V=V, z=z, return.base=FALSE, graph.out=FALSE)

  #calculate the accum cost surface and LCP around and from point b
  res.b <- movecost::movecost(dtm=dtm, origin=b, destin=a, studyplot=studyplot, barrier=barrier, funct=funct, time=time, move=move, field=field, cogn.slp=cogn.slp, sl.crit=sl.crit, W=W, L=L, N=N, V=V, z=z, return.base=FALSE, graph.out=FALSE)

  #combine the 2 accumulated cost surfaces obtained at the preceding steps
  res.corridor <- res.a$accumulated.cost.raster + res.b$accumulated.cost.raster

  }

  #if the input dataset 'a' contains more than 2 features...
  if (length(a) > 2) {

    #create an empty list to store the calculated accumulated cost surfaces
    corridors <- list()

    #create a matrix of all the pair-wise locations index, to be used as indexes
    #inside the following 'for' loop
    pairw.list <- combn(seq(1:length(a)),2)
    pairw.list <-as.matrix(pairw.list)

    message(paste0("Wait...processing ", ncol(pairw.list), " unique pairs of locations..."))

    #set the progress bar
    pb <- txtProgressBar(min = 0, max = length(a), style = 3)

    #loop through all the locations to calculate the accumulated cost surfaces
    #and store them in the list
    for (i in 1:(ncol(pairw.list))) {
      acc.cost.srf.a <- movecost::movecost(dtm=dtm, origin=a[pairw.list[,i][1],], studyplot=studyplot, barrier=barrier, funct=funct, time=time, move=move, field=field, cogn.slp=cogn.slp, sl.crit=sl.crit, W=W, L=L, N=N, V=V, return.base=FALSE, graph.out=FALSE)$accumulated.cost.raster
      acc.cost.srf.b <- movecost::movecost(dtm=dtm, origin=a[pairw.list[,i][2],], studyplot=studyplot, barrier=barrier, funct=funct, time=time, move=move, field=field, cogn.slp=cogn.slp, sl.crit=sl.crit, W=W, L=L, N=N, V=V, return.base=FALSE, graph.out=FALSE)$accumulated.cost.raster
      corridors[[i]] <- acc.cost.srf.a + acc.cost.srf.b
      setTxtProgressBar(pb, i)
    }

    #sum all the rasters stored in the list
    res.corridor <- sum(raster::stack(corridors))
  }

  #if rescale is TRUE...
  if (rescale=="TRUE") {
    #rescale the resulting raster to have values between 0 and 1
    res.corridor <- ((res.corridor - raster::cellStats(res.corridor, "min"))/(raster::cellStats(res.corridor, "max") - raster::cellStats(res.corridor, "min")))
  }

  #define the appropriate text to be used for subsequent plotting
  if (funct=="t") {

    #set the labels to be used within the returned plot
    main.title <- paste0("Least-cost corridor (cost in ", time, ")")
    sub.title <- paste0("Walking-time based on the Tobler's on-path hiking function \n terrain factor N=", N)
    legend.cost <- paste0("walking-time (", time,")")
  }

  if (funct=="tofp") {

    #set the labels to be used within the returned plot
    main.title <- paste0("Least-cost corridor (cost in ", time, ")")
    sub.title <- "Walking-time based on the Tobler's off-path hiking function"
    legend.cost <- paste0("walking-time (", time,")")
  }

  if(funct=="mp") {

    #set the labels to be used within the returned plot
    main.title <- paste0("Least-cost corridor (cost in ", time, ")")
    sub.title <- paste0("Walking-time based on the Marquez-Perez et al.'s modified Tobler hiking function \n terrain factor N=", N)
    legend.cost <- paste0("walking-time (", time,")")
  }

  if(funct=="icmonp") {

    #set the labels to be used within the returned plot
    main.title <- paste0("Least-cost corridor (cost in ", time, ")")
    sub.title <- paste0("Walking-time based on the (male, on-path) Irmischer-Clarke's hiking function \n terrain factor N=", N)
    legend.cost <- paste0("walking-time (", time,")")
  }

  if(funct=="icmoffp") {

    #set the labels to be used within the returned plot
    main.title <- paste0("Least-cost corridor (cost in ", time, ")")
    sub.title <- "Walking-time based on the (male, off-path) Irmischer-Clarke's hiking function"
    legend.cost <- paste0("walking-time (", time,")")
  }

  if(funct=="icfonp") {

    #set the labels to be used within the returned plot
    main.title <- paste0("Least-cost corridor (cost in ", time, ")")
    sub.title <- paste0("Walking-time based on the (female, on-path) Irmischer-Clarke's hiking function \n terrain factor N=", N)
    legend.cost <- paste0("walking-time (", time,")")
  }

  if(funct=="icfoffp") {
    main.title <- paste0("Least-cost corridor (cost in ", time, ")")
    sub.title <- "Walking-time based on the (female, off-path) Irmischer-Clarke's hiking function"
    legend.cost <- paste0("walking-time (", time,")")
  }

  if(funct=="gkrs") {
    main.title <- paste0("Least-cost corridor (cost in ", time, ")")
    sub.title <- paste0("Walking-time based on the Garmy et al.'s hiking function \n terrain factor N=", N)
    legend.cost <- paste0("walking-time (", time,")")
  }

  if(funct=="r") {
    main.title <- paste0("Least-cost corridor (cost in ", time, ")")
    sub.title <- paste0("Walking-time based on the Rees' hiking function \n terrain factor N=", N)
    legend.cost <- paste0("walking-time (", time,")")
  }

  if(funct=="ug") {

    #set the labels to be used within the returned plot
    main.title <- paste0("Least-cost corridor (cost in ", time, ")")
    sub.title <- paste0("Walking-time based on the Uriarte Gonzalez's hiking function \n terrain factor N=", N)
    legend.cost <- paste0("walking-time (", time,")")
  }

  if(funct=="ks") {
    main.title <- paste0("Least-cost corridor (cost in ", time, ")")
    sub.title <- paste0("Walking-time based on the Kondo-Seino's hiking function \n terrain factor N=", N)
    legend.cost <- paste0("walking-time (", time,")")
  }

  if(funct=="ree") {

    #set the labels to be used within the returned plot
    main.title <- "Least-cost corridor"
    sub.title <- paste0("Cost based on the relative energetic expenditure cost function \n terrain factor N=", N)
    legend.cost <- "cost"
  }

  if(funct=="hrz") {

    #set the labels to be used within the returned plot
    main.title <- "Least-cost corridor"
    sub.title <- paste0("Cost based on the Herzog's metabolic cost function \n cost in J / (Kg*m) \n terrain factor N=", N)
    legend.cost <- "metabolic cost J / (Kg*m)"
  }

  if(funct=="wcs") {

    #set the labels to be used within the returned plot
    main.title <- "Least-cost corridor"
    sub.title <- paste0("Cost based on the wheeled-vehicle critical slope cost function \ncritical slope set to ", sl.crit, " percent \n terrain factor N=", N)
    legend.cost <- "cost"
  }

  if(funct=="vl") {

    #set the labels to be used within the returned plot
    main.title <- "Least-cost corridor"
    legend.cost <- "energy expenditure cost (Megawatts)"
    if (V==0) {
      sub.title <- paste0("Cost based on the Van Leusen's metabolic energy expenditure cost function \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V is based on the Tobler on-path hiking function")
    } else {
      sub.title <- paste0("Cost based on the Van Leusen's metabolic energy expenditure cost function \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V: ", V)
    }
  }

  if(funct=="p") {

    #set the labels to be used within the returned plot
    main.title <- "Least-cost corridor"
    legend.cost <- "energy expenditure cost (Megawatts)"
    if (V==0) {
      sub.title <- paste0("Cost based on the Pandolf et al.'s metabolic energy expenditure cost function \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V is based on the Tobler on-path hiking function")
    } else {
      sub.title <- paste0("Cost based on the Pandolf et al.'s metabolic energy expenditure cost function \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V: ", V)
    }
  }

  if(funct=="pcf") {

    #set the labels to be used within the returned plot
    main.title <- "Least-cost corridor"
    legend.cost <- "energy expenditure cost (Megawatts)"
    if (V==0) {
      sub.title <- paste0("Cost based on the Pandolf et al.'s metabolic energy expenditure cost function with correction factor \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V is based on the Tobler on-path hiking function")
    } else {
      sub.title <- paste0("Cost based on the Pandolf et al.'s metabolic energy expenditure cost function with correction factor \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V: ", V)
    }
  }

  if(funct=="m") {

    #set the labels to be used within the returned plot
    main.title <- "Least-cost corridor"
    sub.title <- paste0("Cost based on the Minetti et al.'s metabolic cost function \n cost in J / (Kg*m) \n terrain factor N=", N)
    legend.cost <- "metabolic cost J / (Kg*m)"
  }

  if(funct=="ls") {

    #set the labels to be used within the returned plot
    main.title <- "Least-cost corridor"
    sub.title <- paste0("Cost based on the Llobera-Sluckin's metabolic energy expenditure cost function \n terrain factor N=", N)
    legend.cost <- "energy expenditure cost (KJ/m)"
  }

  if(funct=="alb") {

    #set the labels to be used within the returned plot
    main.title <- paste0("Least-cost corridor (cost in ", time, ")")
    sub.title <- "Walking-time based on the Alberti's modified Tobler hiking function"
    legend.cost <- paste0("walking-time (", time,")")
  }

  if(funct=="b") {

    #set the labels to be used within the returned plot
    main.title <- "Least-cost corridor"
    sub.title <- paste0("Cost based on the Bellavia's cost function \n terrain factor N=", N)
    legend.cost <- "cost"
  }

  if (funct=="ma") {

    #set the labels to be used within the returned plot
    main.title <- paste0("Least-cost corridor (cost in ", time, ")")
    sub.title <- paste0("Walking-time based on the Marin Arroyo's hiking function \n terrain factor N=", N)
    legend.cost <- paste0("walking-time (", time,")")
  }

  if(funct=="a") {

    #set the labels to be used within the returned plot
    main.title <- "Least-cost corridor"
    legend.cost <- "energy expenditure cost J / (Kg*m)"
    if (V==0) {
      sub.title <- paste0("Cost based on the Ardigo et al.'s metabolic energy expenditure cost function \n cost in J / (Kg*m) \nparameters: N: ", N, "; V is based on the Tobler on-path hiking function")
    } else {
      sub.title <- paste0("Cost based on the Ardigo et al.'s metabolic energy expenditure cost function \n cost in J / (Kg*m) \nparameters: N: ", N, "; V: ", V)
    }
  }

  if(funct=="h") {

    #set the labels to be used within the returned plot
    main.title <- "Least-cost corridor"
    sub.title <- paste0("Cost based on the Hare's metabolic cost function \n cost in cal/km \n terrain factor N=", N)
    legend.cost <- "metabolic cost cal/km"
  }

  if(funct=="e") {

    #set the labels to be used within the returned plot
    main.title <- "Least-cost corridor"
    sub.title <- paste0("Cost based on the Eastman's cost function \n terrain factor N=", N)
    legend.cost <- "cost"
  }

  #plot the corridor raster
  raster::plot(res.corridor,
               main=main.title,
               sub=sub.title,
               cex.main=0.95,
               cex.sub=0.75,
               legend.lab=legend.cost)

  #produce the ingredients for the hillshade raster
  #to be used in both the rendered plots
  slope <- raster::terrain(dtm, opt = "slope")
  aspect <- raster::terrain(dtm, opt = "aspect")
  hill <- raster::hillShade(slope, aspect, angle = 45, direction = 0)

  #add the hillshade
  raster::plot(hill,
               col = grey(0:100/100),
               legend = FALSE,
               alpha=transp,
               add=TRUE)


  #if the input dataset a was containing just 1 feature...
  if (length(a) < 2) {
  #plot the LCP from a
  raster::plot(res.a$LCPs,
               col="red",
               add=TRUE)

  #plot the LCP from b
  raster::plot(res.b$LCPs,
               add=TRUE)

  #plot point a
  raster::plot(a,
               pch=20,
               col="red",
               add=TRUE)

  #plot point b
  raster::plot(b,
               pch=20,
               col="black",
               add=TRUE)

  #add label to point a
  raster::text(sp::coordinates(a),
               labels=lab.a,
               col="red",
               pos = 4,
               cex=cex.labs)

  #add label to point b
  raster::text(sp::coordinates(b),
               labels=lab.b,
               pos = 4,
               cex=cex.labs)
  }


  #if the input dataset a was containing more than 2 features...
  if (length(a) > 2) {
    #plot all the point locations on top of the previously
    #plotted hillshade and corridor rasters
    raster::plot(a,
                 pch=20,
                 add=TRUE)
  }


  #if export is TRUE & if the input dataset 'a' was containing just 1 feature
  #export the corridor and the other data related to the two locations a and b
  if(export==TRUE & length(a) < 2){

    raster::writeRaster(res.corridor, paste0("LCcorridor_", funct), format="GTiff")
    raster::writeRaster(res.a$accumulated.cost.raster, paste0("Accum_cost_surf_a_", funct), format="GTiff")
    raster::writeRaster(res.b$accumulated.cost.raster, paste0("Accum_cost_surf_b_", funct), format="GTiff")
    rgdal::writeOGR(res.a$LCPs, ".", paste0("lcp_a_to_b_", funct), driver="ESRI Shapefile")
    rgdal::writeOGR(res.b$LCPs, ".", paste0("lcp_b_to_a_", funct), driver="ESRI Shapefile")

    }

  #if export is TRUE & if the input dataset 'a' was containing more than 2 features
  #export only the corridor
  if(export==TRUE & length(a) > 2){
      raster::writeRaster(res.corridor, paste0("LCcorridor_", funct), format="GTiff")
    }


  #if no DTM was provided (i.e., if 'studyplot' is not NULL), export the downloaded DTM as a raster file
  if(export==TRUE & is.null(studyplot)==FALSE){
    raster::writeRaster(dtm, "dtm", format="GTiff")
  }

  #if the input dataset a was containing just 1 feature...
  if (length(a) < 2) {
  rslt <- list("dtm"=dtm,
                  "lc.corridor"=res.corridor,
                  "lcp_a_to_b"=res.a$LCPs,
                  "lcp_b_to_a"=res.b$LCPs,
                  "accum_cost_surf_a"=res.a$accumulated.cost.raster,
                  "accum_cost_surf_b"=res.b$accumulated.cost.raster)
  }

  #if the input dataset a was containing more than 2 features...
  if (length(a) > 2) {
    rslt<- list("dtm"=dtm,
                "lc.corridor"=res.corridor,
                "corridors"=corridors,
                "locations.matrix"=pairw.list)
  }

  results <- rslt
  }
