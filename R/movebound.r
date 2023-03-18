#' R function for calculating slope-dependant walking cost boundary(ies) around point location(s)
#'
#' The function provides the facility to calculate walking cost boundary(ies) around one or more point locations. Rationale: while \code{\link{movecost}}
#' can calculate and render an accumulated cost surface and corresponding isolines around a point location, the user(s) might want to
#' calculate and plot a boundary (or boundaries) corresponding to a specific walking cost limit around one or more locations,
#' either in terms of walking time or energy expenditure.\cr
#' Visit this \href{https://drive.google.com/file/d/1gLDrkZFh1b_glzCEqKdkPrer72JJ9Ffa/view?usp=sharing}{LINK} to access the package's vignette.\cr
#'
#' The function just requires an input DTM and a dataset ('SpatialPointsDataFrame' class) containing at least one point location.
#' If a DTM is not provided, \code{movebound()} will download elevation data from online sources (see \code{\link{movecost}} for more details).
#' Under the hood, \code{movebound()} relies on the \code{movecost()} function and implements the same
#' cost functions: see the help documentation of \code{movecost()} for further information.\cr
#'
#' The following example uses in-built datasets and calculates 45-minute boundaries around three locations close to Mt Etna (Sicily, Italy), using the
#' Tobler's off-path hiking function (note: elevation data are acquired online for the area enclosed by the polygon fed via the
#' 'studyplot' parameter):\cr
#'
#' result <- movebound(origin=Etna_end_location, cont.value=45, time="m", cont.lab = TRUE, funct="tofp", studyplot = Etna_boundary, add.geom=TRUE)\cr
#'
#' Note that by setting the parameter \code{add.geom} to \code{TRUE}, the function calculates the area enclosed by the boundary represented by each
#' calculated isoline. Needless to say, the unit of measure is the one used by the input layers' coordinate system. The value(s) of the area
#' will be appended as a new variables to a copy of the input 'origin' dataset. The area can only be calculated if the
#' isolines are "complete" and not truncated (i.e., if they do not meet the end of the study area for instance). Therefore, before using this option,
#' the user may want to be sure that all the isolines are actual loops.\cr
#'
#' With reference to the above example, the area of the three 45-minutes boundaries can be retrieved typing what follows:\cr
#'
#' result$origin_w_isolines_geom$area \cr
#'
#' It will return:\cr
#' 17857994 20428575 9172688\cr
#'
#' that are the values of the area of each 45-minute boundary in square meter.\cr
#'Needless to say, if we want to convert to square km we can just:\cr
#' result$origin_w_isolines_geom$area/1000000 \cr
#'
#' which gives\cr
#' 17.857994 20.428575  9.172688 \cr
#'
#' \code{movebound()} produces a plot representing the input DTM overlaid by a slopeshade raster, whose transparency can be adjusted using
#' the 'transp' parameter. On the rendered plot, the calculated isoline(s) is displayed and the label(s) representing the cost limit can be
#' activated or deactivated using the 'cont.lab' parameter. The function also returns the isoline(s) ('SpatialLinesDataFrame' class) corresponding
#' to the selected accumulated cost limit and the copy of the 'origin' dataset (storing information about the boundaries' geometry) (see 'Value' below).
#' The isoline(s) and the copy of the 'origin' dataset can be exported as shapefile by setting the \code{export} parameter to \code{TRUE}. \cr
#'
#' @param dtm Digital Terrain Model (RasterLayer class); if not provided, elevation data will be acquired online for the area enclosed by the 'studyplot' parameter (see \code{\link{movecost}}).
#' @param origin location(s) around which the boundary(ies) is calculated (SpatialPointsDataFrame class).
#' @param studyplot polygon (SpatialPolygonDataFrame class) representing the study area for which online elevation data are acquired (see \code{\link{movecost}}); NULL is default.
#' @param barrier area where the movement is inhibited (SpatialLineDataFrame or SpatialPolygonDataFrame class) (see \code{\link{movecost}}.
#' @param plot.barrier TRUE or FALSE (default) if the user wants or does not want the barrier to be plotted (see \code{\link{movecost}}).
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
#' \strong{trp} uses the Tripcevich's hiking function;\cr
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
#' \strong{a} uses the Ardigo et al.'s metabolic energy expenditure cost function (for all the mentioned cost functions;\cr
#' \strong{h} uses the Hare's metabolic energy expenditure cost function (for all the mentioned cost functions, see \code{\link{movecost}}).\cr
#' @param time time-unit expressed by the isoline(s) if Tobler's and other time-related cost functions are used; h' for hour, 'm' for minutes.
#' @param move number of directions in which cells are connected: 4 (rook's case), 8 (queen's case), 16 (knight and one-cell queen moves; default).
#' @param field value assigned to the cells coincidinng with the barrier (0 by default) (see \code{\link{movecost}}.
#' @param cont.value cost value represented by the calculated isoline(s) (NULL by default); if no value is supplied, it is set to 1/10 of the range of values of the accumulated cost surface.
#' @param cogn.slp  TRUE or FALSE (default) if the user wants or does not want the 'cognitive slope' to be used in place of the real slope (see \code{\link{movecost}}).
#' @param sl.crit critical slope (in percent), typically in the range 8-16 (10 by default) (used by the wheeled-vehicle cost function; see \code{\link{movecost}}).
#' @param W walker's body weight (in Kg; 70 by default; used by the Pandolf's and Van Leusen's cost function; see \code{\link{movecost}}).
#' @param L carried load weight (in Kg; 0 by default; used by the Pandolf's and Van Leusen's cost function; see \code{\link{movecost}}).
#' @param N coefficient representing ease of movement (1 by default) (see \code{\link{movecost}}).
#' @param V speed in m/s (1.2 by default) (used by the Pandolf et al.'s, Pandolf et al.s with correction factor, Van Leusen's, and Ardigo et al.'s cost function; if set to 0, it is internally worked out on the basis of Tobler on-path hiking function (see \code{\link{movecost}}).
#' @param z zoom level for the elevation data downloaded from online sources (from 0 to 15; 9 by default) (see \code{\link{movecost}} and \code{\link[elevatr]{get_elev_raster}}).
#' @param cont.lab TRUE (default) or FALSE if the usuer wants or does not want labels to be attached to the isolines.
#' @param transp set the transparency of the slopeshade raster that is plotted over the DTM (0.5 by default).
#' @param add.geom TRUE or FALSE (default) if the user wants or does not want the area enclosed by each isolines to be calculated (see Details).
#' @param export TRUE or FALSE (default) if the user wants or does not want the isoline(s) and the copy of the input 'origin' dataset (storing boundaries' geometry information)
#'  to be exported; if TRUE, they will be exported as a shapefile; the exported file will bear a suffix corresponding to the cost function selected by the user.
#'  The DTM is exported only if it was not provided by the user and downloaded by the function from online sources.
#'
#' @return The function returns a list storing the following components \itemize{
##'  \item{dtm: }{Digital Terrain Model ('RasterLayer' class)}
##'  \item{isolines: }{contour line(s) representing the selected cost limit ('SpatialLinesDataFrame' class)}
##'  \item{origin_w_isolines_geom: }{copy of the input origin location(s) dataset with a new variable ('area') storing the
##'    area values of the boundary calculated around each location}
##' }
##'
#' @keywords movebound
#' @export
#' @importFrom raster ncell mask crop stack cellStats raster terrain crs
#' @importFrom terra writeVector vect
#' @importFrom elevatr get_elev_raster
#' @importFrom grDevices terrain.colors topo.colors grey
#' @importFrom sf st_area st_as_sf
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
#' # calculate the 5minute walking time boundary around a location
#' # using the Tobler's off-path hiking function
#'
#' result <- movebound(dtm=volc, origin=volc.loc, funct="tofp", move=8, time="m", cont.val=5)
#'
#'
#' # same as above, but around multiple locations; contours' labels are turned off
#'
#' result <- movebound(dtm=volc, origin=destin.loc, funct="tofp", move=8, time="m",
#' cont.val=2, cont.lab=FALSE)
#'
#'
#' @seealso \code{\link{movecost}}
#'
#'
movebound <- function (dtm=NULL, origin, studyplot=NULL, barrier=NULL, plot.barrier=FALSE, funct="t", time="h", move=16, field=0, cont.value=NULL, cogn.slp=FALSE, sl.crit=10, W=70, L=0, N=1, V=1.2, z=9, cont.lab=TRUE, transp=0.5, add.geom=FALSE, export=FALSE){

  #if no dtm is provided
  if (is.null(dtm)==TRUE) {
    #get the elvation data using the elevatr's get_elev_raster() function, using the studyplot dataset (SpatialPolygonDataFrame)
    #to select the area whose elevation data are to be downloaded;
    #z sets the resolution of the elevation datataset
    elev.data <- elevatr::get_elev_raster(studyplot, z = z, verbose=FALSE, override_size_check = TRUE)
    #crop the elevation dataset to the exact boundary of the studyplot dataset
    dtm <- raster::crop(elev.data, studyplot)
  }

  result <- movecost::movecost(dtm=dtm, origin=origin, studyplot=studyplot, barrier=barrier, funct=funct, time=time, move=move, field=field, cogn.slp=cogn.slp, sl.crit=sl.crit, W=W, L=L, N=N, V=V, z=z, graph.out=FALSE)$accumulated.cost.raster

  #if no contour value is entered, set it to one tenth of the range of the values of the final accumulated cost surface
  if(is.null(cont.value)==TRUE){
    #exclude the inf values from the calculation
    cont.value <- round((max(result[][is.finite(result[])]) - min(result[][is.finite(result[])])) / 10,2)
  }

  #define the appropriate text to be used for subsequent plotting
  if (funct=="t") {

    #set the labels to be used within the returned plot
    main.title <- paste0("Boundary(ies) corresponding to ", cont.value, time, " walking limit")
    sub.title <- paste0("Walking-time based on the Tobler's on-path hiking function \n terrain factor N=", N)
    legend.cost <- paste0("walking-time (", time,")")
  }

  if (funct=="tofp") {

    #set the labels to be used within the returned plot
    main.title <- paste0("Boundary(ies) corresponding to ", cont.value, time, " walking limit")
    sub.title <- "Walking-time based on the Tobler's off-path hiking function"
    legend.cost <- paste0("walking-time (", time,")")
  }

  if(funct=="mp") {

    #set the labels to be used within the returned plot
    main.title <-  paste0("Boundary(ies) corresponding to ", cont.value, time, " walking limit")
    sub.title <- paste0("Walking-time based on the Marquez-Perez et al.'s modified Tobler hiking function \n terrain factor N=", N)
    legend.cost <- paste0("walking-time (", time,")")
  }

  if(funct=="icmonp") {

    #set the labels to be used within the returned plot
    main.title <- paste0("Boundary(ies) corresponding to ", cont.value, time, " walking limit")
    sub.title <- paste0("Walking-time based on the (male, on-path) Irmischer-Clarke's hiking function \n terrain factor N=", N)
    legend.cost <- paste0("walking-time (", time,")")
  }

  if(funct=="icmoffp") {

    #set the labels to be used within the returned plot
    main.title <-  paste0("Boundary(ies) corresponding to ", cont.value, time, " walking limit")
    sub.title <- "Walking-time based on the (male, off-path) Irmischer-Clarke's hiking function"
    legend.cost <- paste0("walking-time (", time,")")
  }

  if(funct=="icfonp") {

    #set the labels to be used within the returned plot
    main.title <-  paste0("Boundary(ies) corresponding to ", cont.value, time, " walking limit")
    sub.title <- paste0("Walking-time based on the (female, on-path) Irmischer-Clarke's hiking function \n terrain factor N=", N)
    legend.cost <- paste0("walking-time (", time,")")
  }

  if(funct=="icfoffp") {
    main.title <-  paste0("Boundary(ies) corresponding to ", cont.value, time, " walking limit")
    sub.title <- "Walking-time based on the (female, off-path) Irmischer-Clarke's hiking function"
    legend.cost <- paste0("walking-time (", time,")")
  }

  if(funct=="gkrs") {
    main.title <-  paste0("Boundary(ies) corresponding to ", cont.value, time, " walking limit")
    sub.title <- paste0("Walking-time based on the Garmy et al.'s hiking function \n terrain factor N=", N)
    legend.cost <- paste0("walking-time (", time,")")
  }

  if(funct=="r") {
    main.title <-  paste0("Boundary(ies) corresponding to ", cont.value, time, " walking limit")
    sub.title <- paste0("Walking-time based on the Rees' hiking function \n terrain factor N=", N)
    legend.cost <- paste0("walking-time (", time,")")
  }

  if(funct=="ug") {

    #set the labels to be used within the returned plot
    main.title <- paste0("Boundary(ies) corresponding to ", cont.value, time, " walking limit")
    sub.title <- paste0("Walking-time based on the Uriarte Gonzalez's hiking function \n terrain factor N=", N)
    legend.cost <- paste0("walking-time (", time,")")
  }

  if(funct=="ks") {

    #set the labels to be used within the returned plot
    main.title <- paste0("Boundary(ies) corresponding to ", cont.value, time, " walking limit")
    sub.title <- paste0("Walking-time based on the Kondo-Seino's hiking function \n terrain factor N=", N)
    legend.cost <- paste0("walking-time (", time,")")
  }

  if(funct=="ree") {

    #set the labels to be used within the returned plot
    main.title <-  paste0("Boundary(ies) corresponding to ", cont.value, " cost lmit")
    sub.title <- paste0("Cost based on the relative energetic expenditure cost function \n terrain factor N=", N)
    legend.cost <- "cost"
  }

  if(funct=="hrz") {

    #set the labels to be used within the returned plot
    main.title <-  paste0("Boundary(ies) corresponding to ", cont.value, " J/(Kg*m) cost limit")
    sub.title <- paste0("Cost based on the Herzog's metabolic cost function \n cost in J / (Kg*m) \n terrain factor N=", N)
    legend.cost <- "metabolic cost J / (Kg*m)"
  }

  if(funct=="wcs") {

    #set the labels to be used within the returned plot
    main.title <-  paste0("Boundary(ies) corresponding to ", cont.value, " cost limit")
    sub.title <- paste0("Cost based on the wheeled-vehicle critical slope cost function \ncritical slope set to ", sl.crit, " percent \n terrain factor N=", N)
    legend.cost <- "cost"
  }

  if(funct=="vl") {

    #set the labels to be used within the returned plot
    main.title <-  paste0("Boundary(ies) corresponding to ", cont.value, " Megawatts cost limit")
    legend.cost <- "energy expenditure cost (Megawatts)"
    if (V==0) {
      sub.title <- paste0("Cost based on the Van Leusen's metabolic energy expenditure cost function \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V is based on the Tobler on-path hiking function")
    } else {
      sub.title <- paste0("Cost based on the Van Leusen's metabolic energy expenditure cost function \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V: ", V)
    }
  }

  if(funct=="p") {

    #set the labels to be used within the returned plot
    main.title <- paste0("Boundary(ies) corresponding to ", cont.value, " Megawatts cost limit")
    legend.cost <- "energy expenditure cost (Megawatts)"
    if (V==0) {
      sub.title <- paste0("Cost based on the Pandolf et al.'s metabolic energy expenditure cost function \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V is based on the Tobler on-path hiking function")
    } else {
      sub.title <- paste0("Cost based on the Pandolf et al.'s metabolic energy expenditure cost function \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V: ", V)
    }
  }

  if(funct=="pcf") {

    #set the labels to be used within the returned plot
    main.title <- paste0("Boundary(ies) corresponding to ", cont.value, " Megawatts cost limit")
    legend.cost <- "energy expenditure cost (Megawatts)"
    if (V==0) {
      sub.title <- paste0("Cost based on the Pandolf et al.'s metabolic energy expenditure cost function with correction factor \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V is based on the Tobler on-path hiking function")
    } else {
      sub.title <- paste0("Cost based on the Pandolf et al.'s metabolic energy expenditure cost function with correction factor \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V: ", V)
    }
  }

  if(funct=="ls") {

    #set the labels to be used within the returned plot
    main.title <- paste0("Boundary(ies) corresponding to ", cont.value, " KJ/m cost limit")
    sub.title <- paste0("Cost based on the Llobera-Sluckin's metabolic energy expenditure cost function \n terrain factor N=", N)
    legend.cost <- "energy expenditure cost (KJ/m)"
  }

  if(funct=="alb") {

    #set the labels to be used within the returned plot
    main.title <- paste0("Boundary(ies) corresponding to ", cont.value, time, " walking limit")
    sub.title <- "Walking-time based on the Alberti's modified Tobler hiking function"
    legend.cost <- paste0("walking-time (", time,")")
  }

  if(funct=="b") {

    #set the labels to be used within the returned plot
    main.title <- paste0("Boundary(ies) corresponding to ", cont.value, " cost limit")
    sub.title <- paste0("Cost based on the Bellavia's cost function \n terrain factor N=", N)
    legend.cost <- "cost"
  }

  if(funct=="m") {

    #set the labels to be used within the returned plot
    main.title <-  paste0("Boundary(ies) corresponding to ", cont.value, " J/(Kg*m) cost limit")
    sub.title <- paste0("Cost based on the Minetti et al.'s metabolic cost function \n cost in J / (Kg*m) \n terrain factor N=", N)
    legend.cost <- "metabolic cost J / (Kg*m)"
  }

  if (funct=="ma") {

    #set the labels to be used within the returned plot
    main.title <- paste0("Boundary(ies) corresponding to ", cont.value, time, " walking limit")
    sub.title <- paste0("Walking-time based on the Marin Arroyo's hiking function \n terrain factor N=", N)
    legend.cost <- paste0("walking-time (", time,")")
  }

  if(funct=="a") {

    #set the labels to be used within the returned plot
    main.title <- paste0("Boundary(ies) corresponding to ", cont.value, " J/(Kg*m) cost limi")
    legend.cost <- "metabolic cost J / (Kg*m)"
    if (V==0) {
      sub.title <- paste0("Cost based on the Ardigo et al.'s metabolic cost function \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V is based on the Tobler on-path hiking function")
    } else {
      sub.title <- paste0("Cost based on the Ardigo et al.'s metabolic cost function \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V: ", V)
    }
  }

  if(funct=="h") {

    #set the labels to be used within the returned plot
    main.title <-  paste0("Boundary(ies) corresponding to ", cont.value, " cal/km cost limit")
    sub.title <- paste0("Cost based on the Hare's metabolic cost function \n cost in cal/km \n terrain factor N=", N)
    legend.cost <- "metabolic cost cal/km"
  }

  if(funct=="e") {

    #set the labels to be used within the returned plot
    main.title <- paste0("Boundary(ies) corresponding to ", cont.value, " cost limit")
    sub.title <- paste0("Cost based on the Eastman's cost function \n terrain factor N=", N)
    legend.cost <- "cost"
  }

  if (funct=="trp") {

    #set the labels to be used within the returned plot
    main.title <- paste0("Boundary(ies) corresponding to ", cont.value, time, " walking limit")
    sub.title <- paste0("Walking-time based on the Tripcevich's hiking function \n terrain factor N=", N)
    legend.cost <- paste0("walking-time (", time,")")
  }

  #extract and store the contours as a SpatialLinesDataFrame
  isolines <- raster::rasterToContour(result, levels=cont.value)

  #if add.geom is set to TRUE
  if(add.geom==TRUE) {
    #extract the coordinates of the isolines
    contour_coord <- sp::coordinates(isolines)

    #create two empty vectors where isolines' geometry will be stored
    contour_area <- as.numeric()

    #loop through the number of origin (i.e., the number of extraced isolines)
    for (i in 1:length(origin)) {
      #trasform the isolines into 'sp' polygon objects that will be read by some sf's functions later on
      p <- sp::Polygon(contour_coord[[1]][[i]])
      ps <- sp::Polygons(list(p),1)
      sps <- sp::SpatialPolygons(list(ps))

      #give the polygons the coordinate system of the input dtm
      #raster::crs(sps) <- raster::crs(dtm)

      #use functions out from sf to calculate the area and perimeter of the isolines
      contour_area[i] <- sf::st_area(st_as_sf(sps))
    }
    #create a copy of the 'origin' dataset
    origin_w_geom <- origin

    #the same as above for the area values
    origin_w_geom$area <- contour_area
  }

  #produce the ingredient for the slopeshade raster
  slope <- raster::terrain(dtm, opt = "slope")

  #plot the DTM
  raster::plot(dtm,
               main=main.title,
               sub=sub.title,
               cex.main=0.95,
               cex.sub=0.75,
               legend.lab="Elevation (masl)")

  #plot the slopeshade
  raster::plot(slope,
               col = rev(grey(0:100/100)),
               legend = FALSE,
               alpha=transp,
               add=TRUE)

  #add the origin(s)
  raster:: plot(origin,
                pch=20,
                add=TRUE)

  #add the contours
  raster::contour(result,
                  add=TRUE,
                  levels=cont.value,
                  drawlabels = cont.lab)

  #if the barrier is provided AND if plot.barrier is TRUE, add the barrier
  if(is.null(barrier)==FALSE & plot.barrier==TRUE) {
    raster::plot(barrier, col="blue", add=TRUE)
  }

  #if add.geom is TRUE...
  if(add.geom==TRUE) {
    #the object origin_w_geom is "confirmed" so that can be returned
    origin_w_geom <- origin_w_geom
  } else {
    #otherwise, origin_w_gem is set to NULL
    origin_w_geom <- NULL
  }

  #if export is TRUE, export the isolines and the copy of the origin dataset as shapefiles
  if(export==TRUE){
    terra::writeVector(vect(isolines), filename=paste0("isolines_", funct), filetype="ESRI Shapefile")
    if(is.null(origin_w_geom)==FALSE){
      terra::writeVector(vect(origin_w_geom), filename=paste0("origin_w_geom_", funct), filetype="ESRI Shapefile")
    }
  }

  #if no DTM was provided (i.e., if 'studyplot' is not NULL), export the downloaded DTM as a raster file
  if(export==TRUE & is.null(studyplot)==FALSE){
    raster::writeRaster(dtm, "dtm", format="GTiff")
  }

  #store the results in a list
  results <- list("dtm"=dtm,
                  "isolines" = isolines,
                  "origin_w_isolines_geom" = origin_w_geom)

}
