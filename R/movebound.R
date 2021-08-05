#' R function for calculating slope-dependant walking cost boundary(ies) around point location(s)
#'
#' The function provides the facility to calculate walking cost boundary(ies) around one or more point locations. Rationale: while \code{\link{movecost}}
#' can calculate and render an accumulated cost surface and corresponding isolines around a point location, the user(s) might want to
#' calculate and plot a boundary (or boundaries) corresponding to a specific walking cost limit around one or more locations,
#' either in terms of walking time or energy expenditure.\cr
#'
#' The following example uses in-built datasets and calculates 45minute boundaries around three locations close to Mt Etna (Sicily, Italy), using the
#' Tobler's off-path hiking function (note: elevation data are acquired online for the area enclosed by the polygon fed via the
#' 'studyplot' parameter):\cr
#'
#' result <- movebound(origin=Etna_end_location, cont.value=45, time="m", cont.lab = TRUE, funct="tofp", studyplot = Etna_boundary)\cr
#'
#' The function just requires an input DTM and a dataset ('SpatialPointsDataFrame' class) containing at least one point location.
#' If a DTM is not provided, 'movebound()' will download elevation data from online sources (see \code{\link{movecost}} for more details).
#' Under the hood, 'movebound()' relies on the \code{\link{movecost}} function and implements the same
#' cost functions: see the help documentation of 'movecost()' for further information.\cr
#'
#' 'movebound()' produces a plot representing the input DTM overlaid by an hillshade raster, whose transparency can be adjusted uding
#' the 'transp' parameter. On the rendered plot, the calculated isoline(s) is displayed and the label(s) representing the cost limit can be
#' activated or deactivated using the 'cont.lab' parameter. The function also returns the isoline(s) ('SpatialLinesDataFrame' class) corresponding
#' to the selected accumulated cost limit (see 'Value' below). The isoline(s) can be exported as shapefile using the 'export' parameter. \cr
#'
#' @param dtm Digital Terrain Model (RasterLayer class); if not provided, elevation data will be acquired online for the area enclosed by the 'studyplot' parameter (see \code{\link{movecost}}).
#' @param origin location(s) around which the boundary(ies) is calculated (SpatialPointsDataFrame class).
#' @param studyplot polygon (SpatialPolygonDataFrame class) representing the study area for which online elevation data are aquired (see \code{\link{movecost}}); NULL is default.
#' @param funct cost function to be used (for details on each of the following, see \code{\link{movecost}}):\cr
#' \strong{t} (default) uses the on-path Tobler's hiking function;\cr
#' \strong{tofp} uses the off-path Tobler's hiking function;\cr
#' \strong{mp} uses the Marquez-Perez et al.'s modified Tobler's function;\cr
#' \strong{icmonp} uses the Irmischer-Clarke's hiking function (male, on-path);\cr
#' \strong{icmoffp} uses the Irmischer-Clarke's hiking function (male, off-path);\cr
#' \strong{icfonp} uses the Irmischer-Clarke's hiking function (female, on-path);\cr
#' \strong{icfoffp} uses the Irmischer-Clarke's hiking function (female, off-path);\cr
#' \strong{ug} uses the Uriarte Gonzalez's slope-dependant walking-time cost function;\cr
#' \strong{alb} uses the Alberti's Tobler hiking function modified for pastoral foraging excursions;\cr
#' \strong{gkrs} uses the Garmy, Kaddouri, Rozenblat, and Schneider's hiking function;\cr
#' \strong{r} uses the Rees' hiking function;\cr
#' \strong{ks} uses the Kondo-Seino's hiking function;\cr
#' \strong{ree} uses the relative energetic expenditure cost function;\cr
#' \strong{hrz} uses the Herzog's metabolic cost function;\cr
#' \strong{wcs} uses the wheeled-vehicle critical slope cost function;\cr
#' \strong{p} uses the Pandolf et al.'s metabolic energy expenditure cost function;\cr
#' \strong{vl} uses the Van Leusen's metabolic energy expenditure cost function;\cr
#' \strong{ls} uses the Llobera-Sluckin's metabolic energy expenditure cost function.\cr
#' \strong{b} uses the Bellavia's cost function.\cr
#' @param time time-unit expressed by the isoline(s) if Tobler's and other time-related cost functions are used;
#' 'h' for hour, 'm' for minutes.
#' @param move number of directions in which cells are connected: 4 (rook's case), 8 (queen's case), 16 (knight and one-cell queen moves; default).
#' @param cont.value cost value represented by the calculated isoline(s).
#' @param cogn.slp  TRUE or FALSE (default) if the user wants or does not want the 'cognitive slope' to be used in place of the real slope (see \code{\link{movecost}}).
#' @param sl.crit critical slope (in percent), typically in the range 8-16 (10 by default) (used by the wheeled-vehicle cost function; see \code{\link{movecost}}).
#' @param W walker's body weight (in Kg; 70 by default; used by the Pandolf's and Van Leusen's cost function; see \code{\link{movecost}}).
#' @param L carried load weight (in Kg; 0 by default; used by the Pandolf's and Van Leusen's cost function; see \code{\link{movecost}}).
#' @param N coefficient representing ease of movement (1 by default) (see \code{\link{movecost}}).
#' @param V speed in m/s (1.2 by default) (used by the Pandolf's and Van Leusen's cost function; if set to 0, it is internally worked out on the basis of Tobler on-path hiking function;see \code{\link{movecost}}).
#' @param z zoom level for the elevation data downloaded from online sources (from 0 to 15; 9 by default) (see \code{\link{movecost}} and \code{\link[elevatr]{get_elev_raster}}).
#' @param cont.lab TRUE (default) or FALSE if the usuer wants or does not want labels to be attached to the isolines.
#' @param transp set the transparency of the hillshade raster that is plotted over the DTM (0.5 by default).
#' @param export TRUE or FALSE (default) if the user wants or does not want the isoline(s) to be exported; if TRUE, the isoline(s) will be exported as a shapefile;
#'  the exported file will bear a suffix corresponding to the cost function selected by the user. The DTM is exported only if it was not provided by the user and downloaded
#'  by the function from online sources.
#'
#' @return The function returns a list storing the following components \itemize{
##'  \item{dtm: }{Digital Terrain Model ('RasterLayer' class)}
##'  \item{isolines: }{contour line(s) representing the selected walking cost limit ('SpatialLinesDataFrame' class)}
##' }
##'
#' @keywords movebound
#' @export
#' @importFrom raster ncell mask crop stack cellStats raster hillShade terrain
#' @importFrom elevatr get_elev_raster
#' @importFrom grDevices terrain.colors topo.colors grey
#'
#' @examples
#' # load a sample Digital Terrain Model
#' volc <- raster::raster(system.file("external/maungawhau.grd", package="gdistance"))
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
movebound <- function (dtm=NULL, origin, studyplot=NULL, funct="t", time="h", move=16, cont.value, cogn.slp=FALSE, sl.crit=10, W=70, L=0, N=1, V=1.2, z=9, cont.lab=TRUE, transp=0.5, export=FALSE){

  #if no dtm is provided
  if (is.null(dtm)==TRUE) {
    #get the elvation data using the elevatr's get_elev_raster() function, using the studyplot dataset (SpatialPolygonDataFrame)
    #to select the area whose elevation data are to be downloaded;
    #z sets the resolution of the elevation datataset
    elev.data <- elevatr::get_elev_raster(studyplot, z = z, verbose=FALSE, override_size_check = TRUE)
    #crop the elevation dataset to the exact boundary of the studyplot dataset
    dtm <- raster::crop(elev.data, studyplot)
  }

  result <- movecost::movecost(dtm=dtm, origin=origin, studyplot=studyplot, funct=funct, time=time, move=move, cogn.slp=cogn.slp, sl.crit=sl.crit, W=W, L=L, N=N, V=V, z=z, graph.out=FALSE)$accumulated.cost.raster

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


  #calculate and store the contours as a SpatialLinesDataFrame
  isolines <- raster::rasterToContour(result, levels=cont.value)

  #produce the ingredients for the hillshade raster
  slope <- raster::terrain(dtm, opt = "slope")
  aspect <- raster::terrain(dtm, opt = "aspect")
  hill <- raster::hillShade(slope, aspect, angle = 45, direction = 0)

  #plot the DTM
  raster::plot(dtm,
               main=main.title,
               sub=sub.title,
               cex.main=0.95,
               cex.sub=0.75,
               legend.lab="Elevation (masl)")

  #plot the hillshade
  raster::plot(hill,
               col = grey(0:100/100),
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


  #if export is TRUE, export the isolines as a shapefile
  if(export==TRUE){
    rgdal::writeOGR(isolines, ".", paste0("isolines_", funct), driver="ESRI Shapefile")
  }


  #if no DTM was provided (i.e., if 'studyplot' is not NULL), export the downloaded DTM as a raster file
  if(export==TRUE & is.null(studyplot)==FALSE){
    raster::writeRaster(dtm, "dtm", format="GTiff")
  }

  #store the results in a list
  results <- list("dtm"=dtm,
                   "isolines" = isolines)

}
