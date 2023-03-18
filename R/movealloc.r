#' R function for calculating slope-dependant walking-cost allocation to origins
#'
#' The function provides the facility to carry out a cost allocation analysis. Given a number of origin locations,
#' a cost allocation raster is produced; each cell of the cost allocation raster is given an integer indicating to which origin
#' a cell is closer in terms of cost. Needless to say, the cost can be conceptualized in terms of either walking time or
#' energy expenditure, and is function of the terrain slope.\cr
#' Visit this \href{https://drive.google.com/file/d/1gLDrkZFh1b_glzCEqKdkPrer72JJ9Ffa/view?usp=sharing}{LINK} to access the package's vignette.\cr
#'
#' The function requires an input DTM ('RasterLayer' class) and a dataset ('SpatialPointsDataFrame' class) containing the
#' origin locations. If a DTM is not provided, \code{movealloc()} will download elevation data from online sources (see \code{\link{movecost}} for more details).
#' Under the hood, \code{movealloc()} relies on the \code{movecost()} function and implements the same cost functions:
#' see the help documentation of \code{movecost()} for further information.\cr
#'
#' Internally, what \code{movealloc()} does is producing an accumulated cost surface around each individual origin location; those
#' accumulated cost surfaces are then stacked together, and then the function looks at each pixel in the stack of surfaces and
#' returns 1 if the first stacked surface has the smallest pixel value, or 2 if the second stacked surface has the smallest pixel value, and so on for bigger stacks.\cr
#'
#' \code{movealloc()} produces a plot featuring a slopeshade image that is overlaid by the cost allocation raster and by a polygon layer
#' where each polygon represents the limits of each allocation zone. A legend can be optionally added to the plot via the \code{leg.alloc} parameter (FALSE by default).
#' Isolines (i.e., contour lines) around each origin location can be optionally plotted via the 'isolines' parameter (FALSE by default).\cr
#'
#' The DTM, the cost allocation raster, the cost allocation polygons, and the isolines (if requested by the user by setting the  \code{isolines} parameter to TRUE), can
#' be exported by setting the 'export' parameter to TRUE. All the exported files (excluding the DTM) will bear a suffix corresponding to the cost function selected by the user.
#'
#'
#' @param dtm Digital Terrain Model (RasterLayer class); if not provided, elevation data will be acquired online for the area enclosed by the 'studyplot' parameter (see \code{\link{movecost}}).
#' @param origin locations (two at least) in relation to which the cost allocation is carried out (SpatialPointsDataFrame class).
#' @param studyplot polygon (SpatialPolygonDataFrame class) representing the study area for which online elevation data are acquired (see \code{\link{movecost}}); NULL is default.
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
#' \strong{a} uses the Ardigo et al.'s metabolic energy expenditure cost function;\cr
#' \strong{h} uses the Hare's metabolic energy expenditure cost function (for all the mentioned cost functions, see \code{\link{movecost}}).\cr
#' @param time time-unit expressed by the isoline(s) if Tobler's and other time-related cost functions are used; h' for hour, 'm' for minutes.
#' @param move number of directions in which cells are connected: 4 (rook's case), 8 (queen's case), 16 (knight and one-cell queen moves; default).
#' @param cogn.slp  TRUE or FALSE (default) if the user wants or does not want the 'cognitive slope' to be used in place of the real slope (see \code{\link{movecost}}).
#' @param sl.crit critical slope (in percent), typically in the range 8-16 (10 by default) (used by the wheeled-vehicle cost function; see \code{\link{movecost}}).
#' @param W walker's body weight (in Kg; 70 by default; used by the Pandolf's and Van Leusen's cost function; see \code{\link{movecost}}).
#' @param L carried load weight (in Kg; 0 by default; used by the Pandolf's and Van Leusen's cost function; see \code{\link{movecost}}).
#' @param N coefficient representing ease of movement (1 by default) (see \code{\link{movecost}}).
#' @param V speed in m/s (1.2 by default) (used by the Pandolf et al.'s, Pandolf et al.s with correction factor, Van Leusen's, and Ardigo et al.'s cost function; if set to 0, it is internally worked out on the basis of Tobler on-path hiking function (see \code{\link{movecost}}).
#' @param z zoom level for the elevation data downloaded from online sources (from 0 to 15; 9 by default) (see \code{\link{movecost}} and \code{\link[elevatr]{get_elev_raster}}).
#' @param isolines TRUE or FALSE (default) is the user wants or does not want cost isolines/contours around the origins to be calculated and plotted.
#' @param breaks contours' (i.e., isolines') interval; if no value is supplied, the interval is set by default to 1/10 of the range of values of the cost surface accumulated around the origins.
#' @param cont.lab if set to TRUE (default) display the labels of the cost contours.
#' @param cex.breaks set the size of the labels attached to the cost contours (0.6 by default).
#' @param leg.alloc if set to TRUE, display the legend in the plotted cost allocation raster; FALSE by default.
#' @param leg.pos set the position of the legend in the plotted cost allocation raster; 'topright' by default (other options: "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right", "center").
#' @param cex.leg set the size of the labels used in the legend displayed in the plotted allocation raster (0.75 by default).
#' @param transp set the transparency of the slopeshade raster that is plotted over the cost allocation raster (0.5 by default).
#' @param export TRUE or FALSE (default) if the user wants or does not want the output to be exported; if TRUE, the isolines (i.e. the contours) and the allocation boundaries will be exported as a shapefile;
#' the cost allocation raster will be exported as 'GeoTiff'; the DTM is exported only if it was not provided by the user and downloaded by the function from online sources; all the exported files (excluding the DTM)
#' will bear a suffix corresponding to the cost function selected by the user.
#'
#' @return The function returns a list storing the following components \itemize{
##'  \item{dtm: }{Digital Terrain Model ('RasterLayer' class)}
##'  \item{cost.allocation.raster: }{raster of the cost allocation ('RasterLayer' class)}
##'  \item{isolines: }{contour lines representing the accumulated cost around the origins ('SpatialLinesDataFrame' class); returned if the 'isolines' parameter is set to TRUE}
##'  \item{alloc.boundaries: }{polygons representing the allocation zones ('SpatialPolygonsDataFrame' class)}
##' }
##'
#' @keywords movealloc
#' @export
#' @importFrom raster ncell mask crop stack cellStats raster terrain which.min
#' @importFrom terra writeVector vect
#' @importFrom elevatr get_elev_raster
#' @importFrom grDevices terrain.colors topo.colors grey
#' @importFrom graphics legend
#'
#' @examples
#' # load a sample Digital Terrain Model
#' data(volc)
#'
#'
#' # load the sample locations on the above DTM
#' data(destin.loc)
#'
#'
#' #carry out a cost allocation analysis using the Tobler's off-path hiking function,
#' #setting the time to minutes and the isolines (i.e., contours) interval to 1 minute;
#' #only 3 locations are used
#'
#' #result <- movealloc(dtm=volc, origin=destin.loc[c(3,7,9),], funct="tofp", time="m",
#' #breaks=1, isolines=TRUE)
#'
#'
#' #same as above, using all the locations and the Kondo-Seino's cost function
#'
#' #result <- movealloc(dtm=volc, origin=destin.loc, funct="ks", time="m",
#' #breaks=1, isolines=TRUE)
#'
#'
#' @seealso \code{\link{movecost}}
#'
#'
movealloc <- function (dtm=NULL, origin, studyplot=NULL, funct="t", time="h", move=16, cogn.slp=FALSE, sl.crit=10, W=70, L=0, N=1, V=1.2, z=9, isolines=FALSE, breaks=NULL, cont.lab=TRUE, cex.breaks=0.6, leg.alloc=FALSE, leg.pos="topright", cex.leg=0.75, transp=0.5, export=FALSE){

  #if no dtm is provided
  if (is.null(dtm)==TRUE) {
    #get the elvation data using the elevatr's get_elev_raster() function, using the studyplot dataset (SpatialPolygonDataFrame)
    #to select the area whose elevation data are to be downloaded;
    #z sets the resolution of the elevation datataset
    elev.data <- elevatr::get_elev_raster(studyplot, z = z, verbose=FALSE, override_size_check = TRUE)
    #crop the elevation dataset to the exact boundary of the studyplot dataset
    dtm <- raster::crop(elev.data, studyplot)
  }

  #create an empty list to store the accumulated cost rasters around the origins,
  #which will be eventually stacked
  rasters_to_stack <- list()

  for (i in 1:length(origin)) {
    rasters_to_stack[[i]] <- movecost::movecost(dtm=dtm, origin=origin[i,], studyplot=studyplot, funct=funct, time=time, move=move, cogn.slp=cogn.slp, sl.crit=sl.crit, W=W, L=L, N=N, V=V, z=z, return.base=FALSE, graph.out=FALSE)$accumulated.cost.raster
  }

  #stack the accumulated cost rasters created in the previous loop
  res_stack <- raster::stack(rasters_to_stack)

  #looks at each pixel and returns 1 if the first stack has the smallest pixel value,
  #or 2 if the second stack has the smallest pixel value (and so on for bigger stacks)
  cost_alloc <- raster::which.min(res_stack)

  #define different types of cost functions and set the appropriate text to be used for subsequent plotting
  if (funct=="t") {
    main.title.accum.cost <- paste0("Accumulated walking-time cost (in ", time, ") around origins")
    main.title.cost.alloc <- paste0("Cost allocation based on walking-time (in ", time, ") around origins")
    sub.title <- paste0("Walking-time based on the Tobler's on-path hiking function \nterrain factor N=", N)
    legend.cost <- paste0("walking-time (", time,")")
  }

  if (funct=="tofp") {
    main.title.accum.cost <- paste0("Accumulated walking-time cost (in ", time, ") around origins")
    main.title.cost.alloc <- paste0("Cost allocation based on walking-time (in ", time, ") around origins")
    sub.title <- "Walking-time based on the Tobler's off-path hiking function"
    legend.cost <- paste0("walking-time (", time,")")
  }

  if(funct=="mp") {
    main.title.accum.cost  <- paste0("Accumulated walking-time cost (in ", time, ") around origins")
    main.title.cost.alloc <- paste0("Cost allocation based on walking-time (in ", time, ") around origins")
    sub.title <- paste0("Walking-time based on the Marquez-Perez et al.'s modified Tobler hiking function \n terrain factor N=", N)
    legend.cost <- paste0("walking-time (", time,")")
  }

  if(funct=="alb") {
    main.title.accum.cost <- paste0("Accumulated walking-time cost (in ", time, ") around origins")
    main.title.cost.alloc <- paste0("Cost allocation based on walking-time (in ", time, ") around origins")
    sub.title <- "Walking-time based on the Alberti's modified Tobler hiking function"
    legend.cost <- paste0("walking-time (", time,")")
  }

  if(funct=="icmonp") {
    main.title.accum.cost <- paste0("Accumulated walking-time cost (in ", time, ") around origins")
    main.title.cost.alloc <- paste0("Cost allocation based on walking-time (in ", time, ") around origins")
    sub.title <- paste0("Walking-time based on the (male, on-path) Irmischer-Clarke's hiking function \n terrain factor N=", N)
    legend.cost <- paste0("walking-time (", time,")")
  }

  if(funct=="icmoffp") {
    main.title.accum.cost <- paste0("Accumulated walking-time cost (in ", time, ") around origins")
    main.title.cost.alloc <- paste0("Cost allocation based on walking-time (in ", time, ") around origins")
    sub.title <- "Walking-time based on the (male, off-path) Irmischer-Clarke's hiking function"
    legend.cost <- paste0("walking-time (", time,")")
  }

  if(funct=="icfonp") {
    main.title.accum.cost <- paste0("Accumulated walking-time cost (in ", time, ") around origins")
    main.title.cost.alloc <- paste0("Cost allocation based on walking-time (in ", time, ") around origins")
    sub.title <- paste0("Walking-time based on the (female, on-path) Irmischer-Clarke's hiking function \n terrain factor N=", N)
    legend.cost <- paste0("walking-time (", time,")")
  }

  if(funct=="icfoffp") {
    main.title.accum.cost <- paste0("Accumulated walking-time cost (in ", time, ") around origins")
    main.title.cost.alloc <- paste0("Cost allocation based on walking-time (in ", time, ") around origins")
    sub.title <- "Walking-time based on the (female, off-path) Irmischer-Clarke's hiking function"
    legend.cost <- paste0("walking-time (", time,")")
  }

  if(funct=="gkrs") {
    main.title.accum.cost <- paste0("Accumulated walking-time cost (in ", time, ") around origin")
    main.title.cost.alloc <- paste0("Cost allocation based on walking-time (in ", time, ") around origins")
    sub.title <- paste0("Walking-time based on the Garmy et al.'s hiking function \n terrain factor N=", N)
    legend.cost <- paste0("walking-time (", time,")")
  }

  if(funct=="ug") {
    main.title.accum.cost <- paste0("Accumulated walking-time cost (in ", time, ") around origins")
    main.title.cost.alloc <- paste0("Cost allocation based on walking-time (in ", time, ") around origins")
    sub.title <- paste0("Walking-time based on the Uriarte Gonzalez's hiking function \n terrain factor N=", N)
    legend.cost <- paste0("walking-time (", time,")")
  }

  if(funct=="r") {
    main.title.accum.cost <- paste0("Accumulated walking-time cost (in ", time, ") around origins")
    main.title.cost.alloc <- paste0("Cost allocation based on walking-time (in ", time, ") around origins")
    sub.title <- paste0("Walking-time based on the Rees' hiking function \n terrain factor N=", N)
    legend.cost <- paste0("walking-time (", time,")")
  }

  if(funct=="ks") {
    main.title.accum.cost <- paste0("Accumulated walking-time cost (in ", time, ") around origins")
    main.title.cost.alloc <- paste0("Cost allocation based on walking-time (in ", time, ") around origins")
    sub.title <- paste0("Walking-time based on the Kondo-Seino's hiking function \n terrain factor N=", N)
    legend.cost <- paste0("walking-time (", time,")")
  }

  if(funct=="ree") {
    main.title.accum.cost <- "Accumulated energetic cost around origins"
    main.title.cost.alloc <- "Allocation based on energetic cost accumulated around origins"
    sub.title <- paste0("Cost based on the relative energetic expenditure cost function \n terrain factor N=", N)
    legend.cost <- "cost"
  }

  if(funct=="hrz") {
    main.title.accum.cost <- "Accumulated metabolic cost around origins"
    main.title.cost.alloc <- "Allocation based on metabolic cost accumulated around origins"
    sub.title <- paste0("Cost based on the Herzog's metabolic cost function \n cost in J / (Kg*m) \n terrain factor N=", N)
    legend.cost <- "metabolic cost J / (Kg*m)"
  }

  if(funct=="wcs") {
    main.title.accum.cost <- "Accumulated cost around origins"
    main.title.cost.alloc <- "Allocation based on cost accumulated around origins"
    sub.title <- paste0("Cost based on the wheeled-vehicle critical slope cost function \ncritical slope set to ", sl.crit, " percent \n terrain factor N=", N)
    legend.cost <- "cost"
  }

  if(funct=="vl") {
    main.title.accum.cost <- "Accumulated metabolic cost around origins"
    main.title.cost.alloc <- "Allocation based on metabolic cost accumulated around origins"
    legend.cost <- "energy expenditure cost (Megawatts)"
    if (V==0) {
      sub.title <- paste0("Cost based on the Van Leusen's metabolic energy expenditure cost function \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V is based on the Tobler on-path hiking function")
    } else {
      sub.title <- paste0("Cost based on the Van Leusen's metabolic energy expenditure cost function \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V: ", V)
    }
  }

  if(funct=="p") {
    main.title.accum.cost <- "Accumulated metabolic cost around origins"
    main.title.cost.alloc <- "Allocation based on metabolic cost accumulated around origins"
    legend.cost <- "energy expenditure cost (Megawatts)"
    if (V==0) {
      sub.title <- paste0("Cost based on the Pandolf et al.'s metabolic energy expenditure cost function \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V is based on the Tobler on-path hiking function")
    } else {
      sub.title <- paste0("Cost based on the Pandolf et al.'s metabolic energy expenditure cost function \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V: ", V)
    }
  }

  if(funct=="pcf") {
    main.title.accum.cost <- "Accumulated metabolic cost around origins"
    main.title.cost.alloc <- "Allocation based on metabolic cost accumulated around origins"
    legend.cost <- "energy expenditure cost (Megawatts)"
    if (V==0) {
      sub.title <- paste0("Cost based on the Pandolf et al.'s metabolic energy expenditure cost function with correction factor \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V is based on the Tobler on-path hiking function")
    } else {
      sub.title <- paste0("Cost based on the Pandolf et al.'s metabolic energy expenditure cost function with correction factor \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V: ", V)
    }
  }

  if(funct=="ls") {
    main.title.accum.cost <- "Accumulated metabolic cost around origins"
    main.title.cost.alloc <- "Allocation based on metabolic cost accumulated around origins"
    sub.title <- paste0("Cost based on the Llobera-Sluckin's metabolic energy expenditure cost function \n terrain factor N=", N)
    legend.cost <- "energy expenditure cost (KJ/m)"
  }

  if(funct=="b") {
    main.title.accum.cost <- "Accumulated cost around origins"
    main.title.cost.alloc <- "Allocation based on cost accumulated around origins"
    sub.title <- paste0("Cost based on the Bellavia's cost function \n terrain factor N=", N)
    legend.cost <- "cost"
  }

  if(funct=="m") {
    main.title.accum.cost <- "Accumulated metabolic cost around origins"
    main.title.cost.alloc <- "Allocation based on metabolic cost accumulated around origins"
    sub.title <- paste0("Cost based on the Minetti et al.'s  metabolic cost function \n cost in J / (Kg*m) \n terrain factor N=", N)
    legend.cost <- "metabolic cost J / (Kg*m)"
  }

  if (funct=="ma") {
    main.title.accum.cost <- paste0("Accumulated walking-time cost (in ", time, ") around origins")
    main.title.cost.alloc <- paste0("Cost allocation based on walking-time (in ", time, ") around origins")
    sub.title <- paste0("Walking-time based on the Marin Arroyo's hiking function \nterrain factor N=", N)
    legend.cost <- paste0("walking-time (", time,")")
  }

  if(funct=="a") {
    main.title.accum.cost <- "Accumulated metabolic cost around origins"
    main.title.cost.alloc <- "Allocation based on metabolic cost accumulated around origins"
    legend.cost <- "energy expenditure cost J / (Kg*m)"
    if (V==0) {
      sub.title <- paste0("Cost based on the Ardigo et al.'s metabolic energy expenditure cost function \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V is based on the Tobler on-path hiking function")
    } else {
      sub.title <- paste0("Cost based on the Ardigo et al.'s metabolic energy expenditure cost function \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V: ", V)
    }
  }

  if(funct=="h") {
    main.title.accum.cost <- "Accumulated metabolic cost around origins"
    main.title.cost.alloc <- "Allocation based on metabolic cost accumulated around origins"
    sub.title <- paste0("Cost based on the Hare's  metabolic cost function \n cost in cal/km \n terrain factor N=", N)
    legend.cost <- "metabolic cost cal/km"
  }

  if(funct=="e") {
    main.title.accum.cost <- "Accumulated cost around origins"
    main.title.cost.alloc <- "Allocation based on cost accumulated around origins"
    sub.title <- paste0("Cost based on the Eastman's cost function \n terrain factor N=", N)
    legend.cost <- "cost"
  }

  if (funct=="trp") {
    main.title.accum.cost <- paste0("Accumulated walking-time cost (in ", time, ") around origins")
    main.title.cost.alloc <- paste0("Cost allocation based on walking-time (in ", time, ") around origins")
    sub.title <- paste0("Walking-time based on the Tripcevich's hiking function \nterrain factor N=", N)
    legend.cost <- paste0("walking-time (", time,")")
  }

  #turn the allocation raster into polygons, so to get
  #boundaries to be used for subsequent plotting;
  #the cells storing the same integer are dissolved
  alloc.boundaries <- raster::rasterToPolygons(cost_alloc, dissolve=TRUE)

  #create the ingredient for the slopeshade raster to be plotted
  slope <- raster::terrain(dtm, opt = "slope")

  #plot the cost allocation raster
  raster::plot(cost_alloc,
               col=topo.colors(length(origin)),
               legend=FALSE,
               main=main.title.cost.alloc,
               cex.main=0.95,
               sub=sub.title,
               cex.sub=0.75)

  #if leg is TRUE, plot the legend
  if (leg.alloc==TRUE) {
    graphics::legend(x=leg.pos,
           legend = as.factor(1:length(origin)),
           fill = topo.colors(length(origin)),
           cex=cex.leg)
  }

  #plot the slopeshade raster
  raster::plot(slope,
               col = rev(grey(0:100/100)),
               legend = FALSE,
               alpha=transp,
               add=TRUE)

  #plot the allocation boundaries (i.e., polygons)
  raster::plot(alloc.boundaries, add=TRUE)

  #plot the origins
  raster::plot(origin,
               add=TRUE,
               pch=20)

  #if isolines is TRUE
  if(isolines==TRUE){
    #calculate the accumulated cost surface considering all the origin locations
    acc_cost_all_orig <- movecost::movecost(dtm=dtm, origin=origin, studyplot=studyplot, funct=funct, time=time, move=move, cogn.slp=cogn.slp, sl.crit=sl.crit, W=W, L=L, N=N, V=V, z=z, return.base=FALSE, graph.out=FALSE)

    #if no break value is entered, set the breaks to one tenth of the range of the values of the final accumulated cost surface
    if(is.null(breaks)==TRUE){
      #exclude the inf values from the calculation
      breaks <- round((max(acc_cost_all_orig$accumulated.cost.raster[][is.finite(acc_cost_all_orig$accumulated.cost.raster[])]) - min(acc_cost_all_orig$accumulated.cost.raster[][is.finite(acc_cost_all_orig$accumulated.cost.raster[])])) / 10,2)
    }

    #set the break values for the isolines, again excluding inf values
    levels <- seq(min(acc_cost_all_orig$accumulated.cost.raster[][is.finite(acc_cost_all_orig$accumulated.cost.raster[])]), max(acc_cost_all_orig$accumulated.cost.raster[][is.finite(acc_cost_all_orig$accumulated.cost.raster[])]), breaks)

    #add the contours (i.e., isolines), using the interval defined at the preceding two steps
    raster::contour(acc_cost_all_orig$accumulated.cost.raster,
                    add=TRUE,
                    col="white",
                    levels=levels,
                    labcex=cex.breaks,
                    drawlabels = cont.lab)
  }

  #if no DTM was provided (i.e., if 'studyplot' is not NULL), export the downloaded DTM as a raster file
  if(export==TRUE & is.null(studyplot)==FALSE){
    raster::writeRaster(dtm, "dtm", format="GTiff")
  }

  #if export is TRUE, export the cost allocation raster and the allocation polygons
  if(export==TRUE){
    raster::writeRaster(cost_alloc, paste0("cost_alloc_raster_", funct), format="GTiff")
    terra::writeVector(vect(alloc.boundaries), filename=paste0("alloc_boundaries_", funct), filetype="ESRI Shapefile")
  }

  #if export is TRUE and the isolines is TRUE, export the isolines
  if(export==TRUE & isolines==TRUE){
    terra::writeVector(vect(acc_cost_all_orig$isolines), filename=paste0("isolines_", funct), filetype="ESRI Shapefile")
  }

  #if isolines is FALSE
  if(isolines==FALSE){
    #create a list to be returned not containing the isolines
    rslt <- list("dtm"=dtm,
                 "cost.allocation.raster"=cost_alloc,
                 "alloc.boundaries" = alloc.boundaries)

  } else {
    #otherwise create a list containing the isolines
    rslt <- list("dtm"=dtm,
                 "cost.allocation.raster"=cost_alloc,
                 "isolines"=acc_cost_all_orig$isolines,
                 "alloc.boundaries" = alloc.boundaries)
  }

  result <- rslt

  return(result)
}

