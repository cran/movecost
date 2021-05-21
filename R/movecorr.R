#' R function for calculating least-cost corridor between two point locations
#'
#' The function provides the facility to calculate the least-cost corridor between two point locations. It just requires an inout DTM,
#' and two point locations ('SpatialPointsDataFrame' class) representing the locations between which the corridor is calculated. Under the
#' hood, 'movecorr' relies on the \code{\link{movecost}} function and, needless to say, implements the same
#' cost functions. See the help documentation of \code{\link{movecost}} for further details.  The function renders a raster representing
#' the least cost corridor (with least-cost paths superimposed), and returns a list containing a number of components (see 'Value' below). The corridor raster can be
#' exported as GeoTiff (see 'Arguments' below). \cr
#'
#' What 'movecorr()' does is to calculate (via the \code{\link{movecost}} function) the accumulated cost surface around location
#' a, and then the accumulated cost surface around b. The two surfaces are eventually combined to produce the least-cost
#' corridor between location a and b. On the produced corridor raster, the cost of a cell is the total cost to reach it
#' from both locations. About least-cost corridors see for instance: \cr
#' Mitchell A. (2012), "The ESRI Guide to GIS Analysis. Vol 3. Modelling Suitability, Movement, and Interaction", New York: Esri Press (257-259).
#'
#'
#' @param dtm Digital Terrain Model (RasterLayer class); if not provided, elevation data will be acquired online for the area enclosed by the 'studyplot' parameter (see \code{\link{movecost}}).
#' @param a first location from which the least-cost corridor is calculated (SpatialPointsDataFrame class).
#' @param b second location from which the least-cost corridor is calculated (SpatialPointsDataFrame class).
#' @param lab.a string to be used to label point a on the outplut plot (A is the default)
#' @param lab.b string to be used to label point a on the outplut plot (B is the default).
#' @param cex.labs scaling factor for the size of the points' labels (0.8 by default)
#' @param studyplot polygon (SpatialPolygonDataFrame class) representing the study area for which online elevation data are aquired (see \code{\link{movecost}}); NULL is default.
#' @param funct cost function to be used (for details on each of the following, see \code{\link{movecost}}):\cr \strong{t} (default) uses the on-path Tobler's hiking function;\cr
#' \strong{tofp} uses the off-path Tobler's hiking function;\cr
#' \strong{mp} uses the Marquez-Perez et al.'s modified Tobler's function;\cr
#' \strong{icmonp} uses the Irmischer-Clarke's modified Tobler's hiking function (male, on-path);\cr
#' \strong{icmoffp} uses the Irmischer-Clarke's modified Tobler's hiking function (male, off-path);\cr
#' \strong{icfonp} uses the Irmischer-Clarke's modified Tobler's hiking function (female, on-path);\cr
#' \strong{icfoffp} uses the Irmischer-Clarke's modified Tobler's hiking function (female, off-path);\cr
#' \strong{ug} uses the Uriarte Gonzalez's slope-dependant walking-time cost function;\cr
#' \strong{alb} uses the Alberti's Tobler hiking function modified for animal foraging excursions;\cr
#' \strong{ree} uses the relative energetic expenditure cost function;\cr
#' \strong{hrz} uses the Herzog's metabolic cost function;\cr
#' \strong{wcs} uses the wheeled-vehicle critical slope cost function;\cr
#' \strong{p} uses the Pandolf et al.'s metabolic energy expenditure cost function;\cr
#' \strong{vl} uses the Van Leusen's metabolic energy expenditure cost function;\cr
#' \strong{ls} uses the Llobera-Sluckin's metabolic energy expenditure cost function.
#' @param time time-unit expressed by the accumulated raster if Tobler's and other time-related cost functions are used;
#' 'h' for hour, 'm' for minutes.
#' @param move number of directions in which cells are connected: 4 (rook's case), 8 (queen's case), 16 (knight and one-cell queen moves; default).
#' @param sl.crit critical slope (in percent), typically in the range 8-16 (10 by default) (used by the wheeled-vehicle cost function; see \code{\link{movecost}}).
#' @param W walker's body weight (in Kg; 70 by default; used by the Pandolf's and Van Leusen's cost function; see \code{\link{movecost}}).
#' @param L carried load weight (in Kg; 0 by default; used by the Pandolf's and Van Leusen's cost function; see \code{\link{movecost}}).
#' @param N coefficient representing ease of movement (1 by default) (used by the Pandolf's and Van Leusen's cost function; see \code{\link{movecost}}).
#' @param V speed in m/s (1.2 by default) (used by the Pandolf's and Van Leusen's cost function; see \code{\link{movecost}}).
#' @param z zoom level for the elevation data downloaded from online sources (from 0 to 15; 9 by default) (see \code{\link{movecost}} and \code{\link[elevatr]{get_elev_raster}}).
#' @param export TRUE or FALSE (default) if the user wants or does not want the output to be exported; if TRUE, the least-cost corridor, the dtm (if not provided by the user but acquired online),
#' and the accumulated cost surface around a and b are exported as a GeoTiff file; the two LCPs (from a to b, and from b to a) as individual shapefiles.
#'
#' @return The function returns a list storing the following components \itemize{
##'  \item{dtm: }{Digital Terrain Model ('RasterLayer' class)}
##'  \item{lc.corridor: }{raster of the least-cost corridor ('RasterLayer' class)}
##'  \item{lcp_a_to_b: }{least-cost past from a to b ('SpatialLines' class)}
##'  \item{lcp_b_to_a: }{least-cost past from b to a ('SpatialLines' class)}
##'  \item{accum_cost_surf_a: }{accumulated cost-surface around a ('RasterLayer' class)}
##'  \item{accum_cost_surf_b: }{accumulated cost-surface around b ('RasterLayer' class)}
##' }
##'
#' @keywords movecorr
#' @export
#' @importFrom raster ncell mask crop
#' @importFrom elevatr get_elev_raster
#' @importFrom graphics layout par
#' @examples
#' # load a sample Digital Terrain Model
#' volc <- raster::raster(system.file("external/maungawhau.grd", package="gdistance"))
#'
#'
#' # load the sample destination locations on the above DTM
#' data(destin.loc)
#'
#' # calculate the least-cost corridor between two locations, using the
#' # relative energetic expenditure cost function, and store the results
#' # in the 'res' object
#'
#' result <- movecorr(dtm=volc, a=destin.loc[1,], b=destin.loc[3,], funct="ree")
#'
#'
#' @seealso \code{\link{movecost}}
#'
#'
movecorr <- function (dtm=NULL, a, b, lab.a="A", lab.b="B", cex.labs=0.8, studyplot=NULL, funct="t", time="h", move=16, sl.crit=10, W=70, L=0, N=1, V=1.2, z=9, export=FALSE){

  #calculate the accum cost surface and LCP around and from point a
  res.a <- movecost(dtm=dtm, origin=a, destin=b, studyplot=studyplot, funct=funct, time=time, move=move, sl.crit=sl.crit, W=W, L=L, N=N, V=V, z=z, graph.out=FALSE)

  #calculate the accum cost surface and LCP around and from point b;
  #uses the dtm stored in the output of the preceding calculation so that the
  #dtm does not have to be downloaded twice in case it was not originally provided by the user
  res.b <- movecost(dtm=res.a$dtm, origin=b, destin=a, studyplot=studyplot, funct=funct, time=time, move=move, sl.crit=sl.crit, W=W, L=L, N=N, V=V, z=z, graph.out=FALSE)

  #combine the 2 accumulated cost surfaces obtained at the preceding steps
  res.corridor <- res.a$accumulated.cost.raster + res.b$accumulated.cost.raster

  #define different types of cost functions and set the appropriate text to be used for subsequent plotting
  if (funct=="t") {

    #set the labels to be used within the returned plot
    main.title <- paste0("Least-cost corridor (cost in ", time, ")")
    sub.title <- "Walking-time based on the Tobler's on-path hiking function"
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
    sub.title <- "Walking-time based on the Marquez-Perez et al.'s modified Tobler hiking function"
    legend.cost <- paste0("walking-time (", time,")")
  }

  if(funct=="icmonp") {

    #set the labels to be used within the returned plot
    main.title <- paste0("Least-cost corridor (cost in ", time, ")")
    sub.title <- "Walking-time based on the (male, on-path) Irmischer-Clarke's modified Tobler hiking function"
    legend.cost <- paste0("walking-time (", time,")")
  }

  if(funct=="icmoffp") {

    #set the labels to be used within the returned plot
    main.title <- paste0("Least-cost corridor (cost in ", time, ")")
    sub.title <- "Walking-time based on the (male, off-path) Irmischer-Clarke's modified Tobler hiking function"
    legend.cost <- paste0("walking-time (", time,")")
  }

  if(funct=="icfonp") {

    #set the labels to be used within the returned plot
    main.title <- paste0("Least-cost corridor (cost in ", time, ")")
    sub.title <- "Walking-time based on the (female, on-path) Irmischer-Clarke's modified Tobler hiking function"
    legend.cost <- paste0("walking-time (", time,")")
  }

  if(funct=="icfoffp") {
    main.title <- paste0("Least-cost corridor (cost in ", time, ")")
    sub.title <- "Walking-time based on the (female, off-path) Irmischer-Clarke's modified Tobler hiking function"
    legend.cost <- paste0("walking-time (", time,")")
  }

  if(funct=="ug") {

    #set the labels to be used within the returned plot
    main.title <- paste0("Least-cost corridor (cost in ", time, ")")
    sub.title <- "Walking-time based on the Uriarte Gonzalez's cost function"
    legend.cost <- paste0("walking-time (", time,")")
  }

  if(funct=="ree") {

    #set the labels to be used within the returned plot
    main.title <- "Least-cost corridor"
    sub.title <- "Cost based on the slope-based relative energetic expenditure cost function"
    legend.cost <- "cost"
  }

  if(funct=="hrz") {

    #set the labels to be used within the returned plot
    main.title <- "Least-cost corridor"
    sub.title <- "Cost based on the Herzog's metabolic cost function \n cost in J / (Kg*m)"
    legend.cost <- "metabolic cost J / (Kg*m)"
  }

  if(funct=="wcs") {

    #set the labels to be used within the returned plot
    main.title <- "Least-cost corridor"
    sub.title <- paste0("Cost based on the wheeled-vehicle critical slope cost function \ncritical slope set to ", sl.crit, " percent")
    legend.cost <- "cost"
  }

  if(funct=="vl") {

    #set the labels to be used within the returned plot
    main.title <- "Least-cost corridor"
    sub.title <- paste0("Cost based on the Van Leusen's metabolic energy expenditure cost function \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V: ", V)
    legend.cost <- "energy expenditure cost (Megawatts)"
  }

  if(funct=="p") {

    #set the labels to be used within the returned plot
    main.title <- "Least-cost corridor"
    sub.title <- paste0("Cost based on the Pandolf et al.'s metabolic energy expenditure cost function \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V: ", V)
    legend.cost <- "energy expenditure cost (Megawatts)"
  }

  if(funct=="ls") {

    #set the labels to be used within the returned plot
    main.title <- "Least-cost corridor"
    sub.title <- paste0("Cost based on the Llobera-Sluckin's metabolic energy expenditure cost function")
    legend.cost <- "energy expenditure cost (KJ/m)"
  }

  if(funct=="alb") {

    #set the labels to be used within the returned plot
    main.title <- paste0("Least-cost corridor (cost in ", time, ")")
    sub.title <- "Walking-time based on the Alberti's modified Tobler hiking function"
    legend.cost <- paste0("walking-time (", time,")")
  }

  #plot the corridor raster
  raster::plot(res.corridor,
               main=main.title,
               sub=sub.title,
               cex.main=0.95,
               cex.sub=0.75,
               legend.lab=legend.cost)

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


  #if export is TRUE, export the LC corridor as a raster file
  #and the LCPs as shapefile
  if(export==TRUE){
    raster::writeRaster(res.corridor, paste0("LCcorridor_", funct), format="GTiff")
    raster::writeRaster(res.a$accumulated.cost.raster, paste0("Accum_cost_surf_a_", funct), format="GTiff")
    raster::writeRaster(res.b$accumulated.cost.raster, paste0("Accum_cost_surf_b_", funct), format="GTiff")
    rgdal::writeOGR(res.a$LCPs, ".", paste0("lcp_a_to_b_", funct), driver="ESRI Shapefile")
    rgdal::writeOGR(res.b$LCPs, ".", paste0("lcp_b_to_a_", funct), driver="ESRI Shapefile")
  }

  #if no DTM was provided (i.e., if 'studyplot' is not NULL), export the downloaded DTM as a raster file
  if(export==TRUE & is.null(studyplot)==FALSE){
    raster::writeRaster(res.a$dtm, "dtm", format="GTiff")
  }

  results <- list("dtm"=res.a$dtm,
                  "lc.corridor"=res.corridor,
                  "lcp_a_to_b"=res.a$LCPs,
                  "lcp_b_to_a"=res.b$LCPs,
                  "accum_cost_surf_a"=res.a$accumulated.cost.raster,
                  "accum_cost_surf_b"=res.b$accumulated.cost.raster)

  }
