#' R function for calculating least-cost path network
#'
#' The function provides the facility to calculate LCPs between multiple origins. The result is a network of LCPs connecting
#' each origin location to all the others locations. Optionally, a raster representing the density of the LCPs network can be produced.
#'
#' Like 'movecost()', the function just requires an input DTM ('RasterLayer' class) and an origin dataset  ('SpatialPointsDataFrame' class).
#' If a DTM is not provided, 'movenetw()' downloads elevation data from online sources for the area enclosed by the polygon fed via
#' the 'studyplot' parameter (see \code{\link{movecost}} for more details). Under the hood, movenetw()' relies on 'movecost()' and implements the same cost functions:
#' see the help documentation of 'movecost()' for further information.\cr
#'
#' 'movenetw()' produces a plot representing the input DTM overlaid by an hillshade raster, whose transparency can be adjusted using
#' the 'transp' parameter. On the rendered plot, the LPCs network ('SpatialLinesDataFrame' class) is represented by black lines. Optionally,
#' by setting the parameter 'lcp.dens' to TRUE, the function produces a raster representing the density of LCPs. The raster, rendered
#' overlaid to an hillshade visualization, expresses the density of LCPs as percentages. The percentages are calculated in relation to the
#' maximum number of LCPs passing through the same cell stored in the raster. A density raster expressing counts is NOT rendered BUT is returned by the function.
#' The density raster retains the cell size and coordinate system of the input DTM.\cr
#'
#' The function returns a list storing the DTM (only in case this has not been fed into the function but aquired online), a list of LCPs
#' split by origin, a SpatialLineDataFrame representing the merged LCPs, and two rasters representing the LCPs network density
#' expressed as counts and percentages respectively.
#'
#' The mentioned data can be exported by setting the 'export' parameter to TRUE. The LCPs network and the density raster will
#' bear a suffix indicating the used cost function.\cr
#'
#' @param dtm Digital Terrain Model (RasterLayer class); if not provided, elevation data will be acquired online for the area enclosed by the 'studyplot' parameter (see \code{\link{movecost}}).
#' @param origin location(s) around which the boundary(ies) is calculated (SpatialPointsDataFrame class).
#' @param studyplot polygon (SpatialPolygonDataFrame class) representing the study area for which online elevation data are aquired (see \code{\link{movecost}}); NULL is default.
#' @param funct cost function to be used (for details on each of the following, see \code{\link{movecost}}):\cr
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
#'
#' \strong{-functions expressing cost as metabolic energy expenditure-}\cr
#' \strong{p} uses the Pandolf et al.'s metabolic energy expenditure cost function;\cr
#' \strong{pcf} uses the Pandolf et al.'s cost function with correction factor for downhill movements;\cr
#' \strong{m} uses the Minetti et al.'s metabolic energy expenditure cost function;\cr
#' \strong{hrz} uses the Herzog's metabolic energy expenditure cost function;\cr
#' \strong{vl} uses the Van Leusen's metabolic energy expenditure cost function;\cr
#' \strong{ls} uses the Llobera-Sluckin's metabolic energy expenditure cost function;\cr
#' \strong{a} uses the Ardigo et al.'s metabolic energy expenditure cost function (for all the mentioned cost functions, see Details);\cr
#' @param move number of directions in which cells are connected: 4 (rook's case), 8 (queen's case), 16 (knight and one-cell queen moves; default).
#' @param cogn.slp  TRUE or FALSE (default) if the user wants or does not want the 'cognitive slope' to be used in place of the real slope (see \code{\link{movecost}}).
#' @param sl.crit critical slope (in percent), typically in the range 8-16 (10 by default) (used by the wheeled-vehicle cost function; see \code{\link{movecost}}).
#' @param W walker's body weight (in Kg; 70 by default; used by the Pandolf's and Van Leusen's cost function; see \code{\link{movecost}}).
#' @param L carried load weight (in Kg; 0 by default; used by the Pandolf's and Van Leusen's cost function; see \code{\link{movecost}}).
#' @param N coefficient representing ease of movement (1 by default) (see \code{\link{movecost}}).
#' @param V speed in m/s (1.2 by default) (used by the Pandolf et al.'s, Pandolf et al.s with correction factor, Van Leusen's, and Ardigo et al.'s cost function; if set to 0, it is internally worked out on the basis of Tobler on-path hiking function (see \code{\link{movecost}}).
#' @param z zoom level for the elevation data downloaded from online sources (from 0 to 15; 9 by default) (see \code{\link{movecost}} and \code{\link[elevatr]{get_elev_raster}}).
#' @param lcp.dens TRUE or FALSE (default) if the user wants or does not want the least-cost paths density raster to be produced.
#' @param transp set the transparency of the hillshade raster that is plotted over the DTM (0.5 by default).
#' @param oneplot TRUE (default) or FALSE if the user wants or does not want the plots displayed in a single window.
#' @param export TRUE or FALSE (default) if the user wants or does not want the LCPs network to be exported as a shapefile, and the LCPs network density as a raster; the DTM is exported only if it was not provided by the user
#' and downloaded by the function from online sources.
#'
#' @return The function returns a list storing the following components \itemize{
##'  \item{dtm: }{Digital Terrain Model ('RasterLayer' class); returned only if aquired online}
##'  \item{LCPs.netw: }{list containing the LCPs ('SpatialLinesDataFrame' class) split by origin}
##'  \item{LCPs.netw.merged: }{'SpatialLinesDataFrame' corresponding to the merged LCPs}
##'  \item{LCPs.density.count: }{raster ('RasterLayer' class) representing the counts of LCPs on each raster's cell}
##'  \item{LCPs.density.perc: }{same as the preceding, but re-expressing the counts as percentages)}
##' }
##'
#' @keywords movenetw
#' @export
#' @importFrom raster raster hillShade terrain rasterize cellStats extend crs
#' @importFrom spatstat.geom pixellate
#' @import maptools
#' @importFrom elevatr get_elev_raster
#' @importFrom grDevices terrain.colors topo.colors grey
#' @importFrom methods as
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
#' # calculate the least-cost path network using the Tobler's hiking
#' # function (for on-path walking)
#'
#' result <- movenetw(dtm=volc, origin=destin.loc[c(1,2,4),], move=8, funct="t")
#'
#'
#' @seealso \code{\link{movecost}}
#'
#'
movenetw <- function (dtm=NULL, origin, studyplot=NULL, funct="t", move=16, cogn.slp=FALSE, sl.crit=10, W=70, L=0, N=1, V=1.2, z=9, lcp.dens=FALSE, transp=0.5, oneplot=TRUE, export=FALSE){
  #if no dtm is provided
  if (is.null(dtm)==TRUE) {
    #get the elvation data using the elevatr's get_elev_raster() function, using the studyplot dataset (SpatialPolygonDataFrame)
    #to select the area whose elevation data are to be downloaded;
    #z sets the resolution of the elevation datataset
    elev.data <- elevatr::get_elev_raster(studyplot, z = z, verbose=FALSE, override_size_check = TRUE)
    #crop the elevation dataset to the exact boundary of the studyplot dataset
    dtm <- raster::crop(elev.data, studyplot)
  }

  #create an empty list where to store the computed LCPs
  paths.netw <- vector(mode = "list", length = length(origin))

  message(paste0("Wait...processing ", length(origin)*(length(origin)-1), " LCPs..."))

  #set the progress bar
  pb <- txtProgressBar(min = 0, max = length(origin), style = 3)

  #loop through the origins and store the LCPs (per each origin) in the previously defined list
  for (i in 1:length(origin)) {
    paths.netw[[i]] <- movecost(dtm = dtm, origin = origin[i,], destin = origin[-i,], funct = funct, move = move, sl.crit = sl.crit, W = W, L = L, N = N, V = V, graph.out = FALSE)$LCPs
    setTxtProgressBar(pb, i)
  }

  #merge the LCPs
  paths.netw.merged<- do.call(rbind, paths.netw)

  #define different types of cost functions and set the appropriate text to be used for subsequent plotting
  if (funct=="t") {
    sub.title <- paste0("LCPs based on the Tobler's on-path hiking function \nterrain factor N=", N)
  }

  if (funct=="tofp") {
    sub.title <- "LCPs based on the Tobler's off-path hiking function"
  }

  if(funct=="mp") {
    sub.title <- paste0("LCPs based on the Marquez-Perez et al.'s modified Tobler hiking function \n terrain factor N=", N)
  }

  if(funct=="alb") {
    sub.title <- "LCPs based on the Alberti's modified Tobler hiking function"
  }

  if(funct=="icmonp") {
    sub.title <- paste0("LCPs based on the (male, on-path) Irmischer-Clarke's hiking function \n terrain factor N=", N)
  }

  if(funct=="icmoffp") {
    sub.title <- "LCPs based on the (male, off-path) Irmischer-Clarke's hiking function"
  }

  if(funct=="icfonp") {
    sub.title <- paste0("LCPs based on the (female, on-path) Irmischer-Clarke's hiking function \n terrain factor N=", N)
  }

  if(funct=="icfoffp") {
    sub.title <- "LCPs based on the (female, off-path) Irmischer-Clarke's hiking function"
  }

  if(funct=="gkrs") {
    sub.title <- paste0("LCPs based on the Garmy et al.'s hiking function \n terrain factor N=", N)
  }

  if(funct=="ug") {
    sub.title <- paste0("LCPs based on the Uriarte Gonzalez's hiking function \n terrain factor N=", N)
  }

  if(funct=="r") {
    sub.title <- paste0("LCPs based on the Rees' hiking function \n terrain factor N=", N)
  }

  if(funct=="ks") {
    sub.title <- paste0("LCPs based on the Kondo-Seino's hiking function \n terrain factor N=", N)
  }

  if(funct=="ree") {
    sub.title <- paste0("LCPs based on the relative energetic expenditure cost function \n terrain factor N=", N)
  }

  if(funct=="hrz") {
    sub.title <- paste0("LCPs based on the Herzog's metabolic cost function \n cost in J / (Kg*m) \n terrain factor N=", N)
  }

  if(funct=="wcs") {
    sub.title <- paste0("LCPs based on the wheeled-vehicle critical slope cost function \ncritical slope set to ", sl.crit, " percent \n terrain factor N=", N)
  }

  if(funct=="vl") {
    if (V==0) {
      sub.title <- paste0("LCPs based on the Van Leusen's metabolic energy expenditure cost function \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V is based on the Tobler on-path hiking function")
    } else {
      sub.title <- paste0("LCPs based on the Van Leusen's metabolic energy expenditure cost function \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V: ", V)
    }
  }

  if(funct=="p") {
    if (V==0) {
      sub.title <- paste0("LCPs based on the Pandolf et al.'s metabolic energy expenditure cost function \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V is based on the Tobler on-path hiking function")
    } else {
      sub.title <- paste0("LCPs based on the Pandolf et al.'s metabolic energy expenditure cost function \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V: ", V)
    }
  }

  if(funct=="pcf") {
    if (V==0) {
      sub.title <- paste0("LCPs based on the Pandolf et al.'s metabolic energy expenditure cost function with correction factor \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V is based on the Tobler on-path hiking function")
    } else {
      sub.title <- paste0("LCPs based on the Pandolf et al.'s metabolic energy expenditure cost function with correction factor \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V: ", V)
    }
  }

  if(funct=="ls") {
    sub.title <- paste0("LCPs based on the Llobera-Sluckin's metabolic energy expenditure cost function \n terrain factor N=", N)
  }

  if(funct=="b") {
    sub.title <- paste0("LCPs based on the Bellavia's cost function \n terrain factor N=", N)
  }

  if(funct=="m") {
    sub.title <- paste0("LCPs based on the Minetti et al.'s  metabolic cost function \n cost in J / (Kg*m) \n terrain factor N=", N)
  }

  if (funct=="ma") {
    sub.title <- paste0("LCPs based on the Marin Arroyo's hiking function \nterrain factor N=", N)
  }

  if(funct=="a") {
    if (V==0) {
      sub.title <- paste0("LCPs based on the Ardigo et al.'s metabolic energy expenditure cost function \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V is based on the Tobler on-path hiking function")
    } else {
      sub.title <- paste0("LCPs based on the Ardigo et al.'s metabolic energy expenditure cost function \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V: ", V)
    }
  }

  #produce the ingredients for the hillshade raster
  slope <- raster::terrain(dtm, opt = "slope")
  aspect <- raster::terrain(dtm, opt = "aspect")
  hill <- raster::hillShade(slope, aspect, angle = 45, direction = 0)

  #conditionally set the layout in just one visualization
  if(lcp.dens==TRUE & oneplot==TRUE){
    m <- rbind(c(1,2))
    layout(m)
  }

  #plot the DTM
  raster::plot(dtm,
               main="Network of least-cost paths among multiple origins",
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

  #add the path network
  raster:: plot(paths.netw.merged,
                add=TRUE)

  #add the origin(s)
  raster:: plot(origin,
                pch=20,
                col="red",
                add=TRUE)


  if(lcp.dens==TRUE){
    #convert the merged LCPs to a spatstat object
    #using the as() function out of maptool package
    merged.lines.obj <- as(paths.netw.merged, "psp")

    #extract the resolution of the DTM, to be used for the eps parameter in 'pixellate'
    dtm.resolution <- raster::res(dtm)[1]

    #create a pixel image that counts the LCPs
    #the spatstat's parameter eps define the window determining the pixel resolution;
    #this is used so that the count raster will retain the dtm resolution
    lcps_count <- spatstat.geom::pixellate(merged.lines.obj, eps= dtm.resolution, what="number")

    #convert the spatstat's pixel image to a raster
    pathcounts <- raster::raster(lcps_count)

    #convert the path counts to percentage
    lcp.density.perc <- pathcounts / max(pathcounts[]) * 100

    #since the density rasters (produced by the spatstat's pixellate function)
    #has an extent equal to the extent of the LCPs (and not of the input DTM)
    #use the raster's 'extend' function to enlarge the extent of the density rasters
    #to the extent of the DTM, and set the new cells' value to 0
    pathcounts <-raster::extend(pathcounts, dtm, value=0)
    lcp.density.perc <-raster::extend(lcp.density.perc, dtm, value=0)

    #assign the crs of the dtm to the density rasters
    raster::crs(pathcounts) <- raster::crs(dtm)
    raster::crs(lcp.density.perc) <- raster::crs(dtm)

    #plot the DTM
    raster::plot(lcp.density.perc,
                 main="Density of least-cost paths among multiple origins",
                 sub=sub.title,
                 cex.main=0.95,
                 cex.sub=0.75,
                 legend.lab="% of LCPs passing in each raster's cell")

    #plot the hillshade
    raster::plot(hill,
                 col = grey(0:100/100),
                 legend = FALSE,
                 alpha=transp,
                 add=TRUE)

    #add the origin(s)
    raster:: plot(origin,
                  pch=20,
                  col="red",
                  add=TRUE)
  }


  #if export is TRUE, export the LCPs network as a shapefile, and the density rasters as geotiff
  if(export==TRUE){
    rgdal::writeOGR(paths.netw.merged, ".", paste0("LCPs.netw.merged_", funct), driver="ESRI Shapefile")
    raster::writeRaster(pathcounts, paste0("LCPs_density_count_", funct, format="GTiff"))
    raster::writeRaster(lcp.density.perc, paste0("LCPs_density_perc_", funct, format="GTiff"))
  }

  #if no DTM was provided (i.e., if 'studyplot' is not NULL), export the downloaded DTM as a raster file
  if(export==TRUE & is.null(studyplot)==FALSE){
    raster::writeRaster(dtm, "dtm", format="GTiff")
  }

  #if studyplot is NULL (i.e., if the DTM has been provided)..
  if(is.null(studyplot)==TRUE){
    #set dtm to NULL, i.e. there is no DTM to return
    dtm <- NULL
  } else{
    #otherwise (i.e., if studyplot was not NULL), set the dtm to the downloaded dtm (i.e., there is
    #something to return)
    dtm <- dtm
  }

  #if the user set lcp.dens to TRUE
  if(lcp.dens==TRUE){
    #there are data to return
    pathcounts <- pathcounts
    lcp.density.perc <- lcp.density.perc
  } else{
    #otherwise there is nothing to return
    pathcounts <- NULL
    lcp.density.perc <- NULL
  }

  #restore the original graphical device's settings if previously modified
  if(lcp.dens==TRUE & oneplot==TRUE){
    par(mfrow = c(1,1))
  }

  results <- list("dtm"=dtm,
                  "LCPs.netw"=paths.netw,
                  "LCPs.netw.merged"=paths.netw.merged,
                  "LCPs.density.count"=pathcounts,
                  "LCPs.density.perc"=lcp.density.perc)

  return(results)
}
