#' R function for calculating least-cost path network
#'
#' The function provides the facility to calculate LCPs between multiple origins.
#' Two types of networks are produced: one where each origin location is connected to all the others locations; one where
#' only pairs of neighboring locations are connected. In other words, in the latter case, each location is connected to the location that is the
#' nearest in terms of walking cost, either in terms of time or energy (or abstract cost), according to the selected cost function.
#' Optionally, a raster representing the density of the first type of network can be produced.\cr
#' Visit this \href{https://drive.google.com/file/d/1gLDrkZFh1b_glzCEqKdkPrer72JJ9Ffa/view?usp=sharing}{LINK} to access the package's vignette.\cr
#'
#' Like \code{movecost()}, the function just requires an input DTM ('RasterLayer' class) and an origin dataset  ('SpatialPointsDataFrame' class).
#' If a DTM is not provided, \code{movenetw()} downloads elevation data from online sources for the area enclosed by the polygon fed via
#' the \code{studyplot} parameter (see \code{\link{movecost}} for more details). Under the hood, \code{movenetw()} relies on \code{movecost()} and implements the same cost functions:
#' see the help documentation of \code{movecost()} for further information.\cr
#'
#' \code{movenetw()} produces a plot representing the input DTM overlaid by a slopeshade raster, whose transparency can be adjusted using
#' the \code{transp} parameter. On the rendered plot, the LPCs network ('SpatialLinesDataFrame' class) is represented by black lines.
#' To calculate the network connecting all the locations, the user may want to set the \code{netw.type} parameter to \code{allpairs} (which is the default value).
#' If the user wants to calculate the network connecting neighbouring locations, the \code{netw.type} parameter is to be set to \code{neigh}.
#' Optionally, by setting the \code{lcp.dens} parameter to \code{TRUE}, the function produces a raster representing the density of the LCPs connecting each location to all the
#' other locations. The raster, which is rendered overlaid to a slopeshade visualization, expresses the density of LCPs as percentages.
#' The percentages are calculated in relation to the maximum number of LCPs passing through the same cell stored in the raster.
#' A density raster expressing counts is NOT rendered BUT is returned by the function. The density raster retains the cell size and coordinate system of the input DTM.\cr
#'
#' The function returns a list storing the DTM (only in case this has not been fed into the function but acquired online), a list of LCPs
#' split by origin, a SpatialLineDataFrame representing the merged LCPs, two rasters representing the LCPs network density
#' expressed as counts and percentages respectively, and cost matrices. As for the latter, if the selected cost function defines cost as
#' walking time, two matrices are returned, one expressing time in minutes, one in hours (\strong{note} that the values are in decimal format).
#' If the selected cost function expresses cost differently (i.e., energy or abstract cost), the two above mentioned cost matrices will be set to NULL,
#' and a third cost matrix will store all the pair-wise costs.\cr
#'
#' The above mentioned data (DTM, LCPs, network density) can be exported by setting the \code{export} parameter to \code{TRUE}. The LCPs network (exported as a shapefile)
#' and the density raster (as a GeoTiff) will bear a suffix indicating the used cost function.\cr
#'
#' @param dtm Digital Terrain Model (RasterLayer class); if not provided, elevation data will be acquired online for the area enclosed by the 'studyplot' parameter (see \code{\link{movecost}}).
#' @param origin locations from which the network of least-cost paths is calculated (SpatialPointsDataFrame class).
#' @param netw.type type of network to be calculated: 'allpairs' (default) or 'neigh' (see 'Details').
#' @param studyplot polygon (SpatialPolygonDataFrame class) representing the study area for which online elevation data are acquired (see \code{\link{movecost}}); NULL is default.
#' @param barrier area where the movement is inhibited (SpatialLineDataFrame or SpatialPolygonDataFrame class) (see \code{\link{movecost}}).
#' @param plot.barrier TRUE or FALSE (default) if the user wants or does not want the barrier to be plotted (see \code{\link{movecost}}).
#' @param irregular.dtm TRUE or FALSE (default) if the input DTM features irregular margins (see \code{\link{movecost}}).
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
#' @param move number of directions in which cells are connected: 4 (rook's case), 8 (queen's case), 16 (knight and one-cell queen moves; default).
#' @param field value assigned to the cells coinciding with the barrier (0 by default) (see \code{\link{movecost}}.
#' @param cogn.slp  TRUE or FALSE (default) if the user wants or does not want the 'cognitive slope' to be used in place of the real slope (see \code{\link{movecost}}).
#' @param sl.crit critical slope (in percent), typically in the range 8-16 (10 by default) (used by the wheeled-vehicle cost function; see \code{\link{movecost}}).
#' @param W walker's body weight (in Kg; 70 by default; used by the Pandolf's and Van Leusen's cost function; see \code{\link{movecost}}).
#' @param L carried load weight (in Kg; 0 by default; used by the Pandolf's and Van Leusen's cost function; see \code{\link{movecost}}).
#' @param N coefficient representing ease of movement (1 by default) (see \code{\link{movecost}}).
#' @param V speed in m/s (1.2 by default) (used by the Pandolf et al.'s, Pandolf et al.s with correction factor, Van Leusen's, and Ardigo et al.'s cost function; if set to 0, it is internally worked out on the basis of Tobler on-path hiking function (see \code{\link{movecost}}).
#' @param z zoom level for the elevation data downloaded from online sources (from 0 to 15; 9 by default) (see \code{\link{movecost}} and \code{\link[elevatr]{get_elev_raster}}).
#' @param lcp.dens TRUE or FALSE (default) if the user wants or does not want the least-cost paths density raster to be produced.
#' @param transp set the transparency of the slopeshade raster that is plotted over the DTM (0.5 by default).
#' @param export TRUE or FALSE (default) if the user wants or does not want the LCPs network to be exported as a shapefile, and the LCPs network density as a GeoTiff; the DTM is exported only if it was not provided by the user
#' and downloaded by the function from online sources.
#'
#' @return The function returns a list storing the following components \itemize{
##'  \item{dtm: }{Digital Terrain Model ('RasterLayer' class); returned only if acquired online}
##'  \item{LCPs.netw: }{list containing the LCPs ('SpatialLinesDataFrame' class) split by origin}
##'  \item{LCPs.netw.merged: }{'SpatialLinesDataFrame' corresponding to the merged LCPs}
##'  \item{LCPs.netw.neigh: }{list containing the LCPs between neighboring locations ('SpatialLinesDataFrame' class) split by origin}
##'  \item{LCPs.netw.merged: }{'SpatialLinesDataFrame' corresponding to the merged LCPs between neighboring locations}
##'  \item{LCPs.density.count: }{raster ('RasterLayer' class) representing the counts of LCPs on each raster's cell}
##'  \item{LCPs.density.perc: }{same as the preceding, but re-expressing the counts as percentages}
##'  \item{cost.matrix.min: }{matrix of cost between locations, expressing cost in minutes}
##'  \item{cost.matrix.hr: }{matrix of cost between locations, expressing cost in hours}
##'  \item{cost.matrix: }{matrix of cost between locations, expressing cost either in energy or abstract cost, depending on the
##'  used cost function}
##' }
##'
#' @keywords movenetw
#' @export
#' @importFrom raster raster terrain rasterize cellStats extend crs
#' @importFrom terra writeVector vect rast rasterizeGeom
#' @importFrom sp SpatialLinesDataFrame
#' @importFrom elevatr get_elev_raster
#' @importFrom grDevices terrain.colors topo.colors grey
#' @importFrom methods as
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
#' # calculate the least-cost path network between neighboring locations
#' # using the Tobler's hiking function (for on-path walking)
#'
#' result <- movenetw(dtm=volc, origin=destin.loc, move=8, funct="t", netw.type="neigh")
#'
#'
#' @seealso \code{\link{movecost}}
#'
#'
movenetw <- function (dtm=NULL, origin, netw.type="allpairs", studyplot=NULL, barrier=NULL, plot.barrier=FALSE, irregular.dtm=FALSE, funct="t", move=16, field=0, cogn.slp=FALSE, sl.crit=10, W=70, L=0, N=1, V=1.2, z=9, lcp.dens=FALSE, transp=0.5, export=FALSE){
  #if no dtm is provided
  if (is.null(dtm)==TRUE) {
    #get the elvation data using the elevatr's get_elev_raster() function, using the studyplot dataset (SpatialPolygonDataFrame)
    #to select the area whose elevation data are to be downloaded;
    #z sets the resolution of the elevation datataset
    elev.data <- elevatr::get_elev_raster(studyplot, z = z, verbose=FALSE, override_size_check = TRUE)
    #crop the elevation dataset to the exact boundary of the studyplot dataset
    dtm <- raster::crop(elev.data, studyplot)
  }

  #use the 'movecost' function to calculate the 'conductance' transitional layer to be used inside the
  #loop that follows. The conductance layer is calculated because it will be used later on
  # by the gdistance's shortestPath() and costDistance() functions.
  # Only one origin location suffices to calculate the conductance.
  message("Wait...calculating the conductance...")
  conductance.to.use <- movecost(dtm = dtm, origin = origin[1,], barrier=barrier, irregular.dtm=irregular.dtm, funct = funct, move = move, field=field, sl.crit = sl.crit, W = W, L = L, N = N, V = V, graph.out = FALSE)$conductance

  if (netw.type=="allpairs") {
    #create an empty list to store the computed LCPs
    paths.netw <- vector(mode = "list", length = length(origin))

    message(paste0("Wait...processing ", length(origin)*(length(origin)-1), " LCPs..."))

    #set the progress bar
    pb <- txtProgressBar(min = 0, max = length(origin), style = 3)

    #loop through the origins and store the LCPs (per each origin) in the previously defined list
    for (i in 1:length(origin)) {
      paths.netw[[i]] <- gdistance::shortestPath(conductance.to.use, sp::coordinates(origin[i,]), sp::coordinates(origin[-i,]), output="SpatialLines")
      setTxtProgressBar(pb, i)
    }
    #merge the LCPs
    paths.netw.merged<- do.call(rbind, paths.netw)
    #set to NULL variables that would be have been populated if netw.type was equal to 'neigh'
    paths.netw.neigh <- NULL
    paths.netw.neigh.merged <- NULL
  }

  #if a walking-time cost function has been selected...
  if (funct=="t" | funct=="tofp" | funct=="mp" | funct=="icmonp" | funct=="icmoffp" | funct=="icfonp" | funct=="icfoffp" | funct=="ug" | funct=="alb" | funct=="gkrs" | funct=="r" | funct=="ks" | funct=="ma"){
    #...create two matrices, one for minutes, one for hours (note: the original conductance trans layer is in seconds)
    cost.matrix.min <- gdistance::costDistance(conductance.to.use, sp::coordinates(origin)) / 60
    cost.matrix.hr <- gdistance::costDistance(conductance.to.use, sp::coordinates(origin)) / 3600
    cost.matrix <- NULL
  } else {
    #otherwise (i.e., if the cost function is not about time)
    #set the two preceding cost-matrices to NULL, and return another cost matrix
    cost.matrix.min <- NULL
    cost.matrix.hr <- NULL
    cost.matrix <- as.matrix(gdistance::costDistance(conductance.to.use, sp::coordinates(origin)))
  }

  if (netw.type=="neigh") {
    #if either the cost matrix in minutes or in hours is NULL (i.e., if energy functions have been used)...
    if(is.null(cost.matrix.min)==TRUE | is.null(cost.matrix.hr==TRUE)) {
      #...use the cost.matrix object
      #set the diagolnal to NA
      diag(cost.matrix) <- NA
      #extract the column index of the minimum cost; the search is carried out row-wisely
      min.index.col <- apply(cost.matrix, 1, which.min)
    } else {
      #if time functions have been used, do the same for the cost.matrix.min
      diag(cost.matrix.min) <- NA
      #extract the column index of the minimum cost; the search is carried out row-wisely
      min.index.col <- apply(cost.matrix.min, 1, which.min)
    }

    #create an empty list to store the computed LCPs among neighboring locations
    paths.netw.neigh <- vector(mode = "list", length = length(origin))

    message(paste0("Wait...processing the LCPs between neighboring locations..."))

    #set the progress bar
    pb <- txtProgressBar(min = 0, max = length(origin), style = 3)

    #loop through the origins and calculate the LCP between each origin
    #and its nearest neighbor, which is represented by the vector called min.index.col created previouly;
    #the LCPs (per each origin) are stored in the previously defined list
    for (i in 1:length(origin)) {
      paths.netw.neigh[[i]] <- gdistance::shortestPath(conductance.to.use, sp::coordinates(origin[i,]), sp::coordinates(origin[min.index.col[i],]), output="SpatialLines")
      setTxtProgressBar(pb, i)
    }
    #merge the LCPs between neigbhoring locations
    paths.netw.neigh.merged <- do.call(rbind, paths.netw.neigh)
    #set to NULL variables that would be have been populated if netw.type was equal to 'neigh'
    paths.netw <- NULL
    paths.netw.merged <- NULL
  }

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

  if(funct=="h") {
    sub.title <- paste0("LCPs based on the Hare's metabolic cost function \n cost in cal/km \n terrain factor N=", N)
  }

  if(funct=="e") {
    sub.title <- paste0("LCPs based on the Eastman's cost function \n terrain factor N=", N)
  }

  if (funct=="trp") {
    sub.title <- paste0("LCPs based on the Tripcevich's hiking function \nterrain factor N=", N)
  }

  #produce the ingredient for the slopeshade raster
  slope <- raster::terrain(dtm, opt = "slope")

  if (netw.type=="allpairs") {
    #plots for the LCPs among all locations
    #plot the DTM
    raster::plot(dtm,
                 main="Network of least-cost paths between multiple origins",
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

    #add the path network
    raster:: plot(paths.netw.merged,
                  add=TRUE)

    #add the origin(s)
    raster:: plot(origin,
                  pch=20,
                  col="red",
                  add=TRUE)

    #if the barrier is provided AND if plot.barrier is TRUE, add the barrier
    if(is.null(barrier)==FALSE & plot.barrier==TRUE) {
      raster::plot(barrier, col="blue", add=TRUE)
    }

    if(lcp.dens==TRUE){
      #convert the merged LCPs to a terra's SpatVector object
      merged.lines.obj <- terra::vect(paths.netw.merged)

      #convert the dtm to a terra's SpatRaster object
      dtm.spatrast <- terra::rast(dtm)

      #use the terra's function to count how many paths cross a cell
      lcps_count <- terra::rasterizeGeom(merged.lines.obj, dtm.spatrast, fun="crosses")

      #convert the path counts to percentage
      lcp.density.perc <- lcps_count / max(lcps_count[]) * 100

      #convert the lcps_count and lcp.density.perc from terra's SpatRaster to raster's Raster object class
      lcps_count <- raster::raster(lcps_count)
      lcp.density.perc <- raster::raster(lcp.density.perc)

      #plot the DTM
      raster::plot(lcp.density.perc,
                   main="Density of least-cost paths between multiple origins",
                   sub=sub.title,
                   cex.main=0.95,
                   cex.sub=0.75,
                   legend.lab="% of LCPs passing in each raster's cell")

      #plot the slopeshade
      raster::plot(slope,
                   col = rev(grey(0:100/100)),
                   legend = FALSE,
                   alpha=transp,
                   add=TRUE)

      #add the origin(s)
      raster:: plot(origin,
                    pch=20,
                    col="red",
                    add=TRUE)
    }
  }

  if (netw.type=="neigh") {
    #plots for the LCPs among neighboring locations
    #plot the DTM
    raster::plot(dtm,
                 main="Network of least-cost paths between neighboring locations",
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

    #add the path network
    raster:: plot(paths.netw.neigh.merged,
                  add=TRUE)

    #add the origin(s)
    raster:: plot(origin,
                  pch=20,
                  col="red",
                  add=TRUE)

    #if the barrier is provided AND if plot.barrier is TRUE, add the barrier
    if(is.null(barrier)==FALSE & plot.barrier==TRUE) {
      raster::plot(barrier, col="blue", add=TRUE)
    }
  }

  #if export is TRUE, export the LCPs network as a shapefile, and the density rasters as geotiff
  if(export==TRUE){
    #convert the merged LCPs network from SpatialLine to SpatialLineDataFrame via the 'SpatialLinesDataFrame()' function
    #this is needed by 'writeOGR()', and export relevant data
    if (netw.type=="allpairs") {
      paths.netw.merged.SPLDF <- sp::SpatialLinesDataFrame(paths.netw.merged, data.frame(id=1:length(paths.netw.merged)), match.ID = F)
      terra::writeVector(vect(paths.netw.merged.SPLDF), filename=paste0("LCPs.netw.merged_", funct), filetype="ESRI Shapefile")
    }
    if (netw.type=="neigh") {
      paths.netw.neigh.merged.SPLDF <- sp::SpatialLinesDataFrame(paths.netw.neigh.merged, data.frame(id=1:length(paths.netw.neigh.merged)), match.ID = F)
      terra::writeVector(vect(paths.netw.neigh.merged.SPLDF), filename=paste0("LCPs.netw.neigh.merged_", funct), filetype="ESRI Shapefile")
    }
    if(lcp.dens==TRUE & netw.type=="allpairs"){
      raster::writeRaster(lcps_count, paste0("LCPs_density_count_", funct), format="GTiff")
      raster::writeRaster(lcp.density.perc, paste0("LCPs_density_perc_", funct), format="GTiff")
    }
  }

  #if no DTM was provided (i.e., if 'studyplot' is not NULL), export the downloaded DTM as a raster file
  if(export==TRUE & is.null(studyplot)==FALSE){
    raster::writeRaster(dtm, "dtm", format="GTiff")
  }

  #if studyplot is NULL (i.e., if the DTM has been provided)..
  if(is.null(studyplot)==TRUE){
    #set dtm to NULL, i.e. there is no DTM to return
    dtm <- NULL
  } else {
    #otherwise (i.e., if studyplot was not NULL), set the dtm to the downloaded dtm (i.e., there is
    #something to return)
    dtm <- dtm
  }

  #if the user set lcp.dens to TRUE
  if(lcp.dens==TRUE & netw.type=="allpairs"){
    #there are data to return
    lcps_count <- lcps_count
    lcp.density.perc <- lcp.density.perc
  } else {
    #otherwise there is nothing to return
    lcps_count <- NULL
    lcp.density.perc <- NULL
  }

  results <- list("dtm"=dtm,
                  "LCPs.netw"=paths.netw,
                  "LCPs.netw.merged"=paths.netw.merged,
                  "LCPs.netw.neigh"=paths.netw.neigh,
                  "LCPs.netw.neigh.merged"=paths.netw.neigh.merged,
                  "LCPs.density.count"=lcps_count,
                  "LCPs.density.perc"=lcp.density.perc,
                  "cost.matrix.min"=cost.matrix.min,
                  "cost.matrix.hrs"=cost.matrix.hr,
                  "cost.matrix"=cost.matrix)

  return(results)
}
