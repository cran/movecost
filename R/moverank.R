#' R function for calculating sub-optimal least-cost paths bewteen an origin and a destination location
#'
#' The function provides the facility to calculate the LCP between an origin and a destination location and (more importantly) to
#' work out the first six sub-optimal LCPs between those locations. The underlying idea is the following: given two locations, we can
#' calculate the least-costly path between them; but, if we disregard that LCP, what path would be the second least costly?
#' And if we in turn disregard those first two, what the third least costly path would be? The same reasoning holds for all the subsequent n-th LCPs. Under the hood,
#' \code{moverank()} rests on \code{\link{movecost}} and implements the same cost functions. See the help documentation of \code{movecost()} for further details.\cr
#' Visit this \href{https://drive.google.com/file/d/1gLDrkZFh1b_glzCEqKdkPrer72JJ9Ffa/view?usp=sharing}{LINK} to access the package's vignette.\cr
#'
#' Internally, \code{moverank()} uses \code{movecost()} to generate the first (optimal) LCP. In a second iteration,
#' the optimal LCP is internally used as barrier (see \code{\link{movecost}}) when calculating the 2nd LCPs. Then, in a third iteration,
#' the two previously generated LCPs are used as barriers when working out the 3rd LCPs. The process repeats along the same lines until the
#' 6th LCP is calculated. The 1st LCP is deemed to represent the optimal path (cost-wise) between the two locations, while the
#' 2nd-to-5th LCPs are deemed to represent progressively sub-optimal paths.\cr
#'
#' It is worth noting that it may happen that some LCP will cross another one; this cannot be anticipated and is context dependent.
#' In those cases, the user may want to set the \code{move} parameter to 8 (see the section about
#' inhibition of movement in the help documentation of \code{\link{movecost}}). \cr
#'
#' The function provides the facility to render the LCPs either on the input DTM or on the least-cost corridor
#' between the two locations. The second option can be obtained by setting the \code{use.corr} parameter to \code{TRUE}.
#' Also, by setting the \code{add.chart} parameter to \code{TRUE}, the function renders a bubble chart that plots the LCPs length
#' against their rank, while the size of the bubbles is proportional to the cost. All the LCPs will be plotted.\cr
#'
#' By setting the \code{export} parameter to \code{TRUE}, the LCPs, the DTM (if acquired online), and the least-cost corridor
#' (if obtained by setting the \code{use.corr} parameter to \code{TRUE}), will be exported: the DTM and least-cost corridor as a raster layer,
#' the LCPs as shapefile layer. The LCPs and the least-cost corridor files will be given a suffix indicating which cost function has been used.\cr
#'
#' @param dtm Digital Terrain Model (RasterLayer class); if not provided, elevation data will be acquired online for the area enclosed by the 'studyplot' parameter (see \code{\link{movecost}}).
#' @param origin first location from which the least-cost corridor is calculated (SpatialPointsDataFrame class); if it contains more than two locations, see the 'Description' section above.
#' @param destin second location from which the least-cost corridor is calculated (SpatialPointsDataFrame class); if parameter 'a' stores more than two locations, this parameter is disregarded; see the 'Description' section above.
#' @param studyplot polygon (SpatialPolygonDataFrame class) representing the study area for which online elevation data are acquired (see \code{\link{movecost}}); NULL is default.
#' @param barrier area where the movement is inhibited (SpatialLineDataFrame or SpatialPolygonDataFrame class) (see \code{\link{movecost}}).
#' @param plot.barrier TRUE or FALSE (default) if the user wants or does not want the barrier to be plotted (see \code{\link{movecost}}).
#' @param irregular.dtm TRUE or FALSE (default) if the input DTM features irregular margins (see \code{\link{movecost}}).
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
#' @param time time-unit expressed by the accumulated raster if Tobler's and other time-related cost functions are used; h' for hour, 'm' for minutes.
#' @param lcp.n number of LCPs rendered in the output plot (min=1, max=6; 3 by default; the 1st LCP is the optimal one, while the LCPs from the 2nd to the 6th are the sub-optimal ones).
#' @param move number of directions in which cells are connected: 4 (rook's case), 8 (queen's case), 16 (knight and one-cell queen moves; default).
#' @param cogn.slp  TRUE or FALSE (default) if the user wants or does not want the 'cognitive slope' to be used in place of the real slope (see \code{\link{movecost}}).
#' @param sl.crit critical slope (in percent), typically in the range 8-16 (10 by default) (used by the wheeled-vehicle cost function; see \code{\link{movecost}}).
#' @param W walker's body weight (in Kg; 70 by default; used by the Pandolf's and Van Leusen's cost function; see \code{\link{movecost}}).
#' @param L carried load weight (in Kg; 0 by default; used by the Pandolf's and Van Leusen's cost function; see \code{\link{movecost}}).
#' @param N coefficient representing ease of movement (1 by default) (see \code{\link{movecost}}).
#' @param V speed in m/s (1.2 by default) (used by the Pandolf et al.'s, Pandolf et al.s with correction factor, Van Leusen's, and Ardigo et al.'s cost function; if set to 0, it is internally worked out on the basis of Tobler on-path hiking function (see \code{\link{movecost}}).
#' @param z zoom level for the elevation data downloaded from online sources (from 0 to 15; 9 by default) (see \code{\link{movecost}} and \code{\link[elevatr]{get_elev_raster}}).
#' @param use.corr TRUE or FALSE (default) is the user wants or does not want the least-cost corridor raster to be rendered in place of the input DTM.
#' @param leg.pos set the position of the legend in rendered plot; 'topright' by default (other options: "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right", "center").
#' @param leg.cex set the size of the labels used in the legend displayed in the rendered plot (0.55 by default).
#' @param transp set the transparency of the slopeshade raster that is plotted over the least-cost corridor raster (0.5 by default).
#' @param add.chart TRUE or FALSE (default) is the user wants or does not want a bubble chart visualising LCPs length vs rank vs cost to be rendered.
#' @param bubble.cex set the size of the labels reporting the LCPs cost in the bubble chart (0.5 by default).
#' @param export TRUE or FALSE (default) if the user wants or does not want the output to be exported; if TRUE, the least-cost corridor and the DTM (if not provided by the user but acquired online)
#' are expoerted as a GeoTiff file, while the LCPs as a shapefile layer. All the exported files (excluding the DTM) will bear a suffix corresponding to the cost function selected by the user.
#'
#' @return The function returns a list storing the following components \itemize{
##'  \item{dtm: }{Digital Terrain Model ('RasterLayer' class); returned only if not provided by the user and acquired online instead}
##'  \item{LCPs: }{least-cost paths ('SpatialLinesDataFrame' class) ranked from 1 (optimal) to 6 (sub-optimal LCPs)}
##'  \item{lc.corr: }{least-cost corridor between the origin and destination location ('RasterLayer' class); returned if the \code{use.corr} parameter is set to \code{TRUE}}
##' }
##'
#' @keywords moverank
#' @export
#' @importFrom raster crop raster terrain
#' @importFrom elevatr get_elev_raster
#' @importFrom grDevices grey
#' @importFrom graphics symbols text
#'
#' @examples
#' # load a sample Digital Terrain Model
#' data(volc)
#'
#'
#' # load the sample destination locations on the above DTM
#' data(destin.loc)
#'
#' # calculate the optimal and sub-optimals LCPs between two locations
#' #result <- moverank(volc, destin.loc[1,], destin.loc[4,], move=8, funct="t")
#'
#' @seealso \code{\link{movecost}}
#'
#'
moverank <- function (dtm=NULL, origin, destin, studyplot=NULL, barrier=NULL, plot.barrier=FALSE, irregular.dtm=FALSE, funct="t", time="h", lcp.n=3, move=16, cogn.slp=FALSE, sl.crit=10, W=70, L=0, N=1, V=1.2, z=9, use.corr=FALSE, leg.pos="topright", leg.cex=0.55, add.chart=FALSE, bubble.cex=0.50, transp=0.5, export=FALSE){

  #if no dtm is provided
  if (is.null(dtm)==TRUE) {
    #get the elvation data using the elevatr's get_elev_raster() function, using the studyplot dataset (SpatialPolygonDataFrame)
    #to select the area whose elevation data are to be downloaded;
    #z sets the resolution of the elevation datataset
    elev.data <- elevatr::get_elev_raster(studyplot, z = z, verbose=FALSE, override_size_check = TRUE)
    #crop the elevation dataset to the exact boundary of the studyplot dataset
    dtm <- raster::crop(elev.data, studyplot)
  }

  #create an empty list to store the different LCPs
  result <- list()

  #in 'moverank', 'field' is not entered by the user;
  #it is defined here, and is a small value close to 0 but not equal to 0, otherwise the calculation of the LCPs CANNOT be
  #carried out
  field <- 0.01

  #IF the barrier IS provided...
  if(is.null(barrier)==FALSE) {
    #if the barrier is of class SpatialPolygonsDataFrame, convert it to SpatialLinesDataFrame because in the code
    #further below it is not possible to add together layers of different type
    if(class(barrier)[1]=="SpatialPolygonsDataFrame") {
      barrier <- as(barrier, "SpatialLinesDataFrame")
    }
    #...create an empty list of a copy of the barrier layer
    barrier.copy <- list()
    #create 5 slots within the list, each one containing a copy of the barrier;
    #those five elements will be used (one after one) later on, adding each one iteratively to the various
    #calculared LCPs
    for(i in 1:5){
      barrier.copy[[i]] <- barrier
    }

    message("Wait...calculating the optimal LCP...")

    #calculate the optimal LCP; if a barrier is provided, pass the barrier to the 'movecost()' function
    output.data <- movecost(dtm=dtm, origin=origin, destin=destin, barrier=barrier, irregular.dtm=irregular.dtm, funct=funct, time=time, field=field, move=move, cogn.slp=cogn.slp, sl.crit=sl.crit, W=W, L=L, N=N, V=V, z=z, return.base=FALSE, graph.out=FALSE)
    result[[1]] <- output.data$LCPs
    result[[1]]$rank <- 1
    result[[1]]$cost <- output.data$dest.loc.w.cost$cost

    message("Wait...calculating sub-optimal LCP...1/5")

    #calculate the 2nd (sub-optimal) LCP using the previous 1st (optimal) LCP as barrier;
    #note that the first slot of the barrier.copy list is used.
    output.data <- movecost(dtm=dtm, origin=origin, destin=destin, barrier=barrier.copy[[1]]+result[[1]], irregular.dtm=irregular.dtm, funct=funct, time=time, field=field, move=move, cogn.slp=cogn.slp, sl.crit=sl.crit, W=W, L=L, N=N, V=V, z=z, return.base=FALSE, graph.out=FALSE)
    result[[2]] <- output.data$LCPs
    result[[2]]$rank <- 2
    result[[2]]$cost <- output.data$dest.loc.w.cost$cost

    message("Wait...calculating sub-optimal LCP...2/5")

    #same as above, using the 1st and 2nd LCPs as barriers...
    #note that the second slot of the barrier.copy list is used.
    output.data <- movecost(dtm=dtm, origin=origin, destin=destin, barrier=barrier.copy[[2]]+result[[1]]+result[[2]], irregular.dtm=irregular.dtm, funct=funct, time=time, field=field, move=move, cogn.slp=cogn.slp, sl.crit=sl.crit, W=W, L=L, N=N, V=V, z=z, return.base=FALSE, graph.out=FALSE)
    result[[3]] <- output.data$LCPs
    result[[3]]$rank <- 3
    result[[3]]$cost <- output.data$dest.loc.w.cost$cost

    message("Wait...calculating sub-optimal LCP...3/5")

    #same as above, using the 1st, 2nd, and 3rd LCPs as barriers...
    #note that the third slot of the barrier.copy list is used.
    output.data <- movecost(dtm=dtm, origin=origin, destin=destin, barrier=barrier.copy[[3]]+result[[1]]+result[[2]]+result[[3]], irregular.dtm=irregular.dtm, funct=funct, time=time, field=field, move=move, cogn.slp=cogn.slp, sl.crit=sl.crit, W=W, L=L, N=N, V=V, z=z, return.base=FALSE, graph.out=FALSE)
    result[[4]] <- output.data$LCPs
    result[[4]]$rank <- 4
    result[[4]]$cost <- output.data$dest.loc.w.cost$cost

    message("Wait...calculating sub-optimal LCP...4/5")

    output.data <- movecost(dtm=dtm, origin=origin, destin=destin, barrier=barrier.copy[[4]]+result[[1]]+result[[2]]+result[[3]]+result[[4]], irregular.dtm=irregular.dtm, funct=funct, time=time, field=field, move=move, cogn.slp=cogn.slp, sl.crit=sl.crit, W=W, L=L, N=N, V=V, z=z, return.base=FALSE, graph.out=FALSE)
    result[[5]] <- output.data$LCPs
    result[[5]]$rank <- 5
    result[[5]]$cost <- output.data$dest.loc.w.cost$cost

    message("Wait...calculating sub-optimal LCP...5/5")

    output.data <- movecost(dtm=dtm, origin=origin, destin=destin, barrier=barrier.copy[[5]]+result[[1]]+result[[2]]+result[[3]]+result[[4]]+result[[5]], irregular.dtm=irregular.dtm, funct=funct, time=time, field=field, move=move, cogn.slp=cogn.slp, sl.crit=sl.crit, W=W, L=L, N=N, V=V, z=z, return.base=FALSE, graph.out=FALSE)
    result[[6]] <- output.data$LCPs
    result[[6]]$rank <- 6
    result[[6]]$cost <- output.data$dest.loc.w.cost$cost
  }

  #IF the barrier is NOT provided...
  if(is.null(barrier)==TRUE) {

    message("Wait...calculating the optimal LCP...")

    #calculate the optimal LCP
    output.data <- movecost(dtm=dtm, origin=origin, destin=destin, barrier=NULL, irregular.dtm=irregular.dtm, funct=funct, time=time, field=field, move=move, cogn.slp=cogn.slp, sl.crit=sl.crit, W=W, L=L, N=N, V=V, z=z, return.base=FALSE, graph.out=FALSE)
    result[[1]] <- output.data$LCPs
    result[[1]]$rank <- 1
    result[[1]]$cost <- output.data$dest.loc.w.cost$cost

    message("Wait...calculating sub-optimal LCP...1/5")

    #calculate the 2nd (sub-optimal) LCP using the previous 1st (optimal) LCP as barrier;
    output.data <- movecost(dtm=dtm, origin=origin, destin=destin, barrier=result[[1]], irregular.dtm=irregular.dtm, funct=funct, time=time, field=field, move=move, cogn.slp=cogn.slp, sl.crit=sl.crit, W=W, L=L, N=N, V=V, z=z, return.base=FALSE, graph.out=FALSE)
    result[[2]] <- output.data$LCPs
    result[[2]]$rank <- 2
    result[[2]]$cost <- output.data$dest.loc.w.cost$cost

    message("Wait...calculating sub-optimal LCP...2/5")

    output.data <- movecost(dtm=dtm, origin=origin, destin=destin, barrier=result[[1]]+result[[2]], irregular.dtm=irregular.dtm, funct=funct, time=time, field=field, move=move, cogn.slp=cogn.slp, sl.crit=sl.crit, W=W, L=L, N=N, V=V, z=z, return.base=FALSE, graph.out=FALSE)
    result[[3]] <- output.data$LCPs
    result[[3]]$rank <- 3
    result[[3]]$cost <- output.data$dest.loc.w.cost$cost

    message("Wait...calculating sub-optimal LCP...3/5")

    output.data <- movecost(dtm=dtm, origin=origin, destin=destin, barrier=result[[1]]+result[[2]]+result[[3]], irregular.dtm=irregular.dtm, funct=funct, time=time, field=field, move=move, cogn.slp=cogn.slp, sl.crit=sl.crit, W=W, L=L, N=N, V=V, z=z, return.base=FALSE, graph.out=FALSE)
    result[[4]] <- output.data$LCPs
    result[[4]]$rank <- 4
    result[[4]]$cost <- output.data$dest.loc.w.cost$cost

    message("Wait...calculating sub-optimal LCP...4/5")

    output.data <- movecost(dtm=dtm, origin=origin, destin=destin, barrier=result[[1]]+result[[2]]+result[[3]]+result[[4]], irregular.dtm=irregular.dtm, funct=funct, time=time, field=field, move=move, cogn.slp=cogn.slp, sl.crit=sl.crit, W=W, L=L, N=N, V=V, z=z, return.base=FALSE, graph.out=FALSE)
    result[[5]] <- output.data$LCPs
    result[[5]]$rank <- 5
    result[[5]]$cost <- output.data$dest.loc.w.cost$cost

    message("Wait...calculating sub-optimal LCP...5/5")

    output.data <- movecost(dtm=dtm, origin=origin, destin=destin, barrier=result[[1]]+result[[2]]+result[[3]]+result[[4]]+result[[5]], irregular.dtm=irregular.dtm, funct=funct, time=time, field=field, move=move, cogn.slp=cogn.slp, sl.crit=sl.crit, W=W, L=L, N=N, V=V, z=z, return.base=FALSE, graph.out=FALSE)
    result[[6]] <- output.data$LCPs
    result[[6]]$rank <- 6
    result[[6]]$cost <- output.data$dest.loc.w.cost$cost
  }

  #string together the LCPs calculated so far
  merged.lcps <- do.call(rbind, result)

  #define the appropriate text to be used for subsequent plotting
  if (funct=="t") {
    #set the labels to be used within the returned plot
    sub.title <- paste0("LCPs based on the Tobler's on-path hiking function \n terrain factor N=", N)
    legend.lab <- "cost in walking-time"
  }

  if (funct=="tofp") {
    #set the labels to be used within the returned plot
    sub.title <- "LCPs based on the Tobler's off-path hiking function"
    legend.lab <- "cost in walking-time"
  }

  if(funct=="mp") {
    #set the labels to be used within the returned plot
    sub.title <- paste0("LCPs on the Marquez-Perez et al.'s modified Tobler hiking function \n terrain factor N=", N)
    legend.lab <- "cost in walking-time"
  }

  if(funct=="icmonp") {
    #set the labels to be used within the returned plot
    sub.title <- paste0("LCPs on the (male, on-path) Irmischer-Clarke's hiking function \n terrain factor N=", N)
    legend.lab <- "cost in walking-time"
  }

  if(funct=="icmoffp") {
    #set the labels to be used within the returned plot
    sub.title <- "LCPs based on the (male, off-path) Irmischer-Clarke's hiking function"
    legend.lab <- "cost in walking-time"
  }

  if(funct=="icfonp") {
    #set the labels to be used within the returned plot
    sub.title <- paste0("LCPs based on the (female, on-path) Irmischer-Clarke's hiking function \n terrain factor N=", N)
    legend.lab <- "cost in walking-time"
  }

  if(funct=="icfoffp") {
    sub.title <- "LCPs based on the (female, off-path) Irmischer-Clarke's hiking function"
    legend.lab <- "cost in walking-time"
  }

  if(funct=="gkrs") {
    sub.title <- paste0("LCPs based on the Garmy et al.'s hiking function \n terrain factor N=", N)
    legend.lab <- "cost in walking-time"
  }

  if(funct=="r") {
    sub.title <- paste0("LCPs based on the Rees' hiking function \n terrain factor N=", N)
    legend.lab <- "cost in walking-time"
  }

  if(funct=="ug") {
    #set the labels to be used within the returned plot
    sub.title <- paste0("LCPs based on the Uriarte Gonzalez's hiking function \n terrain factor N=", N)
    legend.lab <- "cost in walking-time"
  }

  if(funct=="ks") {
    #set the labels to be used within the returned plot
    sub.title <- paste0("LCPs based on the Kondo-Seino's hiking function \n terrain factor N=", N)
    legend.lab <- "cost in walking-time"
  }

  if(funct=="ree") {
    #set the labels to be used within the returned plot
    sub.title <- paste0("LCPs based on the relative energetic expenditure cost function \n terrain factor N=", N)
    legend.lab <- "cost in relative energetic expenditure"
  }

  if(funct=="hrz") {
    #set the labels to be used within the returned plot
    sub.title <- paste0("LCPs based on the Herzog's metabolic cost function \n terrain factor N=", N)
    legend.lab <- "cost in J / (Kg*m)"
  }

  if(funct=="wcs") {
    #set the labels to be used within the returned plot
    sub.title <- paste0("LCPs based on the wheeled-vehicle critical slope cost function \ncritical slope set to ", sl.crit, " percent; terrain factor N=", N)
    legend.lab <- "abstract cost"
  }

  if(funct=="vl") {
    #set the labels to be used within the returned plot
    if (V==0) {
      sub.title <- paste0("LCPs based on the Van Leusen's metabolic energy expenditure cost function \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V is based on the Tobler on-path hiking function")
    } else {
      sub.title <- paste0("LCPs based on the Van Leusen's metabolic energy expenditure cost function \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V: ", V)
    }
    legend.lab <- "cost in Megawatts"
  }

  if(funct=="p") {
    #set the labels to be used within the returned plot
    if (V==0) {
      sub.title <- paste0("LCPs based on the Pandolf et al.'s metabolic energy expenditure cost function \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V is based on the Tobler on-path hiking function")
    } else {
      sub.title <- paste0("LCPs based on the Pandolf et al.'s metabolic energy expenditure cost function \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V: ", V)
    }
    legend.lab <- "cost in Megawatts"
  }

  if(funct=="pcf") {
    #set the labels to be used within the returned plot
    if (V==0) {
      sub.title <- paste0("LCPs based on the Pandolf et al.'s metabolic energy expenditure cost function with correction factor \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V is based on the Tobler on-path hiking function")
    } else {
      sub.title <- paste0("LCPs based on the Pandolf et al.'s metabolic energy expenditure cost function with correction factor \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V: ", V)
    }
    legend.lab <- "cost in Megawatts"
  }

  if(funct=="ls") {
    #set the labels to be used within the returned plot
    sub.title <- paste0("LCPs based on the Llobera-Sluckin's metabolic energy expenditure cost function \n terrain factor N=", N)
    legend.lab <- "cost in KJ/m"
  }

  if(funct=="alb") {
    #set the labels to be used within the returned plot
    sub.title <- "LCPs based on the Alberti's modified Tobler hiking function"
    legend.lab <- "cost in walking-time"
  }

  if(funct=="b") {
    #set the labels to be used within the returned plot
    sub.title <- paste0("LCPs based on the Bellavia's cost function \n terrain factor N=", N)
    legend.lab <- "abstract cost"
  }

  if(funct=="m") {
    #set the labels to be used within the returned plot
    sub.title <- paste0("LCPs based on the Minetti et al.'s metabolic cost function \n terrain factor N=", N)
    legend.lab <- "cost in J / (Kg*m)"
  }

  if (funct=="ma") {
    #set the labels to be used within the returned plot
    sub.title <- paste0("LCPs based on the Marin Arroyo's hiking function \n terrain factor N=", N)
    legend.lab <- "cost in walking time"
  }

  if(funct=="a") {
    #set the labels to be used within the returned plot
    if (V==0) {
      sub.title <- paste0("LCPs based on the Ardigo et al.'s metabolic cost function \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V is based on the Tobler on-path hiking function")
    } else {
      sub.title <- paste0("LCPs based on the Ardigo et al.'s metabolic cost function \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V: ", V)
    }
    legend.lab <- "cost in J / (Kg*m)"
  }

  if(funct=="h") {
    #set the labels to be used within the returned plot
    sub.title <- paste0("Cost based on the Hare's metabolic cost function \n cost in cal/km \n terrain factor N=", N)
    legend.lab <- "cost in cal/km"
  }

  if(funct=="e") {
    #set the labels to be used within the returned plot
    sub.title <- paste0("Cost based on the Eastman's cost function \n terrain factor N=", N)
    legend.lab <- "abstract cost"
  }

  if (funct=="trp") {
    #set the labels to be used within the returned plot
    sub.title <- paste0("LCPs based on the Tripcevich's hiking function \n terrain factor N=", N)
    legend.lab <- "cost in walking time"
  }

  sub.title <- paste0(sub.title, "\nblack dot=start location; red dot=destination location")

  #produce the ingredient for the slopeshade raster
  #to be used in the rendered plots
  slope <- raster::terrain(dtm, opt = "slope")

  #if use.corr is FALSE...
  if(use.corr==FALSE) {
    #...plot the DTM raster
    raster::plot(dtm,
                 main="Ranked Least-Cost Paths",
                 sub=sub.title,
                 cex.main=0.95,
                 cex.sub=0.75,
                 legend.lab="Elevation (masl)")
  } else {
    #...otherwise calculate the least-cost corridor...
    message("Wait...calculating the least-cost corridor")
    lc.corr <- movecorr(dtm=dtm, a=origin, b=destin, funct=funct, time=time, move=move, barrier=barrier, irregular.dtm=irregular.dtm, cogn.slp=cogn.slp, sl.crit=sl.crit, W=W, L=L, N=N, V=V, z=z, transp=transp, graph.out = FALSE)$lc.corridor
    #...and plot it
    raster::plot(lc.corr,
                 main="Least-cost corridor with ranked Least-Cost Paths",
                 sub=sub.title,
                 cex.main=0.95,
                 cex.sub=0.75,
                 legend.lab=legend.lab)
  }

  #add the slopeshade
  raster::plot(slope,
               col = rev(grey(0:100/100)),
               legend = FALSE,
               alpha=transp,
               add=TRUE)

  #add the (merged) LCPs
  raster::plot(merged.lcps[1:lcp.n,],
               add=T,
               lty=as.numeric(as.factor(merged.lcps$rank)))

  #add the origin
  raster::plot(origin,
               pch=20,
               add=TRUE)

  #add the destinations
  raster::plot(destin,
               pch=20,
               col="red",
               add=TRUE)

  #add the legend
  legend(leg.pos,
         legend=1:lcp.n,
         lty=1:lcp.n,
         cex=leg.cex)

  #if the barrier is provided AND if plot.barrier is TRUE, add the barrier
  if(is.null(barrier)==FALSE & plot.barrier==TRUE) {
    raster::plot(barrier, col="blue", add=TRUE)
  }

  if (add.chart==TRUE) {
    radius <- sqrt(merged.lcps$cost/pi)
    graphics::symbols(merged.lcps$length, merged.lcps$rank, circles=radius,
                      inches=0.35,
                      fg="white",
                      bg="red",
                      xlab="LCPs length",
                      ylab="LCPs rank",
                      cex.lab=0.95,
                      main="Bubble chart of LCPs lenght vs rank", cex.main=0.85,
                      sub=paste0("bubbles size is proportional to ", legend.lab), cex.sub=0.70)
    graphics::text(merged.lcps$length, merged.lcps$rank, round(merged.lcps$cost, 3), cex=bubble.cex)
  }

  #if export is TRUE, export the LCPs as a shapefile
  if(export==TRUE){
    rgdal::writeOGR(merged.lcps, ".", paste0("LCPs_", funct), driver="ESRI Shapefile")
  }


  #if no DTM was provided (i.e., if 'studyplot' is not NULL), export the downloaded DTM as a raster file
  if(export==TRUE & is.null(studyplot)==FALSE){
    raster::writeRaster(dtm, "dtm", format="GTiff")
  }

  #if export is TRUE and use.corr is TRUE, export least-cost corridor
  if(export==TRUE & use.corr==TRUE){
    raster::writeRaster(lc.corr, paste0("lc_corridor_", funct), format="GTiff")
  }

  #if studyplot is NULL (i.e., if the DTM has been provided)...
  if(is.null(studyplot)==TRUE){
    #set dtm to NULL, i.e. there is no DTM to return
    dtm <- NULL
  } else{
    #otherwise (i.e., if studyplot was not NULL), set the dtm to the downloaded dtm (i.e., there is
    #something to return)
    dtm <- dtm
  }

  #same as the step above with regards to the least-cost corridor
  if(use.corr==TRUE){
    lc.corr <- lc.corr
  } else {
    lc.corr  <- NULL
  }

  results <- list("dtm"=dtm,
                  "LCPs"=merged.lcps,
                  "least-cost corridor"=lc.corr)

  return(results)
}
