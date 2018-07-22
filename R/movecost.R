#' R function for calculating accumulated cost of movement across the terrain and least-cost paths from an origin
#'
#' The function provides the facility to calculate the accumulated cost of movement around a starting location and to optionally calculate least-cost path(s) toward
#' one or multiple destinations. It implements different cost estimations related to human movement across the landscape.
#' The function takes as input a Digital Terrain Model ('RasterLayer' class) and a point feature ('SpatialPointsDataFrame' class), the latter representing
#' the starting location, i.e. the location from which the accumulated cost is calculated. \cr
#'
#' If the parameter 'destin' is fed with a dataset representing destination location(s) ('SpatialPointsDataFrame' class), the function will also calculate
#' least-cost path(s) plotted on the input DTM; the length of each path will be saved under the variable 'length' stored in the 'LCPs' dataset ('SpatialLines' class) returned by the function.
#' The red dot(s) representing the destination location(s) will be labelled with numeric values representing
#' the cost value at the location(s). The cost value will be also appended to the updated destination dataset returned by the function and
#' storing a new variable named 'cost'.\cr
#'
#' The function builds on functions out of Jacob van Etten's 'gdistance' package, and by default uses a 16-directions movement in calculating the accumulated cost-surface.
#' The number of movements can be optionally set by the user via the 'moves' parameter.
#'
#' The slope is internally calculated by the function using the 'terrain()' function out of the 'raster' package, using 8 neighbors.
#' The slope is calculated in degrees and, for the sake of its use within the implemented cost functions, it is internally modified (i.e., turned into either gradient or percent)
#' according to the user-selected cost function. The user can input a slope dataset (RasterLayer class) using the parameter 'slope'; this can prove usefull
#' in cases when the slope has to be preliminarily modified, for instance to mask out areas that cannot be crossed or to weight the slope values according
#' to a user-defined weighting scheme.\cr
#'
#' The following cost functions are implemented (\strong{x} stands for \strong{slope}):\cr
#'
#' \strong{Tobler's hiking function (on-path) (speed in kmh)}:\cr
#'
#' \eqn{ 6 * exp(-3.5 * abs(tan(x*pi/180) + 0.05)) }\cr
#'
#'
#' \strong{Tobler's hiking function (off-path) (speed in kmh)}:\cr
#'
#' \eqn{ (6 * exp(-3.5 * abs(tan(x*pi/180) + 0.05))) * 0.6 }\cr
#'
#' as per Tobler's indication, the off-path walking speed is reduced by 0.6.\cr
#'
#'
#' \strong{Marquez-Perez et al.'s modified Tobler hiking function (speed in kmh)}:\cr
#'
#' \eqn{ 4.8 * exp(-5.3 * abs((tan(x*pi/180) * 0.7) + 0.03)) }\cr
#'
#' modified version of the Tobler's hiking function as proposed by Joaquin Marquez-Parez, Ismael Vallejo-Villalta & Jose I. Alvarez-Francoso (2017), "Estimated travel time for walking trails in natural areas",
#' Geografisk Tidsskrift-Danish Journal of Geography, 117:1, 53-62, DOI: 10.1080/00167223.2017.1316212.\cr
#'
#'
#' \strong{Irmischer-Clarke's modified Tobler hiking function (on-path)}:\cr
#'
#' \eqn{ (0.11 + exp(-(tan(x*pi/180)*100 + 5)^2 / (2 * 30)^2)) * 3.6 }\cr
#'
#' modified version of the Tobler's function as proposed for (male) on-path hiking by Irmischer, I. J., & Clarke, K. C. (2018). Measuring and modeling the speed of human navigation.
#' Cartography and Geographic Information Science, 45(2), 177-186. https://doi.org/10.1080/15230406.2017.1292150. \strong{Note}: the function originally expresses speed in m/s; it has been is reshaped (multiplied by 3.6)
#' to turn it into kmh for consistency with the other Tobler-related cost functions; also, originally the slope is in percent; \eqn{tan(x*pi/180)*100} turns the slope from degrees to percent (=rise/run*100).\cr
#'
#'
#'\strong{Irmischer-Clarke's modified Tobler hiking function (off-path)}:\cr
#'
#' \eqn{ (0.11 + 0.67 * exp(-(tan(x*pi/180)*100 + 2)^2 / (2 * 30)^2)) * 3.6 }\cr
#'
#'
#'\strong{Uriarte Gonzalez's slope-dependant walking-time cost function}:\cr
#'
#' \eqn{ 1/ (0.0277 * (tan(x*pi/180)*100) + 0.6115) }\cr
#'
#' proposed by Uriarte Gonzalez;
#' \strong{see}: Chapa Brunet, T., Garcia, J., Mayoral Herrera, V., & Uriarte Gonzalez, A. (2008). GIS landscape models for the study of preindustrial settlement patterns in Mediterranean areas.
#' In Geoinformation Technologies for Geo-Cultural Landscapes (pp. 255-273). CRC Press. https://doi.org/10.1201/9780203881613.ch12.\cr
#' The cost function is originally expressed in seconds; for the purpose of its implementation in this function, it is the reciprocal of time (1/T) that is used in order to eventually get
#' T/1. Also, originally the slope is in percent: \eqn{tan(x*pi/180)*100} turns the slope from degrees to percent (=rise/run*100).
#' Unlike the original cost function, here the pixel resolution is not taken into account since 'gdistance' takes care of the cells' dimension
#' when calculating accumulated costs.
#'
#'
#' \strong{Relative energetic expenditure cost function}:\cr
#'
#' \eqn{ 1 / (tan(x*pi/180) / tan (1*pi/180)) }\cr
#'
#' slope-based cost function expressing change in potential energy expenditure;
#' \strong{see} Conolly, J., & Lake, M. (2006). Geographic Information Systems in Archaeology. Cambridge: Cambridge University Press, p. 220;
#' \strong{see also} Newhard, J. M. L., Levine, N. S., & Phebus, A. D. (2014). The development of integrated terrestrial and marine pathways in the Argo-Saronic region, Greece. Cartography and Geographic Information Science, 41(4), 379-390, with references to studies that use this
#' function; \strong{see also} ten Bruggencate, R. E., Stup, J. P., Milne, S. B., Stenton, D. R., Park, R. W., & Fayek, M. (2016). A human-centered GIS approach to modeling mobility on southern Baffin Island, Nunavut,
#' Canada. Journal of Field Archaeology, 41(6), 684-698. https://doi.org/10.1080/00934690.2016.1234897.\cr
#'
#'
#' \strong{Herzog's metabolic cost function in J/(kg*m)}:\cr
#'
#' \eqn{ 1 / ((1337.8 * tan(x*pi/180)^6 + 278.19 * tan(x*pi/180)^5 - 517.39 * tan(x*pi/180)^4 - 78.199 * tan(x*pi/180)^3 + 93.419 * tan(x*pi/180)^2 + 19.825 * tan(x*pi/180) + 1.64)) }\cr
#'
#' \strong{see} Herzog, I. (2016). Potential and Limits of Optimal Path Analysis. In A. Bevan & M. Lake (Eds.), Computational Approaches to Archaeological Spaces (pp. 179-211). New York: Routledge.\cr
#'
#'
#' \strong{Wheeled-vehicle critical slope cost function}:\cr
#'
#' \eqn{ 1 / (1 + ((tan(x*pi/180)*100) / sl.crit)^2)  }\cr
#'
#' where \eqn{sl.crit} (=critical slope, in percent) is "the transition where switchbacks become more effective than direct uphill or downhill paths" and typically is in the range 8-16;
#' \strong{see} Herzog, I. (2016). Potential and Limits of Optimal Path Analysis. In A. Bevan & M. Lake (Eds.), Computational Approaches to Archaeological Spaces (pp. 179-211). New York: Routledge. \cr
#'
#'
#' \strong{Pandolf et al.'s metabolic energy expenditure cost function (in Watts)}:\cr
#'
#' \eqn{ 1 / (1.5 * W + 2.0 * (W + L) * (L / W)^2 + N * (W + L) * (1.5 * V^2 + 0.35 * V * (tan(x*pi/180)*100))) }\cr
#'
#' where \eqn{W} is the walker's body weight (Kg), \eqn{L} is the carried load (in Kg), \eqn{V} is the velocity in m/s, \eqn{N} is a coefficient representing ease of movement on the terrain.\cr
#'
#' As for the latter, suggested values available in literature are: Asphalt/blacktop=1.0; Dirt road=1.1; Grass=1.1; Light brush=1.2; Heavy brush=1.5; Swampy bog=1.8; Loose sand=2.1; Hard-packed snow=1.6; Ploughed field=1.3;
#' \strong{see} de Gruchy, M., Caswell, E., & Edwards, J. (2017). Velocity-Based Terrain Coefficients for Time-Based Models of Human Movement. Internet Archaeology, 45(45). https://doi.org/10.11141/ia.45.4.\cr
#'
#' For this cost function, \strong{see} Pandolf, K. B., Givoni, B., & Goldman, R. F. (1977). Predicting energy expenditure with loads while standing or walking very slowly. Journal of Applied Physiology, 43(4), 577-581. https://doi.org/10.1152/jappl.1977.43.4.577.\cr
#'
#' For the use of this cost function in a case study, \strong{see} Rademaker, K., Reid, D. A., & Bromley, G. R. M. (2012). Connecting the Dots: Least Cost Analysis, Paleogeography, and the Search for Paleoindian Sites in Southern Highland Peru. In D. A. White & S. L. Surface-Evans (Eds.), Least Cost Analysis of Social Landscapes. Archaeological Case Studies (pp. 32-45). University of Utah Press;
#' \strong{see also} Herzog, I. (2013). Least-cost Paths - Some Methodological Issues, Internet Archaeology 36 (http://intarch.ac.uk/journal/issue36/index.html) with references.\cr
#'
#' \strong{Note}: in the returned charts, the cost is transposed from Watts to Megawatts (see, e.g., Rademaker et al 2012 cited above).\cr
#'
#'
#' \strong{Van Leusen's metabolic energy expenditure cost function (in Watts)}:\cr
#'
#' \eqn{ 1 / (1.5 * W + 2.0 * (W + L) * (L / W)^2 + N * (W + L) * (1.5 * V^2 + 0.35 * V * (tan(x*pi/180)*100) + 10))  }\cr
#'
#' which modifies the Pandolf et al.'s equation; \strong{see} Van Leusen, P. M. (2002). Pattern to process: methodological investigations into the formation and interpretation of spatial patterns in archaeological landscapes. University of Groningen.\cr
#' \strong{Note} that, as per Herzog, I. (2013). Least-cost Paths - Some Methodological Issues, Internet Archaeology 36 (http://intarch.ac.uk/journal/issue36/index.html) and
#' unlike Van Leusen (2002), in the above equation slope is expressed in percent and speed in m/s; also, in the last bit of the equantion, 10 replaces
#' the value of 6 used by Van Leusen (as per Herzog 2013).\cr
#' \strong{Note}: in the returned charts, the cost is transposed from Watts to Megawatts.\cr
#'
#' \strong{Note} that the walking-speed-related cost functions listed above are used as they are, while the other functions are reciprocated.
#' This is done since "gdistance works with conductivity rather than the more usual approach using costs"; therefore
#' "we need inverse cost functions" (Nakoinz-Knitter (2016). "Modelling Human Behaviour in Landscapes". New York: Springer, p. 183).
#'  As a consequence, if we want to estimate time, we have to use the walking-speed functions as they are since the final accumulated values will correspond to the
#'  reciprocal of speed, i.e. pace. In the other cases, we have to use 1/cost-function to eventually get cost-function/1.\cr
#'
#' When using the Tobler-related cost functions, the time unit can be selected by the user setting the 'time' parameter to 'h' (hour) or to 'm' (minutes).\cr
#'
#' In general, the user can also select which type of visualization the function has to produce; this is achieved setting the 'outp' parameter to either 'r' (=raster)
#' or to 'c' (=contours). The former will produce a raster image with a colour scale and contour lines representing the accumulated cost surface; the latter parameter will only
#' produce contour lines.\cr
#'
#' The contour lines' interval is set using the parameter 'breaks'; if no value is passed to the parameter, the interval will be set by default to
#' 1/10 of the range of values of the accumulated cost surface.\cr
#'
#' @param dtm digital terrain model (RasterLayer class).
#' @param slope slope dataset (RasterLayer class); if NULL (default), the slope (in degree) is internally calculated (see Details).
#' @param origin location from which the walking time is computed (SpatialPointsDataFrame class).
#' @param destin location(s) to which least-cost path(s) is calculated (SpatialPointsDataFrame class).
#' @param funct cost function to be used: \strong{t} (default) uses the on-path Tobler's hiking function;
#' \strong{tofp} uses the off-path Tobler's hiking function; \strong{mt} uses the modified Tobler's function;
#' \strong{ic} uses the Irmischer-Clarke's modified Tobler's hiking function (on-path);
#' \strong{icofp} uses the Irmischer-Clarke's modified Tobler's hiking function (off-path);
#' \strong{ug} uses the Uriarte Gonzalez's slope-dependant walking-time cost function;
#' \strong{ree} uses the relative energetic expenditure cost function;
#' \strong{hrz} uses the Herzog's metabolic cost function;
#' \strong{wcs} uses the wheeled-vehicle critical slope cost function;
#' \strong{p} uses the Pandolf et al.'s metabolic energy expenditure cost function;
#' \strong{vl} uses the Van Leusen's metabolic energy expenditure cost function (see Details).
#' @param time time-unit expressed by the accumulated raster and by the isolines if Tobler's and Tobler-related cost functions are used;
#' 'h' for hour, 'm' for minutes.
#' @param outp type of output: 'raster' or 'contours' (see Details).
#' @param sl.crit critical slope (in percent), typically in the range 8-16 (10 by default) (used by the wheeled-vehicle cost function; see Details).
#' @param W walker's body weight (in Kg; used by the Pandolf's and Van Leusen's cost function; see Details).
#' @param L carried load weight (in Kg; used by the Pandolf's and Van Leusen's cost function; see Details).
#' @param N coefficient representing ease of movement (1 by default) (used by the Pandolf's and Van Leusen's cost function; see Details).
#' @param V speed in m/s (1.2 by default) (used by the Pandolf's and Van Leusen's cost function; see Details).
#' @param moves number of directions used when computing the accumulated cost-surface (16 by default).
#' @param breaks isolines interval; if no value is supplied, the interval is set by default to 1/10 of the range of values of the accumulated cost surface.
#' @param cont.lab if set to TRUE (default) display the labels of the contours over the accumulated cost surface.
#' @param destin.lab if set to TRUE (default) display the label(s) indicating the cost at the destination location(s).
#' @param cex.breaks set the size of the time labels used in the isochrones plot (0.6 by default).
#' @param cex.lcp.lab set the size of the labels used in least-cost path(s) plot (0.6 by default).
#' @param oneplot TRUE (default) or FALSE if the user wants or does not want the plots displayed in a single window.
#' @param export TRUE or FALSE (default) if the user wants or does not want the outputs to be exported; if TRUE, the accumulated cost surface will be
#' exported as a GeoTiff file, while the isolines and the least-cost path(s) will be exported as shapefile; all the exported files will bear a suffix corresponding
#' to the cost function selected by the user.
#' @return The function returns a list storing the following components \itemize{
##'  \item{accumulated.cost.raster: }{raster representing the accumualted cost ('RasterLayer' class)}
##'  \item{isolines: }{contour lines derived from the accumulated cost surface ('SpatialLinesDataFrame' class)}
##'  \item{LCPs: }{estimated least-cost paths ('SpatialLines' class)}
##'  \item{LCPs$length: }{length of each least-cost path (units depend on the unit used in the input DTM)}
##'  \item{dest.loc.w.cost: }{copy of the input destination location(s) dataset with a new variable ('cost') added}
##' }
#' @keywords movecost
#' @export
#' @importFrom grDevices terrain.colors
#' @importFrom graphics layout par
#' @examples
#' # load a sample Digital Terrain Model
#' volc <- raster::raster(system.file("external/maungawhau.grd", package="gdistance"))
#'
#' # load a sample start location on the above DTM
#' data(volc.loc)
#'
#' # load the sample destination locations on the above DTM
#' data(destin.loc)
#'
#' # calculate walking-time isochrones based on the on-path Tobler's hiking function,
#' # setting the time unit to hours and the isochrones interval to 0.05 hour;
#' # also, since destination locations are provided,
#' # least-cost paths from the origin to the destination locations will be calculated
#' # and plotted
#' result <- movecost(dtm=volc,origin=volc.loc, destin=destin.loc, breaks=0.05)
#'
movecost <- function (dtm, slope=NULL, origin, destin=NULL, funct="t", time="h", outp="r", sl.crit=10, W=70, L=0, N=1, V=1.2, moves=16, breaks=NULL, cont.lab=TRUE, destin.lab=TRUE, cex.breaks=0.6, cex.lcp.lab=0.6, oneplot=TRUE, export=FALSE){

  #calculate the terrain slope in degrees
  if(is.null(slope)==TRUE){
    slope <- raster::terrain(dtm, opt='slope', unit='degrees', neighbors=8)
  } else {
    slope <- slope
  }

  #define different types of cost functions and set the appropriate text to be used for subsequent plotting
  if (funct=="t") {
    #Tobler's hiking function; kmh; tan(x*pi/180) turns the slope from degrees to rise over run
    cost_function <- function(x){6 * exp(-3.5 * abs(tan(x*pi/180) + 0.05))}

    #set the labels to be used within the returned plot
    main.title <- paste0("Walking-time isochrones (in ", time, ") around origin")
    sub.title <- "Walking-time based on the Tobler's on-path hiking function"
    legend.cost <- paste0("walking-time (", time,")")
    sub.title.lcp.plot <- paste0("LCP(s) and walking-time distance(s) based on the Tobler's on-path hiking function (time in ", time, ") \nblack dot=start location\n red dot(s)=destination location(s)")
  }

  if (funct=="tofp") {
    #Tobler's hiking function off-path routes; kmh; tan(x*pi/180) turns the slope from degrees to rise over run
    #note that the multiplier 0.6 suggested by Tobler is meant to reduce the off-path walking speed
    cost_function <- function(x){(6 * exp(-3.5 * abs(tan(x*pi/180) + 0.05))) * 0.6}

    #set the labels to be used within the returned plot
    main.title <- paste0("Walking-time isochrones (in ", time, ") around origin")
    sub.title <- "Walking-time based on the Tobler's off-path hiking function"
    legend.cost <- paste0("walking-time (", time,")")
    sub.title.lcp.plot <- paste0("LCP(s) and walking-time distance(s) based on the Tobler's off-path hiking function (time in ", time, ") \nblack dot=start location\n red dot(s)=destination location(s)")
  }

  if(funct=="mt") {
    #Marquez-Perez et al.'s modified Tobler hiking function; kmh; tan(x*pi/180) turns the slope from degrees to rise over run
    cost_function <- function(x){4.8 * exp(-5.3 * abs((tan(x*pi/180) * 0.7) + 0.03))}

    #set the labels to be used within the returned plot
    main.title <- paste0("Walking-time isochrones (in ", time, ") around origin")
    sub.title <- "Walking-time based on the Marquez-Perez et al.'s modified Tobler hiking function"
    legend.cost <- paste0("walking-time (", time,")")
    sub.title.lcp.plot <- paste0("LCP(s) and walking-time distance(s) based on the Marquez-Perez et al.'s modified Tobler hiking function (time in ", time, ") \nblack dot=start location\n red dot(s)=destination location(s)")
  }

  if(funct=="ic") {
    #Irmischer-Clarke's modified Tobler hiking function; originally in m/s (males, on-path);
    # the formula is reshaped (multiplied by 3.6) below to turn it into kmh for consistency with the other Tobler-related cost functions;
    # also, originally the slope is in percent; tan(x*pi/180)*100 turns the slope from degrees to percent (=rise/run*100)
    cost_function <- function(x){(0.11 + exp(-(tan(x*pi/180)*100 + 5)^2 / (2 * 30)^2)) * 3.6}

    #set the labels to be used within the returned plot
    main.title <- paste0("Walking-time isochrones (in ", time, ") around origin")
    sub.title <- "Walking-time based on the (on-path) Irmischer-Clarke's modified Tobler hiking function"
    legend.cost <- paste0("walking-time (", time,")")
    sub.title.lcp.plot <- paste0("LCP(s) and walking-time distance(s) based on the (on-path) Irmischer-Clarke's modified Tobler hiking function (time in ", time, ") \nblack dot=start location\n red dot(s)=destination location(s)")
  }

  if(funct=="icofp") {
    #Irmischer-Clarke's modified Tobler hiking function; originally in m/s (males, off-path);
    # the formula is reshaped (multiplied by 3.6) below to turn it into kmh for consistency with the other Tobler-related cost functions;
    # also, originally the slope is in percent; tan(x*pi/180)*100 turns the slope from degrees to percent (=rise/run*100)
    cost_function <- function(x){(0.11 + 0.67 * exp(-(tan(x*pi/180)*100 + 2)^2 / (2 * 30)^2)) * 3.6}

    #set the labels to be used within the returned plot
    main.title <- paste0("Walking-time isochrones (in ", time, ") around origin")
    sub.title <- "Walking-time based on the (off-path) Irmischer-Clarke's modified Tobler hiking function"
    legend.cost <- paste0("walking-time (", time,")")
    sub.title.lcp.plot <- paste0("LCP(s) and walking-time distance(s) based on the (off-path) Irmischer-Clarke's modified Tobler hiking function (time in ", time, ") \nblack dot=start location\n red dot(s)=destination location(s)")
  }

  if(funct=="ug") {
    #Antonio Uriarte Gonzalez's slope-dependant walking-time cost function;
    # the cost function is originally expressed in seconds; here it is the reciprocal of time (1/T) that is used in order to eventually get
    #T/1. Also, originally the slope is in percent; tan(x*pi/180)*100 turns the slope from degrees to percent (=rise/run*100).
    #Note: unlike the original formula, here the pixel resolution is not taken into account since 'gdistance' takes care of the cells' dimension
    #when calculating accumulated costs.
    cost_function <- function(x){ 1/ (0.0277 * (tan(x*pi/180)*100) + 0.6115) }

    #set the labels to be used within the returned plot
    main.title <- paste0("Walking-time isochrones (in ", time, ") around origin")
    sub.title <- "Walking-time based on the Uriarte Gonzalez's cost function"
    legend.cost <- paste0("walking-time (", time,")")
    sub.title.lcp.plot <- paste0("LCP(s) and walking-time distance(s) based on the Uriarte Gonzalez's cost function (time in ", time, ") \nblack dot=start location\n red dot(s)=destination location(s)")
  }

  if(funct=="ree") {
    #relative energetic expenditure;
    # to calculate tangent of degrees (as requested by the cost function) we must first convert degrees to radians by multypling by pi/180
    cost_function <- function(x){ 1 / (tan(x*pi/180) / tan (1*pi/180)) }

    #set the labels to be used within the returned plot
    main.title <- "Accumulated cost isolines around origin"
    sub.title <- "Cost based on the slope-based relative energetic expenditure cost function"
    legend.cost <- "cost"
    sub.title.lcp.plot <- paste0("LCP(s) and cost distance(s) based on the slope-based relative energetic expenditure cost function \nblack dot=start location\n red dot(s)=destination location(s)")
  }

  if(funct=="hrz") {
    #Herzog metabolic cost function in J/(kg*m);
    #tan(x*pi/180) turns slope from degrees to rise/run, which is requested by the cost function
    cost_function <- function(x){ 1 / ((1337.8 * tan(x*pi/180)^6 + 278.19 * tan(x*pi/180)^5 - 517.39 * tan(x*pi/180)^4 - 78.199 * tan(x*pi/180)^3 + 93.419 * tan(x*pi/180)^2 + 19.825 * tan(x*pi/180) + 1.64)) }

    #set the labels to be used within the returned plot
    main.title <- "Accumulated cost isolines around origin"
    sub.title <- "Cost based on the Herzog's metabolic cost function \n cost in J / (Kg*m)"
    legend.cost <- "metabolic cost J / (Kg*m)"
    sub.title.lcp.plot <- paste0("LCP(s) and cost distance(s) based on the Herzog's metabolic cost function \ncost in J / (Kg*m) \nblack dot=start location\n red dot(s)=destination location(s)")
  }

  if(funct=="wcs") {
    #wheel critical slope cost function; tan(x*pi/180)*100 turns the slope from degrees to percent; the latter is requested by the cost function
    cost_function <- function(x){ 1 / (1 + ((tan(x*pi/180)*100) / sl.crit)^2) }

    #set the labels to be used within the returned plot
    main.title <- "Accumulated cost isolines around origin"
    sub.title <- paste0("Cost based on the wheeled-vehicle critical slope cost function \ncritical slope set to ", sl.crit, " percent")
    legend.cost <- "cost"
    sub.title.lcp.plot <- paste0("LCP(s) and cost distance(s) based on the wheeled-vehicle critical slope cost function \ncritical slope set to ", sl.crit, " percent \nblack dot=start location\n red dot(s)=destination location(s)")
  }

  if(funct=="vl") {
    #Van Leusen's metabolic energy expenditure cost function
    #note: V is velocity in m/s; tan(x*pi/180)*100 turns the slope from degrees to percent
    cost_function <- function(x){ 1 / (1.5 * W + 2.0 * (W + L) * (L / W)^2 + N * (W + L) * (1.5 * V^2 + 0.35 * V * (tan(x*pi/180)*100) + 10)) }

    #set the labels to be used within the returned plot
    main.title <- "Accumulated cost isolines around origin"
    sub.title <- paste0("Cost based on the Van Leusen's metabolic energy expenditure cost function \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V: ", V)
    legend.cost <- "energy expenditure cost (Megawatts)"
    sub.title.lcp.plot <- paste0("LCP(s) and cost distance(s) based on the Van Leusen's metabolic energy expenditure cost function \n cost in Megawatts; parameters: W: ", W, "; L: ", L, "; N: ", N, "; V: ", V, "\nblack dot=start location\n red dot(s)=destination location(s)")
  }

  if(funct=="p") {
    #Pandolf et al.'s metabolic energy expenditure cost function
    #note: V is velocity in m/s; tan(x*pi/180)*100 turns the slope from degrees to percent
    cost_function <- function(x){ 1 / (1.5 * W + 2.0 * (W + L) * (L / W)^2 + N * (W + L) * (1.5 * V^2 + 0.35 * V * (tan(x*pi/180)*100))) }

    #set the labels to be used within the returned plot
    main.title <- "Accumulated cost isolines around origin"
    sub.title <- paste0("Cost based on the Pandolf et al.'s metabolic energy expenditure cost function \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V: ", V)
    legend.cost <- "energy expenditure cost (Megawatts)"
    sub.title.lcp.plot <- paste0("LCP(s) and cost distance(s) based on the Pandolf et al.'s metabolic energy expenditure cost function \n cost in Megawatts; parameters: W: ", W, "; L: ", L, "; N: ", N, "; V: ", V, "\nblack dot=start location\n red dot(s)=destination location(s)")
  }

  #cost calculation for walking-speed-based cost functions
  if (funct=="t" | funct=="tofp" | funct=="mt" | funct=="ic" | funct=="icofp") {
    #calculate the walking speed for each cell of the slope raster
    speed.kmh <- raster::calc(slope, cost_function)

    #turn the walking speed from kmh to ms (0.278=1000/3600)
    speed.ms <- speed.kmh * 0.278

    #create a transitional matrix for the speed in ms
    cost.TRANS <- gdistance::transition(speed.ms, transitionFunction=mean, directions=moves)
  }

  #cost calculation for other types of cost functions;
  #note the Uriarte Gonzalez's slope-dependant walking-time cost function is in this group since (unlike the above functions)
  #it expresses cost as time NOT speed
  if (funct=="ree" | funct=="hrz" | funct=="wcs" | funct=="vl" | funct=="p" | funct=="ug") {
    #calculate the cost for each cell of the slope raster
    cost <- raster::calc(slope, cost_function)

    #create a transitional matrix for the cost
    cost.TRANS <- gdistance::transition(cost, transitionFunction=mean, directions=moves)
  }

  #geocorrection of the cost matrix
  Conductance <- gdistance::geoCorrection(cost.TRANS)

  #accumulate the pace outwards from the origin
  accum_final <- gdistance::accCost(Conductance, sp::coordinates(origin))

  #if user select the Tobler's, the modified Tobler's, or the Uriarte Gonzalez's function, turn seconds into the user-defined time-scale
  if (funct=="t" | funct=="tofp" | funct=="mt" | funct=="ic" | funct=="ug") {
    if (time=="h") {
      #turn seconds into hours
      accum_final <- accum_final / 3600
    } else {
      #turn seconds into minutes
      accum_final <- accum_final / 60
    }
  }

  #if user select the Val Leusen's or the Pandolf et al.'s function, turn the cost from Watts to Megawatts
  if (funct=="vl" | funct=="p") {
    accum_final <- accum_final / 1000000
  }

  #if no break value is entered, set the breaks to one tenth of the range of the values of the final accumulated cost surface
  if(is.null(breaks)==TRUE){
    #exclude the inf values from the calculation
    breaks <- round((max(accum_final[][is.finite(accum_final[])]) - min(accum_final[][is.finite(accum_final[])])) / 10,2)
  }

  #set the break values for the isolines, again excluding inf values
  levels <- seq(min(accum_final[][is.finite(accum_final[])]), max(accum_final[][is.finite(accum_final[])]), breaks)

  #conditionally set the layout in just one visualization
  if(is.null(destin)==FALSE & oneplot==TRUE){
    m <- rbind(c(1,2))
    layout(m)
  }

  #produce the output
  if (outp=="r") {
    #produce a raster with contours
    raster::plot(accum_final,
         main=main.title,
         sub=sub.title,
         cex.main=0.95,
         cex.sub=0.75,
         legend.lab=legend.cost,
         col = terrain.colors(255))
    raster::contour(accum_final,
            add=TRUE,
            levels=levels,
            labcex=cex.breaks,
            drawlabels = cont.lab)
    raster:: plot(origin,
         pch=20,
         add=TRUE)

  } else {
    #only produce contours
    raster::contour(accum_final,
            levels=levels,
            main=main.title,
            sub=sub.title,
            cex.main=0.95,
            cex.sub=0.75,
            labcex=cex.breaks,
            drawlabels = cont.lab)
    raster::plot(origin,
         pch=20,
         add=TRUE)
  }

  #calculate and store the contours as a SpatialLinesDataFrame
  isolines <- raster::rasterToContour(accum_final, levels=levels)

  #if 'destin' is NOT NULL, calculate the least-cost path(s) from the origin to the destination(s);
  #the 'Conductance' transitional layer is used
  if(is.null(destin)==FALSE){
    #calculate the least-cost path(s)
    sPath <- gdistance::shortestPath(Conductance, sp::coordinates(origin), sp::coordinates(destin), output="SpatialLines")

    #plot the dtm
    raster::plot(dtm, main="Digital Terrain Model with Least-cost Path(s)",
         sub=sub.title.lcp.plot,
         cex.main=0.90,
         cex.sub=0.7,
         legend.lab="Elevation (masl)")

    #add the origin
    raster::plot(origin, add=TRUE, pch=20)

    #add the destination(s)
    raster::plot(destin,
         add=TRUE,
         pch=20,
         col="red")

    #add the LCPs
    graphics::lines(sPath)

    #calculate the length of the least-cost paths and store the values by appending them to a new variable of the sPath object
    sPath$length <- rgeos::gLength(sPath, byid=TRUE)

    #extract the cost from the accum_final to the destination location(s), appending the data to a new column
    destin$cost <- raster::extract(accum_final, destin)

    #if destin.lab is TRUE, add the point(s)'s labels
    if(destin.lab==TRUE){
      raster::text(sp::coordinates(destin),
           labels=round(destin$cost,2),
           pos = 4,
           cex=cex.lcp.lab)

      #if export is TRUE, export the LPCs as a shapefile
      if(export==TRUE){
        rgdal::writeOGR(sPath, ".", paste0("LCPs_", funct), driver="ESRI Shapefile")
      }
    }

  } else {
    sPath=NULL
    dest.loc.w.cost=NULL
  }

  #if export is TRUE, export the accumulated cost surface as a raster file
  #and the isolines as a shapefile
  if(export==TRUE){
    raster::writeRaster(accum_final, paste0("accum_cost_surf_", funct), format="GTiff")
    rgdal::writeOGR(isolines, ".", paste0("isolines_", funct), driver="ESRI Shapefile")
  }

  #restore the original graphical device's settings if previously modified
  if(is.null(destin)==FALSE & oneplot==TRUE){
    par(mfrow = c(1,1))
  }

  results <- list("accumulated.cost.raster"=accum_final,
                  "isolines" = isolines,
                  "LCPs"=sPath,
                  "dest.loc.w.cost"=destin)
}
