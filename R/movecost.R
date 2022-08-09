#' R function for calculating accumulated anisotropic slope-dependant cost of movement across the terrain and least-cost paths from a point origin
#'
#' The function provides the facility to calculate the anisotropic accumulated cost of movement around a starting location and to optionally calculate least-cost path(s) toward
#' one or multiple destinations. It implements different cost estimations related to human movement across the landscape.
#' The function takes as input a Digital Terrain Model ('RasterLayer' class) and a point feature ('SpatialPointsDataFrame' class), the latter representing
#' the starting location, i.e. the location from which the accumulated cost is calculated. Besides citing this package,
#' you may want to refer to the following journal article: \strong{Alberti (2019) <doi:10.1016/j.softx.2019.100331>}.\cr
#' Visit this \href{https://drive.google.com/file/d/1gLDrkZFh1b_glzCEqKdkPrer72JJ9Ffa/view?usp=sharing}{LINK} to access the package's vignette.\cr
#'
#' If the parameter 'destin' is fed with a dataset representing destination location(s) ('SpatialPointsDataFrame' class), the function also calculates
#' least-cost path(s) plotted on the input DTM; the length of each path is saved under the variable 'length' stored in the 'LCPs' dataset ('SpatialLines' class) returned by the function.
#' In the produced plot, the red dot(s) representing the destination location(s) are labelled with numeric values representing
#' the cost value at the location(s). \cr
#'
#' The cost value is also appended to the updated destination locations dataset returned by the function, which
#' stores a new variable labelled 'cost'. If the cost is expressed in terms of walking time, the labels accompaining each destinaton location will
#' express time in sexagesimal numbers (hours, minutes, seconds). In this case, the variable 'cost' appended to the returned destination location datset
#' will store the time figures in decimal numbers, while another variable named 'cost_hms' will store the corresponding value in sexagesimal numbers.
#' When interpreting the time values stored in the 'cost' variable, the user may want to bear in mind the selected time unit (see right below).\cr
#'
#' When using cost functions expressing cost in terms of time, the time unit can be selected by the user setting the 'time' parameter to 'h' (hours) or to 'm' (minutes).\cr
#'
#' In general, the user can also select which type of visualization the function has to produce; this is achieved setting the 'outp' parameter to either 'r' (=raster)
#' or to 'c' (=contours). The former will produce a raster with a colour scale and contour lines representing the accumulated cost surface; the latter parameter will only
#' produce contour lines.\cr
#'
#' The contour lines' interval is set using the parameter 'breaks'; if no value is passed to the parameter, the interval will be set by default to
#' 1/10 of the range of values of the accumulated cost surface.\cr
#'
#' It is worth reminding the user(s) that all the input layers (i.e., DTM, start location, and destination locations) must use the same projected coordinate system.\cr
#'
#'
#' \strong{Cost surface calculation}:\cr
#' for the cost-surface and LCPs calculation, 'movecost' builds on functions from Jacob van Etten's
#' \href{https://cran.r-project.org/package=gdistance}{gdistance} package.
#' Under the hood, 'movecost' calculates the slope as rise over run, following the procedure described
#' by van Etten, "R Package gdistance: Distances and Routes on Geographical Grids" in Journal of Statistical Software 76(13), 2017, pp. 14-15.
#' The number of directions in which cells are connected in the cost calculation can be set to 4 (rook's case), 8 (queen's case), or
#' 16 (knight and one-cell queen moves) using the 'move' parameter (see 'Arguments').\cr
#'
#'
#' \strong{Inibition of movement (barrier)}:\cr
#' areas where the movement is inhibited can be fed into the analysis via the \code{barrier} parameter; SpatialLineDataFrame or SpatialPolygonDataFrame can be used.
#' The barrier is assigned a conductance value of 0 (i.e., movement is inhibited) by default, but the user can assign any other value via the
#' \code{field} parameter. Internally, the barrier creation rests on the \code{\link[leastcostpath]{create_barrier_cs}} function
#' from the \emph{leastcostpath} package.\cr
#'
#' To test this facility, consider the following example:\cr
#'
#' First, we use in-built data to come up with a linear feature (i.e., a LCP) that we will later use as barrier:\cr
#'
#' result1 <- movecost(volc, destin.loc[1,], destin.loc[4,])\cr
#'
#' After, we calculate the LCP between two other locations, first not using any barrier (result2), then using the mentioned LCP (from result1)
#' as a barrier (result3):\cr
#'
#' result2 <- movecost(volc, destin.loc[3,], destin.loc[6,], move=8)\cr
#' result3 <- movecost(volc, destin.loc[3,], destin.loc[6,], barrier=result1$LCPs, move=8)\cr
#'
#' As apparent by comparing result2 to result3, when the barrier is used (result3), the LCP does not cross the barrier but is "forced"
#' to make a long detour. \strong{Note} that the \code{move} parameter has been set to 8; if set to 16, the LCP will be "able" to jump the barrier.\cr
#'
#'
#' \strong{Acquiring online elevation data}:\cr
#' if a DTM is not provided, 'movecost()' will download elevation data from online sources.
#' Elevation data will be aquired for the area enclosed  by the  polygon supplied by the 'studyplot' parameter (SpatialPolygonDataFrame class).
#' To tap online elevation data, 'movecost' internally builds on the
#' \code{\link[elevatr]{get_elev_raster}} function from the \emph{elevatr} package.\cr
#'
#' The zoom level of the downloaded DTM (i.e., its resolution) is controlled by the parameter 'z', which is
#' set to 9 by default (a trade off between resolution and download time).\cr
#'
#' To test this facility, the user may want to try the following code, that will generate a least-cost surface and least-cost paths
#' in an area close the Mount Etna (Sicily, Italy), whose elevation data are acquired online; the start and end locations, and the
#' polygon defining the study area, are provided in this same package:\cr
#'
#' result <- movecost(origin=Etna_start_location, destin=Etna_end_location, studyplot=Etna_boundary) \cr
#'
#' The LCPs back to the origin can be calculated and plotted setting the parameter 'return.base' to TRUE:\cr
#'
#' result <- movecost(origin=Etna_start_location, destin=Etna_end_location, studyplot=Etna_boundary, return.base=TRUE) \cr
#'
#' To know more about what elevation data are tapped from online
#' sources, visit: https://cran.r-project.org/web/packages/elevatr/vignettes/introduction_to_elevatr.html. \cr
#'
#' For more information about the elevation data resolution per zoom level, visit
#' https://github.com/tilezen/joerd/blob/master/docs/data-sources.md#what-is-the-ground-resolution.\cr
#'
#' To know what is sourced at what zoom level, visit
#' https://github.com/tilezen/joerd/blob/master/docs/data-sources.md#what-is-sourced-at-what-zooms. \cr
#'
#'
#' \strong{Terrain slope and cognitive slope}:\cr
#' when it comes to the terrain slope, the function provides the facility to use the so-called 'cognitive slope',
#' following Pingel TJ (2013), Modeling Slope as a Contributor to Route Selection in Mountainous Areas, in Cartography and Geographic Information Science, 37(2), 137-148.
#' According to Pingel, "Humans tend to overestimate geographic slopes by a surprisingly high margin...This analysis indicates downhill slopes are overestimated
#' at approximately 2.3 times the vertical, while uphill slopes are overestimated at 2 times the vertical.". As a result,
#' if the parameter 'cogn.slp' is set to TRUE, positive slope values are preliminarily multiplied by 1.99, while negative slope values are multiplied by 2.31.
#'
#'
#' \strong{Terrain factor (N)}:\cr
#' virtually all the implemented cost functions (with few exceptions) can take into account a 'terrain factor' (N parameter; 1 by default), which
#' represents the easiness/difficulty of moving on different terrain types. According to the type of terrain, the energy spent when walking
#' increases. The same holds true for time, which increases because on a loose terrain (for instance) the walking speed decreases.
#' While a terrain factor is 'natively' part of the Van Leusen's, Pandolf et al.'s, and Bellavia's cost function,
#' it has been integrated into the other cost functions as well (when/if relevant).\cr
#'
#' Note that the terrain factor has NOT been implemented in the Alberti's, Tobler's off-path, and Irmischer-Clarke's off-path cost function.
#' As for the latter two, they already natively feature a terrain factor. Therefore, it has been implemented only in their on-path version.
#' Needless to say, if we use a terrain factor of 1.67 with the Tobler's (on-path) hiking function, the results
#' will be equal to those obtained using the Tobler's off-path function (the reciprocal of 1.67, i.e. 0.60, is in fact
#' natively used by the Tobler's function for off-path hiking). In fact, compare the results of the following two runs
#' of 'movecost()' (first some data are loaded):\cr
#'
#' volc <- raster::raster(system.file("external/maungawhau.grd", package="gdistance"))\cr
#' data(volc.loc)\cr
#' data(destin.loc)\cr
#'
#' result1 <- movecost(dtm=volc, origin=volc.loc, destin=destin.loc, breaks=0.05, funct="t", N=1.67)\cr
#' result2 <- movecost(dtm=volc, origin=volc.loc, destin=destin.loc, breaks=0.05, funct="tofp")\cr
#'
#' The user may want to refer to the following \strong{list of terrain factors}, which is based on the data collected in Herzog, I. (2020).
#' Spatial Analysis Based on Cost Functions. In Gillings M, Haciguzeller P, Lock G (eds), "Archaeological Spatial Analysis. A Methodological Guide.",
#' Routledge: New York, 340 (with previous references). The list is divided into two sections (a and b), the first reporting the terrain
#' factors to be used for cost functions measuring time, the second for functions measuring cost other than time:\cr
#'
#' (a)\cr
#' \itemize{
#'   \item Blacktop roads, improved dirt paths, cement = 1.00
#'   \item Lawn grass = 1.03
#'   \item Loose beach sand = 1.19
#'   \item Disturbed ground (former stone quarry) = 1.24
#'   \item Horse riding path, flat trails and meadows = 1.25
#'   \item Tall grassland (with thistle and nettles) = 1.35
#'   \item Open space above the treeline (i.e., 2000 m asl) = 1.50
#'   \item Bad trails, stony outcrops and river beds = 1.67
#'   \item Off-paths = 1.67
#'   \item Bog = 1.79
#'   \item Off-path areas below the treeline (pastures, forests, heathland) = 2.00
#'   \item Rock = 2.50
#'   \item Swamp, water course = 5.00
#' }
#'
#' (b)\cr
#' \itemize{
#'   \item Asphalt/blacktop = 1.00
#'   \item Dirt road or grass = 1.10
#'   \item Hard-surface road = 1.20
#'   \item Light brush = 1.20
#'   \item Ploughed field = 1.30 or 1.50
#'   \item Heavy brush = 1.50
#'   \item Hard-packed snow = 1.60
#'   \item Swampy bog = 1.80
#'   \item Sand dunes = 1.80
#'   \item Loose sand = 2.10
#' }
#'
#'
#'
#' \strong{Implemented cost functions}:\cr
#' note that in what follows \strong{x[adj]} stands for slope as rise/run calculated for adjacent cells:\cr
#'
#' \strong{Tobler's hiking function (on-path) (speed in kmh)}:\cr
#'
#' \eqn{ (6 * exp(-3.5 * abs(x[adj] + 0.05))) * (1/N) }\cr
#'
#'
#' \strong{Tobler's hiking function (off-path) (speed in kmh)}:\cr
#'
#' \eqn{(6 * exp(-3.5 * abs(x[adj] + 0.05))) * 0.6}\cr
#'
#' as per Tobler's indication, the off-path walking speed is reduced by 0.6.\cr
#'
#'
#' \strong{Marquez-Perez et al.'s modified Tobler hiking function (speed in kmh)}:\cr
#'
#' \eqn{ (4.8 * exp(-5.3 * abs((x[adj] * 0.7) + 0.03))) * (1/N) }\cr
#'
#' modified version of the Tobler's hiking function as proposed by Joaquin Marquez-Parez, Ismael Vallejo-Villalta & Jose I. Alvarez-Francoso (2017), "Estimated travel time for walking trails in natural areas",
#' Geografisk Tidsskrift-Danish Journal of Geography, 117:1, 53-62, DOI: 10.1080/00167223.2017.1316212.\cr
#'
#'
#' \strong{Irmischer-Clarke's modified Tobler hiking function (male, on-path; speed in kmh)}:\cr
#'
#' \eqn{ ((0.11 + exp(-(abs(x[adj])*100 + 5)^2 / (2 * 30^2))) * 3.6) * (1/N) }\cr
#'
#' modified version of the Tobler's function as proposed for (male) on-path hiking by Irmischer, I. J., & Clarke, K. C. (2018). Measuring and modeling the speed of human navigation.
#' Cartography and Geographic Information Science, 45(2), 177-186. https://doi.org/10.1080/15230406.2017.1292150.
#' It is interesting to note that the hiking speed predicted by this and by the other functions proposed by the authors is slower than the one
#' modelled by Tobler's hiking function. This is attributed to the cognition involved in wayfinding
#', such as map reading, analyzing the terrain, decision making, determining routes, etc.
#' \strong{Note}: all the all the Irmischer-Clarke's functions originally express speed in m/s; they have been reshaped (multiplied by 3.6) to turn m/s into km/h for consistency
#' with the other Tobler-related cost functions; slope is in percent.\cr
#'
#'\strong{Irmischer-Clarke's modified Tobler hiking function (male, off-path; speed in kmh)}:\cr
#'
#' \eqn{(0.11 + 0.67 * exp(-(abs(x[adj])*100 + 2)^2 / (2 * 30)^2)) * 3.6}\cr
#'
#' \strong{Irmischer-Clarke's modified Tobler hiking function (female, on-path; speed in kmh)}:\cr
#'
#' \eqn{ ((0.95 * (0.11 + exp(-(abs(x[adj]) * 100 + 5)^2/(2 * 30^2)))) * 3.6) * (1/N) }\cr
#'
#' \strong{Irmischer-Clarke's modified Tobler hiking function (female, off-path; speed in kmh)}:\cr
#'
#' \eqn{(0.95 * (0.11 + 0.67 * exp(-(abs(x[adj]) * 100 + 2)^2/(2 * 30^2)))) * 3.6}\cr
#'
#'
#'\strong{Uriarte Gonzalez's walking-time cost function}:\cr
#'
#' \eqn{ 1 / ((0.0277 * (abs(x[adj])*100) + 0.6115) * N) }\cr
#'
#' proposed by Uriarte Gonzalez;
#' \strong{see}: Chapa Brunet, T., Garcia, J., Mayoral Herrera, V., & Uriarte Gonzalez, A. (2008). GIS landscape models for the study of preindustrial settlement patterns in Mediterranean areas.
#' In Geoinformation Technologies for Geo-Cultural Landscapes (pp. 255-273). CRC Press. https://doi.org/10.1201/9780203881613.ch12.\cr
#' The cost function originally expresses walking time in seconds; for the purpose of its implementation in this function, it is the reciprocal of time (1/T) that is used in order to eventually get
#' T/1. Unlike in the original cost function, here the pixel resolution is not taken into account since 'gdistance' takes care of the cells' dimension
#' when calculating accumulated costs.\cr
#'
#'
#'
#' \strong{Marin Arroyo's walking-time cost function}:\cr
#'
#' \eqn{ ifelse((abs(x[adj])*100) < 0, 1 / ((0.6 * ((abs(x[adj])*100)/23+1))*N), 1 / ((0.6 * ((abs(x[adj])*100)/11+1))*N)) }\cr
#'
#' used by Marin Arroyo A.B. (2009), The use of optimal foraging theory to estimate Late Glacial site catchments areas from a central place:
#' the case of eastern Cantabria, Span, in Journal of Anthropological Archaeology 28, 27-36.
#' The cost function originally expresses walking time in seconds; here it is the reciprocal of time (1/T) that is used in order to eventually get
#' T/1. Slope is in percent. Note: unlike in the original equation, here d (distance travelled in meter) is not taken into account since 'gdistance'
#' takes care of the cells' dimension when calculating accumulated costs.\cr
#'
#'
#'
#'\strong{Alberti's Tobler hiking function modified for pastoral foraging excursions (speed in kmh)}:\cr
#'
#' \eqn{(6 * exp(-3.5 * abs(x[adj] + 0.05))) * 0.25}\cr
#'
#' proposed by Gianmarco Alberti;
#' \strong{see}: \href{https://www.um.edu.mt/library/oar/bitstream/123456789/64172/1/Chapter_9_Locating_potential_pastoral_foraging_routes.pdf}{Locating potential pastoral foraging routes in Malta through the use of a Geographic Information System}.
#' The Tobler's function has been rescaled to fit animal walking speed during foraging excursions. The distribution of the latter, as empirical data show, turns out to be right-skewed
#' and to vary along a continuum. It ranges from very low speed values (corresponding to slow grazing activities grazing while walking) to comparatively higher values
#' (up to about 4.0 km/h) corresponding to travels without grazing (directional travel toward feeding stations).
#' The function consider 1.5 km/h as the average flock speed, which roughly corresponds to the average speed recorded in some studies, and
#' can be considered the typical speed of flocks during excursions in which grazing takes place while walking  (typical form of grazing in most situations).
#' Tobler's hiking function has been rescaled by a factor of 0.25 to represent the walking pace of a flock instead of humans.
#' The factor corresponds to the ratio between the flock average speed (1.5 km/h) and the maximum human walking speed (about 6.0 km/h) on a favourable slope.\cr
#'
#'
#'
#' \strong{Garmy, Kaddouri, Rozenblat, and Schneider's hiking function (speed in kmh)}:\cr
#'
#' \eqn{ (4 * exp(-0.008 * ((atan(abs(x[adj]))*180/pi)^2))) * (1/N) }\cr
#'
#' slope in degrees;
#' \strong{see}: Herzog, I. (2020). Spatial Analysis Based on Cost Functions. In Gillings M, Haciguzeller P, Lock G (eds), "Archaeological Spatial Analysis. A Methodological Guide.",
#' Routledge: New York, 333-358 (with previous references).\cr
#'
#'
#'
#' \strong{Rees' hiking function (speed in kmh)}:\cr
#'
#' \eqn{ ((1 / (0.75 + 0.09 * abs(x[adj]) + 14.6 * (abs(x[adj]))^2)) * 3.6) * (1/N) }\cr
#'
#' Rees' slope-dependant cost function; it is originally expressed in terms of time (1/v in Rees' publication);
#' here it is the reciprocal of time (i.e. speed) that is used in order to eventually get the reciprocal of speed (i.e. time).
#' Slope is dealt with here as originally expressed in Rees' publication (i.e. rise over run). The speed, which is originally expressed in m/s,
#' has been here transposed to kmh (i.e., multiplied by 3.6) for consistency with other hiking functions.\cr
#' For this cost function \strong{see}: Rees, WG (2004). Least-cost paths in mountainous terrain.
#' Computers & Geosciences, 30(3), 203-209. See also: Campbell MJ, Dennison PE, Butler BW, Page WG (2019). Using crowdsourced
#' fitness tracker data to model the relationship between slope and travel rates. Applied Geography 106, 93-107 (with previous references).\cr
#'
#'
#'
#'\strong{Kondo-Seino's modified Tobler hiking function (speed in kmh)}:\cr
#'
#' \eqn{ ifelse(abs(x[adj]) >= -0.07, (5.1 * exp(-2.25 * abs(x[adj] + 0.07))) * (1/N), (5.1 * exp(-1.5 * abs(x[adj] + 0.07)))) * (1/N) }\cr
#'
#' Kondo-Seino's modified Tobler hiking function; it expresses walking speed in Kmh; slope as rise/run;
#' \strong{see} Kondo Y., Seino Y. (2010). GPS-aided Walking Experiments and Data-driven Travel Cost Modelingon the Historical Road of Nakasendo-Kisoji
#' (Central Highland Japan), in: Frischer B., Webb Crawford J., Koller D. (eds.), Making History Interactive.
#' Computer Applications and Quantitative Methods in Archaeology (CAA). Proceedings of the 37th International Conference, Williamsburg, Virginia, United States of America,
#' March 22-26 (BAR International Series S2079). Archaeopress, Oxford, 158-165.\cr
#'
#'
#'
#' \strong{Wheeled-vehicle critical slope cost function}:\cr
#'
#' \eqn{ 1 / ((1 + ((abs(x[adj])*100) / sl.crit)^2) * N) }\cr
#'
#' where \eqn{sl.crit} (=critical slope, in percent) is "the transition where switchbacks become more effective than direct uphill or downhill paths" and typically is in the range 8-16;
#' \strong{see} Herzog, I. (2016). Potential and Limits of Optimal Path Analysis. In A. Bevan & M. Lake (Eds.), Computational Approaches to Archaeological Spaces (pp. 179-211). New York: Routledge. \cr
#'
#'
#'
#'\strong{Relative energetic expenditure cost function}:\cr
#'
#' \eqn{ 1 / ((tan((atan(abs(x[adj]))*180/pi)*pi/180) / tan (1*pi/180)) * N) }\cr
#'
#' slope-based cost function expressing change in potential energy expenditure;
#' \strong{see} Conolly, J., & Lake, M. (2006). Geographic Information Systems in Archaeology. Cambridge: Cambridge University Press, p. 220;
#' \strong{see also} Newhard, J. M. L., Levine, N. S., & Phebus, A. D. (2014). The development of integrated terrestrial and marine pathways in the Argo-Saronic region, Greece. Cartography and Geographic Information Science, 41(4), 379-390, with references to studies that use this
#' function; \strong{see also} ten Bruggencate, R. E., Stup, J. P., Milne, S. B., Stenton, D. R., Park, R. W., & Fayek, M. (2016). A human-centered GIS approach to modeling mobility on southern Baffin Island, Nunavut,
#' Canada. Journal of Field Archaeology, 41(6), 684-698. https://doi.org/10.1080/00934690.2016.1234897.\cr
#'
#'
#'
#' \strong{Bellavia's cost function}:\cr
#'
#' \eqn{1 / (N * ((atan(abs(x[adj]))*180/pi)+1))}\cr
#'
#' proposed by G. Bellavia, it measures abstract cost. Slope in degrees; N is a terrain factor (see above).
#' \strong{See}: Herzog I. (2020). Spatial Analysis Based on Cost Functions. In Gillings M, Haciguzeller P, Lock G (eds), "Archaeological Spatial Analysis. A Methodological Guide.",
#' Routledge: New York, 333-358 (with previous references).\cr
#'
#'
#'
#' \strong{Eastman's cost function}:\cr
#'
#' \eqn{ 1 / ((0.031*(atan(abs(x[adj]))*180/pi)^2-0.025*(atan(abs(x[adj]))*180/pi)+1)*N) }\cr
#'
#' proposed by J.R. Eastman, it measures abstract cost; slope in degrees.
#' \strong{See}: Vaissie E., Mobility of Paleolithic Populations: Biomechanical Considerations and Spatiotemporal Modelling,
#' in PaleoAnthropology 2021 (1): 120-144 (with previous reference to Eastman 1999).\cr
#'
#'
#'
#' \strong{Pandolf et al.'s metabolic energy expenditure cost function (in Watts)}:\cr
#'
#' \eqn{ 1 / ((1.5 * W + 2.0 * (W + L) * (L / W)^2 + N * (W + L) * (1.5 * (V^2) + 0.35 * V * (abs(x[adj])*100)))*N) }\cr
#'
#' where \eqn{W} is the walker's body weight (Kg), \eqn{L} is the carried load (in Kg), \eqn{V} is the velocity in m/s, \eqn{N} is a coefficient representing ease of movement on the terrain (see above).
#' \strong{Note} that if \eqn{V} is set to 0 by the user, it is internally worked out on the basis of the Tobler function for on-path hiking; therefore, \eqn{V} will not be
#' considered constant throughout the analysed area, but will vary as function of the slope.
#'
#' For this cost function, \strong{see} Pandolf, K. B., Givoni, B., & Goldman, R. F. (1977). Predicting energy expenditure with loads while standing or walking very slowly. Journal of Applied Physiology,
#' 43(4), 577-581. https://doi.org/10.1152/jappl.1977.43.4.577.\cr
#'
#' For the use of this cost function in a case study, \strong{see} Rademaker, K., Reid, D. A., & Bromley, G. R. M. (2012). Connecting the Dots: Least Cost Analysis, Paleogeography, and
#' the Search for Paleoindian Sites in Southern Highland Peru. In White D.A. & Surface-Evans S.L. (Eds.), Least Cost Analysis of Social Landscapes. Archaeological Case Studies (pp. 32-45).
#' University of Utah Press; \strong{see also} Herzog, I. (2013). Least-cost Paths - Some Methodological Issues, Internet Archaeology 36 (http://intarch.ac.uk/journal/issue36/index.html) with references.
#' For the idea of using an hiking function inside an energetic cost function, \strong{see} for instance White D.A., Prehistoric Trail Networks of the Western Papaguarie.
#' A Multifaceted Least Cost Graph Theory Analysis. In White D.A. & Surface-Evans S.L. (Eds.), Least Cost Analysis of Social Landscapes. Archaeological Case Studies (pp. 188-206).
#' University of Utah Press.\cr
#'
#' \strong{Note}: in the returned charts, the cost is transposed from Watts to Megawatts (see, e.g., Rademaker et al 2012 cited above).\cr
#'
#'
#'
#' \strong{Pandolf et al.'s metabolic energy expenditure cost function with correction factor for downhill movements (in Watts)}:\cr
#'
#' \eqn{ ifelse(abs(x[adj])*100 > 0 , 1 / (1.5 * W + 2.0 * (W+L) * (L/W)^2 + N * (W+L) * (1.5 * V^2 + 0.35 * V * (abs(x[adj])*100))),
#' 1 / ((1.5 * W + 2.0 * (W+L) * (L/W)^2 + N * (W+L) * (1.5 * V^2 + 0.35 * V * (abs(x[adj])*100))) - (N * ((abs(x[adj])*100) * (W+L) * V/3.5) - ((W+L) * ((abs(x[adj])*100)+6)^2/W) + (25-V^2)))) }\cr
#'
#' for the parameters \eqn{W}, \eqn{L}, \eqn{V}, and \eqn{N}, see above. If \eqn{V} is set to 0 by the user, it is internally worked out on the basis of the Tobler function for on-path hiking; therefore, \eqn{V} will not be
#' considered constant throughout the analysed area, but will vary as function of the slope.
#' For the correction factor applied to the Pandolf et al.'s cost function, \strong{see} Yokota M., Berglund L.G., Santee W.R., Buller M.J., Hoyt R.W. (2004),
#' Predicting Individual Physiological Responsed During Marksmanship Field Training Using an Updates Scenario-J Model. U.S. Army Research Institute of Environmental Medicine
#' Technical Report T04-09. For an archaeological application of the Pandol et al's cost function with correction factor, \strong{see} White D.A., Prehistoric Trail Networks of the Western Papaguarie.
#' A Multifaceted Least Cost Graph Theory Analysis. In White D.A. & Surface-Evans S.L. (Eds.), Least Cost Analysis of Social Landscapes. Archaeological Case Studies (pp. 188-206).
#' University of Utah Press.\cr
#'
#' \strong{Note}: in the returned charts, the cost is transposed from Watts to Megawatts (see, e.g., Rademaker et al 2012 cited above).\cr
#'
#'
#'
#' \strong{Minetti et al.'s metabolic energy cost function (in J/(kg*m))}:\cr
#'
#' \eqn{ 1 / (((280.5 * abs(x[adj])^5) - (58.7 * abs(x[adj])^4) - (76.8 * abs(x[adj])^3) + (51.9 * abs(x[adj])^2) + (19.6 * abs(x[adj])) + 2.5) * N) }\cr
#'
#' \strong{see} Minetti A.E., Moia C., Roi G.S., Susta D., Ferretti G. (2002), Enery cost of walking and running at extreme uphilland downhill slopes,
#' in Journal of Applied Physiology 93, 1039-1046. \strong{Note} that this equation is valid for slopes in the range -0.5/0.5; outside this range,
#' the function's output becomes counterintuitive, as noted in Paez et al. in Journal of Transport Geography 82 (2020) and Herzog I. (2013), Theory and practice
#' of cost functions, in Contreras F., Farjas M., Melero F.J. (eds), "Fusion of cultures. Proceedings of the 38th annual conference on computer
#' applications and quantitative methods in archaeology". BAR IS, 2494, 375-382. Oxford: Archaeopress. In the latter work, Herzog proposes to
#' replace Minetti et al.'s equation with its 6th degrees polynomial approximation (see the Herzog's metabolic cost function below).\cr
#'
#'
#'
#' \strong{Herzog's metabolic cost function in J/(kg*m)}:\cr
#'
#' \eqn{1 / (((1337.8 * abs(x[adj])^6) + (278.19 * abs(x[adj])^5) - (517.39 * abs(x[adj])^4) - (78.199 * abs(x[adj])^3) + (93.419 * abs(x[adj])^2) + (19.825 * abs(x[adj])) + 1.64) * N)  }\cr
#'
#' \strong{see} Herzog, I. (2016). Potential and Limits of Optimal Path Analysis. In A. Bevan & M. Lake (Eds.), Computational Approaches to Archaeological Spaces (pp. 179-211). New York: Routledge.
#' Herzog suggests to use this as a 6th degree polynomial approximation of the Minetti et al's cost function (see above).\cr
#'
#'
#'
#' \strong{Van Leusen's metabolic energy expenditure cost function (in Watts)}:\cr
#'
#' \eqn{ 1 / ((1.5 * W + 2.0 * (W + L) * (L / W)^2 + N * (W + L) * (1.5 * (V^2) + 0.35 * V * ((abs(x[adj])*100) + 10)))*N) }\cr
#'
#' which modifies the Pandolf et al.'s equation; \strong{see} Van Leusen, P. M. (2002). Pattern to process: methodological investigations into the formation and interpretation of spatial patterns in archaeological landscapes. University of Groningen.
#' \strong{Note} that, as per Herzog, I. (2013). Least-cost Paths - Some Methodological Issues, Internet Archaeology 36 (http://intarch.ac.uk/journal/issue36/index.html) and
#' unlike Van Leusen (2002), in the above equation slope is expressed in percent and speed in m/s; also, in the last bit of the equantion, 10 replaces
#' the value of 6 used by Van Leusen (as per Herzog 2013).\cr
#' As explained above, if \eqn{V} is set to 0 by the user, it is internally worked out on the basis of the Tobler function for on-path hiking; therefore, \eqn{V} will not be considered constant
#' throughout the analysed area, but will vary as function of the slope.\cr
#'
#' \strong{Note}: in the returned charts, the cost is transposed from Watts to Megawatts.\cr
#'
#'
#'
#' \strong{Llobera-Sluckin's metabolic energy expenditure cost function (in KJ/m)}:\cr
#'
#' \eqn{ 1 / ((2.635 + (17.37 * abs(x[adj])) + (42.37 * abs(x[adj])^2) - (21.43 * abs(x[adj])^3) + (14.93 * abs(x[adj])^4)) * N)  }\cr
#'
#' for which \strong{see} Llobera M.-Sluckin T.J. (2007). Zigzagging: Theoretical insights on climbing strategies, Journal of Theoretical Biology 249, 206-217.\cr
#'
#'
#'
#' \strong{Ardigo et al.'s metabolic energy expenditure cost function (in J/(kg*m))}:\cr
#'
#' \eqn{ 1 / ((1.866 * exp(4.911*abs(x[adj])) * V^2 - 3.773 * exp(3.416*abs(x[adj])) * V + (45.71*abs(x[adj]^2) + 18.90*abs(x[adj])) + 4.456) * N) }\cr
#'
#' \strong{see} Ardigo L.P., Saibene F., Minetti A.E. (2003), The optimal locomotion on gradients: walking, running or cycling?, in
#' Eur J Appl Physiol 90, 365-371. If \eqn{V} is set to 0 by the user, it is internally worked out on the basis of the Tobler function for on-path hiking; therefore,
#' \eqn{V} will not be considered constant throughout the analysed area, but will vary as function of the slope.\cr
#'
#'
#'
#' \strong{Hare's metabolic energy expenditure cost function (in cal/km)}:\cr
#'
#' \eqn{ 1 / ((48+30/(1/(6 * exp(-3.5 * abs(x[adj] + 0.05)))))*N) }\cr
#'
#' \strong{see} Hare T.S., Using Measures of Cost Distance in the Estimation of Polity Bounrdaries in the Post Classic
#' Yautepec valley, Mexico, in Journal of Archaeological Science 31 (2004). Energetic expenditure is expressed in calories per km; walking speed
#' is internally worked out from the DTM using the on-path Tobler's hiking function, which is expressed as its reciprocal; walking speed in km/h as per original Tobler's equation
#' and as requested by Hare's function.\cr
#'
#'
#'
#' \strong{Note} that the walking-speed-related cost functions listed above are used as they are, while the other functions are reciprocated.
#' This is done since "gdistance works with conductivity rather than the more usual approach using costs"; therefore
#' "we need inverse cost functions" (Nakoinz-Knitter (2016). "Modelling Human Behaviour in Landscapes". New York: Springer, p. 183).
#'  As a consequence, if we want to estimate time, we have to use the walking-speed functions as they are since the final accumulated values will correspond to the
#'  reciprocal of speed, i.e. pace. In the other cases, we have to use 1/cost-function to eventually get cost-function/1.\cr
#'
#'
#'
#' @param dtm Digital Terrain Model (RasterLayer class); if not provided, elevation data will be acquired online for the area enclosed by the 'studyplot' parameter (see Details).
#' @param origin location from which the cost surface is calculated (SpatialPointsDataFrame class).
#' @param destin location(s) to which least-cost path(s) is calculated (SpatialPointsDataFrame class).
#' @param studyplot polygon (SpatialPolygonDataFrame class) representing the study area for which online elevation data are aquired (see Details); NULL is default.
#' @param barrier area where the movement is inhibited (SpatialLineDataFrame or SpatialPolygonDataFrame class).
#' @param funct cost function to be used:\cr
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
#' \strong{h} uses the Hare's metabolic energy expenditure cost function (for all the mentioned cost functions, see Details);\cr
#' @param time time-unit expressed by the accumulated raster and by the isolines if Tobler's and other time-related cost functions are used;
#' 'h' for hour, 'm' for minutes.
#' @param outp type of output: 'raster' or 'contours' (see Details).
#' @param move number of directions in which cells are connected: 4 (rook's case), 8 (queen's case), 16 (knight and one-cell queen moves; default).
#' @param field value assigned to the cells coincidinng with the barrier (0 by default).
#' @param cogn.slp  TRUE or FALSE (default) if the user wants or does not want the 'cognitive slope' to be used in place of the real slope (see Details).
#' @param sl.crit critical slope (in percent), typically in the range 8-16 (10 by default) (used by the wheeled-vehicle cost function; see Details).
#' @param W walker's body weight (in Kg; 70 by default; used by the Pandolf's and Van Leusen's cost function; see Details).
#' @param L carried load weight (in Kg; 0 by default; used by the Pandolf's and Van Leusen's cost function; see Details).
#' @param N coefficient representing ease of movement (1 by default) (see Details).
#' @param V speed in m/s (1.2 by default) (used by the Pandolf et al.'s, Pandolf et al.s with correction factor, Van Leusen's, and Ardigo et al.'s cost function; if set to 0, it is internally worked out on the basis of Tobler on-path hiking function; see Details).
#' @param z zoom level for the elevation data downloaded from online sources (0 to 15; 9 by default) (see Details and \code{\link[elevatr]{get_elev_raster}}).
#' @param return.base TRUE or FALSE (default) if the user wants or does not want the least-cost paths back to the origin to be calculated and plotted (as dashed lines).
#' @param rb.lty line type used to represent the least-cost paths back to the origin in the returned plot (2 by default; dashed line; see 'lty' parameter in \code{\link[graphics]{par}}).
#' @param breaks contour interval; if no value is supplied, the interval is set by default to 1/10 of the range of values of the accumulated cost surface.
#' @param cont.lab if set to TRUE (default) display the labels of the contours over the accumulated cost surface.
#' @param destin.lab if set to TRUE (default) display the label(s) indicating the cost at the destination location(s).
#' @param cex.breaks set the size of the cost labels used in the contour plot (0.6 by default).
#' @param cex.lcp.lab set the size of the labels used in least-cost path(s) plot (0.6 by default).
#' @param graph.out TRUE (default) or FALSE if the user wants or does not want a graphical output to be generated.
#' @param transp set the transparency of the hillshade raster that is plotted over the rendered plots (0.5 by default).
#' @param oneplot TRUE (default) or FALSE if the user wants or does not want the plots displayed in a single window.
#' @param export TRUE or FALSE (default) if the user wants or does not want the outputs to be exported; if TRUE, the DTM, the cost-surface, and the accumulated cost surface are
#' exported as a GeoTiff file, while the isolines, the least-cost path(s), and a copy of the input destination locations (storing the cost measured at each location)
#' are exported as shapefile; all the exported files (excluding the DTM) will bear a suffix corresponding to the cost function selected by the user.
#' Note that the DTM is exported only if it was not provided by the user and downloaded by the function from online sources.
#'
#' @return The function returns a list storing the following components \itemize{
##'  \item{dtm: }{Digital Terrain Model ('RasterLayer' class)}
##'  \item{cost.surface: }{raster representing the cost-surface ('RasterLayer' class)}
##'  \item{accumulated.cost.raster: }{raster representing the accumualted cost ('RasterLayer' class)}
##'  \item{isolines: }{contour lines derived from the accumulated cost surface ('SpatialLinesDataFrame' class)}
##'  \item{LCPs: }{estimated least-cost paths ('SpatialLinesDataFrame' class)}
##'  \item{LCPs.back: }{estimated least-cost paths back to the origin ('SpatialLinesDataFrame' class)}
##'  \item{LCPs$length: }{length of each least-cost path (units depend on the unit used in the input DTM)}
##'  \item{LCPs.back$length: }{length of each least-cost path back to the origin (units depend on the unit used in the input DTM)}
##'  \item{dest.loc.w.cost: }{copy of the input destination location(s) dataset with a new variable ('cost') added; if
##'  the cost is expressed in terms of time, the 'cost' variable will store the time values in decimal numbers, while another variable named
##'  'cost_hms' will store the time values in sexagesimmal numbers (hours, minutes, seconds)}
##'  \item{conductance: }{conductance 'Transitional Layer', returned because internally used by the 'movenetw()' function}
##' }
##'
##'
#' @keywords movecost
#' @export
#' @importFrom raster ncell mask crop raster hillShade terrain
#' @importFrom elevatr get_elev_raster
#' @importFrom chron times
#' @importFrom grDevices terrain.colors topo.colors grey
#' @importFrom graphics layout par
#' @importFrom leastcostpath create_barrier_cs
#'
#'
#' @examples
#' # load a sample Digital Terrain Model
#' data(volc)
#'
#' # load a sample start location on the above DTM
#' data(volc.loc)
#'
#' # load the sample destination locations on the above DTM
#' data(destin.loc)
#'
#' # calculate walking-time isochrones based on the on-path Tobler's hiking function (default),
#' # setting the time unit to hours and the isochrones interval to 0.05 hour;
#' # also, since destination locations are provided,
#' # least-cost paths from the origin to the destination locations will be calculated
#' # and plotted; 8-directions move is used
#'
#' result <- movecost(dtm=volc, origin=volc.loc, destin=destin.loc, move=8, breaks=0.05)
#'
#'
#' # same as above, but using the Irmischer-Clarke's hiking function (male, on-path)
#'
#' result <- movecost(dtm=volc, origin=volc.loc, destin=destin.loc, funct="icmonp",
#' move=8, breaks=0.05)
#'
#'
#' # same as above, but using the 'cognitive slope'
#'
#' result <- movecost(dtm=volc, origin=volc.loc, destin=destin.loc, funct="icmonp",
#' move=8, breaks=0.05, cogn.slp=TRUE)
#'
#'
#' # calculate accumulated cost surface and the least-cost path between the
#' # origin and one destination, and also calculate the LCP back to the origin
#'
#' results <- movecost(dtm=volc, origin=volc.loc, destin=destin.loc[2,], move=8, return.base = TRUE)
#'
#'
#' @seealso \code{\link[elevatr]{get_elev_raster}}, \code{\link{movecorr}}, \code{\link{movebound}}, \code{\link{movealloc}},  \code{\link{movecomp}}, \code{\link{movenetw}}
#'
#'
movecost <- function (dtm=NULL, origin, destin=NULL, studyplot=NULL, barrier=NULL, funct="t", time="h", outp="r", move=16, field=0, cogn.slp=FALSE, sl.crit=10, W=70, L=0, N=1, V=1.2, z=9, return.base=FALSE, rb.lty=2, breaks=NULL, cont.lab=TRUE, destin.lab=TRUE, cex.breaks=0.6, cex.lcp.lab=0.6, graph.out=TRUE, transp=0.5, oneplot=TRUE, export=FALSE){

  #deactivate the warning messages because a warning that can be safely ignored will be produced by the procedure
  #used to get slope as rise over run
  options(warn = -1)

  #if no dtm is provided
  if (is.null(dtm)==TRUE) {
    #get the elvation data using the elevatr's get_elev_raster() function, using the studyplot dataset (SpatialPolygonDataFrame)
    #to select the area whose elevation data are to be downloaded;
    #z sets the resolution of the elevation datataset
    elev.data <- elevatr::get_elev_raster(studyplot, z = z, verbose=FALSE, override_size_check = TRUE)
    #crop the elevation dataset to the exact boundary of the studyplot dataset
    dtm <- raster::crop(elev.data, studyplot)
  }

  #calculate the altitudinal difference between adjacent cells
  altDiff <- function(x){x[2] - x[1]}
  hd <- gdistance::transition(dtm, altDiff, directions=move, symm=FALSE)

  #use the geoCorrection function to divide the altitudinal difference by the distance between cells
  #so getting slope as rise over run
  slope <- gdistance::geoCorrection(hd)

  #if 'cogn.slp' is set to TRUE, positive slope values (ie. up-hill movement) multiplied by 1.99
  # and negative slope values (ie. down-hill movement) multiplied by 2.31
  if (cogn.slp==TRUE) {
    slope@transitionMatrix@x <- ifelse(slope@transitionMatrix@x > 0, slope@transitionMatrix@x * 1.99, slope@transitionMatrix@x * 2.31)
  }

  #define different types of cost functions and set the appropriate text to be used for subsequent plotting
  if (funct=="t") {
    #Tobler's hiking function; kmh
    cost_function <- function(x){ (6 * exp(-3.5 * abs(x[adj] + 0.05))) * (1/N) }

    #set the labels to be used within the returned plot
    main.title <- paste0("Walking-time isochrones (in ", time, ") around origin")
    sub.title <- paste0("Walking-time based on the Tobler's on-path hiking function \nterrain factor N=", N)
    legend.cost <- paste0("walking-time (", time,")")
    sub.title.lcp.plot <- paste0("LCP(s) and walking-time distance(s) based on the Tobler's on-path hiking function \nterrain factor N=", N, "\nblack dot=start location\n red dot(s)=destination location(s)")
  }

  if (funct=="tofp") {
    #Tobler's hiking function off-path routes; kmh
    #note that the multiplier 0.6 suggested by Tobler is meant to reduce the off-path walking speed
    cost_function <- function(x){(6 * exp(-3.5 * abs(x[adj] + 0.05))) * 0.6}

    #set the labels to be used within the returned plot
    main.title <- paste0("Walking-time isochrones (in ", time, ") around origin")
    sub.title <- "Walking-time based on the Tobler's off-path hiking function"
    legend.cost <- paste0("walking-time (", time,")")
    sub.title.lcp.plot <- paste0("LCP(s) and walking-time distance(s) based on the Tobler's off-path hiking function \nblack dot=start location\n red dot(s)=destination location(s)")
  }

  if(funct=="mp") {
    #Marquez-Perez et al.'s modified Tobler hiking function; kmh
    cost_function <- function(x){ (4.8 * exp(-5.3 * abs((x[adj] * 0.7) + 0.03))) * (1/N) }

    #set the labels to be used within the returned plot
    main.title <- paste0("Walking-time isochrones (in ", time, ") around origin")
    sub.title <- paste0("Walking-time based on the Marquez-Perez et al.'s modified Tobler hiking function \n terrain factor N=", N)
    legend.cost <- paste0("walking-time (", time,")")
    sub.title.lcp.plot <- paste0("LCP(s) and walking-time distance(s) based on the Marquez-Perez et al.'s modified Tobler hiking function \n terrain factor N=", N, "\nblack dot=start location\n red dot(s)=destination location(s)")
  }

  if(funct=="alb") {
    #Alberti's modified Tobler hiking function, adapted for animal foraging excursions
    cost_function <- function(x){(6 * exp(-3.5 * abs(x[adj] + 0.05))) * 0.25}

    #set the labels to be used within the returned plot
    main.title <- paste0("Walking-time isochrones (in ", time, ") around origin")
    sub.title <- "Walking-time based on the Alberti's modified Tobler hiking function"
    legend.cost <- paste0("walking-time (", time,")")
    sub.title.lcp.plot <- paste0("LCP(s) and walking-time distance(s) based on the Alberti's modified Tobler hiking function \nblack dot=start location\n red dot(s)=destination location(s)")
  }

  if(funct=="icmonp") {
    #Irmischer-Clarke's modified Tobler hiking function; originally in m/s (male, on-path);
    # the formula is reshaped (multiplied by 3.6) below to turn it into kmh for consistency with the other Tobler-related cost functions;
    # Slope in percent.
    cost_function <- function(x){ ((0.11 + exp(-(abs(x[adj])*100 + 5)^2 / (2 * 30^2))) * 3.6) * (1/N) }

    #set the labels to be used within the returned plot
    main.title <- paste0("Walking-time isochrones (in ", time, ") around origin")
    sub.title <- paste0("Walking-time based on the (male, on-path) Irmischer-Clarke's hiking function \n terrain factor N=", N)
    legend.cost <- paste0("walking-time (", time,")")
    sub.title.lcp.plot <- paste0("LCP(s) and walking-time distance(s) based on the (male, on-path) Irmischer-Clarke's hiking function \n terrain factor N=", N, "\nblack dot=start location\n red dot(s)=destination location(s)")
  }

  if(funct=="icmoffp") {
    #Irmischer-Clarke's modified Tobler hiking function; originally in m/s (male, off-path);
    # the formula is reshaped (multiplied by 3.6) below to turn it into kmh for consistency with the other Tobler-related cost functions;
    # Slope in percent.
    cost_function <- function(x){ (0.11 + 0.67 * exp(-(abs(x[adj])*100 + 2)^2 / (2 * 30^2))) * 3.6 }

    #set the labels to be used within the returned plot
    main.title <- paste0("Walking-time isochrones (in ", time, ") around origin")
    sub.title <- "Walking-time based on the (male, off-path) Irmischer-Clarke's hiking function"
    legend.cost <- paste0("walking-time (", time,")")
    sub.title.lcp.plot <- paste0("LCP(s) and walking-time distance(s) based on the (male, off-path) Irmischer-Clarke's hiking function \nblack dot=start location\n red dot(s)=destination location(s)")
  }

  if(funct=="icfonp") {
    #Irmischer-Clarke's modified Tobler hiking function; originally in m/s (female, on-path);
    # the formula is reshaped (multiplied by 3.6) below to turn it into kmh for consistency with the other Tobler-related cost functions;
    # Slope in percent.
    cost_function <- function(x){ ((0.95 * (0.11 + exp(-(abs(x[adj]) * 100 + 5)^2/(2 * 30^2)))) * 3.6) * (1/N) }

    #set the labels to be used within the returned plot
    main.title <- paste0("Walking-time isochrones (in ", time, ") around origin")
    sub.title <- paste0("Walking-time based on the (female, on-path) Irmischer-Clarke's hiking function \n terrain factor N=", N)
    legend.cost <- paste0("walking-time (", time,")")
    sub.title.lcp.plot <- paste0("LCP(s) and walking-time distance(s) based on the (female, on-path) Irmischer-Clarke's hiking function \n terrain factor N=", N, "\nblack dot=start location\n red dot(s)=destination location(s)")
  }

  if(funct=="icfoffp") {
    #Irmischer-Clarke's modified Tobler hiking function; originally in m/s (female, off-path);
    # the formula is reshaped (multiplied by 3.6) below to turn it into kmh for consistency with the other Tobler-related cost functions;
    # Slope in percent.
    cost_function <- function(x){ (0.95 * (0.11 + 0.67 * exp(-(abs(x[adj]) * 100 + 2)^2/(2 * 30^2)))) * 3.6 }

    #set the labels to be used within the returned plot
    main.title <- paste0("Walking-time isochrones (in ", time, ") around origin")
    sub.title <- "Walking-time based on the (female, off-path) Irmischer-Clarke's hiking function"
    legend.cost <- paste0("walking-time (", time,")")
    sub.title.lcp.plot <- paste0("LCP(s) and walking-time distance(s) based on the (female, off-path) Irmischer-Clarke's hiking function \nblack dot=start location\n red dot(s)=destination location(s)")
  }

  if(funct=="gkrs") {
    #Garmy, Kaddouri, Rozenblat, Schneider's slope-dependant hiking function; speed in kmh;
    #slope is originally in degrees; (atan(abs(x[adj]))*180/pi) turns rise/run into degrees
    cost_function <- function(x){ (4 * exp(-0.008 * ((atan(abs(x[adj]))*180/pi)^2))) * (1/N) }

    #set the labels to be used within the returned plot
    main.title <- paste0("Walking-time isochrones (in ", time, ") around origin")
    sub.title <- paste0("Walking-time based on the Garmy et al.'s hiking function \n terrain factor N=", N)
    legend.cost <- paste0("walking-time (", time,")")
    sub.title.lcp.plot <- paste0("LCP(s) and walking-time distance(s) based on the Garmy et al.'s hiking function \n terrain factor N=", N, "\nblack dot=start location\n red dot(s)=destination location(s)")
  }

  if(funct=="ug") {
    #Antonio Uriarte Gonzalez's slope-dependant walking-time cost function;
    # the cost function is originally expressed in seconds; here it is the reciprocal of time (1/T) that is used in order to eventually get
    #T/1. Slope is in percent.
    #Note: unlike the original formula, here the pixel resolution is not taken into account since 'gdistance' takes care of the cells' dimension
    #when calculating accumulated costs.
    cost_function <- function(x){ 1 / ((0.0277 * (abs(x[adj])*100) + 0.6115) * N) }

    #set the labels to be used within the returned plot
    main.title <- paste0("Walking-time isochrones (in ", time, ") around origin")
    sub.title <- paste0("Walking-time based on the Uriarte Gonzalez's hiking function \n terrain factor N=", N)
    legend.cost <- paste0("walking-time (", time,")")
    sub.title.lcp.plot <- paste0("LCP(s) and walking-time distance(s) based on the Uriarte Gonzalez's hiking function \n terrain factor N=", N, "\nblack dot=start location\n red dot(s)=destination location(s)")
  }

  if(funct=="r") {
    #Rees' slope-dependant walking-time cost function;
    #the cost function is originally expressed in terms of time (1/v in Rees' publication);
    #here it is the reciprocal of time (i.e. speed) that is used in order to eventually get
    #the reciprocal of speed (i.e. time). Slope is as originally expressed in Rees' publication (i.e. rise over run);
    #the speed is originally expressed in m/s, so here has been transposed to kmh (i.e., multiplied by 3.6)
    cost_function <- function(x){ ((1 / (0.75 + 0.09 * abs(x[adj]) + 14.6 * (abs(x[adj]))^2)) * 3.6) * (1/N) }

    #set the labels to be used within the returned plot
    main.title <- paste0("Walking-time isochrones (in ", time, ") around origin")
    sub.title <- paste0("Walking-time based on the Rees' hiking function \n terrain factor N=", N)
    legend.cost <- paste0("walking-time (", time,")")
    sub.title.lcp.plot <- paste0("LCP(s) and walking-time distance(s) based on the Rees' hiking function \n terrain factor N=", N, "\nblack dot=start location\n red dot(s)=destination location(s)")
  }

  if(funct=="ks") {
    #Kondo-Saino's modifier Tobler's hiking function.
    #It expresses walking speed in KmH; slope is rise/run.
    cost_function <- function(x){ ifelse(abs(x[adj]) >= -0.07, (5.1 * exp(-2.25 * abs(x[adj] + 0.07))) * (1/N), (5.1 * exp(-1.5 * abs(x[adj] + 0.07)))) * (1/N) }

    #set the labels to be used within the returned plot
    main.title <- paste0("Walking-time isochrones (in ", time, ") around origin")
    sub.title <- paste0("Walking-time based on the Kondo-Seino's hiking function \n terrain factor N=", N)
    legend.cost <- paste0("walking-time (", time,")")
    sub.title.lcp.plot <- paste0("LCP(s) and walking-time distance(s) based on the  Kondo-Seino's hiking function \n terrain factor N=", N, "\nblack dot=start location\n red dot(s)=destination location(s)")
  }

  if(funct=="ree") {
    #relative energetic expenditure;
    # to calculate tangent of degrees (as requested by the cost function) we must first convert degrees to radians by multypling by pi/180;
    #(atan(abs(x[adj]))*180/pi) turns rise/run into degrees, which are then converted into radians before calculating the tangent
    cost_function <- function(x){ 1 / ((tan((atan(abs(x[adj]))*180/pi)*pi/180) / tan (1*pi/180)) * N) }

    #set the labels to be used within the returned plot
    main.title <- "Accumulated cost isolines around origin"
    sub.title <- paste0("Cost based on the relative energetic expenditure cost function \n terrain factor N=", N)
    legend.cost <- "cost"
    sub.title.lcp.plot <- paste0("LCP(s) and cost distance(s) based on the slope-based relative energetic expenditure cost function \n terrain factor N=", N, "\nblack dot=start location\n red dot(s)=destination location(s)")
  }

  if(funct=="hrz") {
    #Herzog metabolic cost function in J/(kg*m);
    #rise/run is requested by the cost function
    cost_function <- function(x){ 1 / (((1337.8 * abs(x[adj])^6) + (278.19 * abs(x[adj])^5) - (517.39 * abs(x[adj])^4) - (78.199 * abs(x[adj])^3) + (93.419 * abs(x[adj])^2) + (19.825 * abs(x[adj])) + 1.64) * N) }

    #set the labels to be used within the returned plot
    main.title <- "Accumulated cost isolines around origin"
    sub.title <- paste0("Cost based on the Herzog's metabolic cost function \n cost in J / (Kg*m) \n terrain factor N=", N)
    legend.cost <- "metabolic cost J / (Kg*m)"
    sub.title.lcp.plot <- paste0("LCP(s) and cost distance(s) based on the Herzog's metabolic cost function \ncost in J / (Kg*m); terrain factor N=", N, "\nblack dot=start location\n red dot(s)=destination location(s)")
  }

  if(funct=="wcs") {
    #wheeled-vehicle critical slope cost function; the slope is expressed in percent
    cost_function <- function(x){ 1 / ((1 + ((abs(x[adj])*100) / sl.crit)^2) * N) }

    #set the labels to be used within the returned plot
    main.title <- "Accumulated cost isolines around origin"
    sub.title <- paste0("Cost based on the wheeled-vehicle critical slope cost function \ncritical slope set to ", sl.crit, " percent \n terrain factor N=", N)
    legend.cost <- "cost"
    sub.title.lcp.plot <- paste0("LCP(s) and cost distance(s) based on the wheeled-vehicle critical slope cost function \ncritical slope set to ", sl.crit, " percent; terrain factor N=", N, "\nblack dot=start location\n red dot(s)=destination location(s)")
  }

  if(funct=="vl") {
    #Van Leusen's metabolic energy expenditure cost function
    #note: V is velocity in m/s; the slope is in percent
    #if V is set to 0 by the user, it is worked out from the DTM using the on-path Tobler hiking function
    #which is expressed as its reciprocal and is multiplied by 0.278 to turn kmh to m/s
    if (V==0) {
      cost_function <- function(x){ 1 / ((1.5 * W + 2.0 * (W + L) * (L / W)^2 + N * (W + L) * (1.5 * ((1 / ((6 * exp(-3.5 * abs(x[adj] + 0.05))) * 0.278))^2) + 0.35 * (1 / ((6 * exp(-3.5 * abs(x[adj] + 0.05))) * 0.278)) * ((abs(x[adj])*100) + 10)))*N) }
    } else {
      cost_function <- function(x){ 1 / ((1.5 * W + 2.0 * (W + L) * (L / W)^2 + N * (W + L) * (1.5 * (V^2) + 0.35 * V * ((abs(x[adj])*100) + 10)))*N) }
    }

    #set the labels to be used within the returned plot
    main.title <- "Accumulated cost isolines around origin"
    legend.cost <- "energy expenditure cost (Megawatts)"
    if (V==0) {
      sub.title <- paste0("Cost based on the Van Leusen's metabolic energy expenditure cost function \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V is based on the Tobler on-path hiking function")
      sub.title.lcp.plot <- paste0("LCP(s) and cost distance(s) based on the Van Leusen's metabolic energy expenditure cost function \n cost in Megawatts; parameters: W: ", W, "; L: ", L, "; N: ", N, "; V is based on the Tobler on-path hiking function \nblack dot=start location\n red dot(s)=destination location(s)")
    } else {
      sub.title <- paste0("Cost based on the Van Leusen's metabolic energy expenditure cost function \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V: ", V)
      sub.title.lcp.plot <- paste0("LCP(s) and cost distance(s) based on the Van Leusen's metabolic energy expenditure cost function \n cost in Megawatts; parameters: W: ", W, "; L: ", L, "; N: ", N, "; V: ", V, "\nblack dot=start location\n red dot(s)=destination location(s)")
    }
  }

  if(funct=="p") {
    #Pandolf et al.'s metabolic energy expenditure cost function
    #note: V is velocity in m/s; the slope is expressed in percent
    #if V is set to 0 by the user, it is worked out from the DTM using the on-path Tobler hiking function
    #which is expressed as its reciprocal and is multiplied by 0.278 to turn kmh to m/s
    if (V==0) {
      cost_function <- function(x){ 1 / ((1.5 * W + 2.0 * (W + L) * (L / W)^2 + N * (W + L) * (1.5 * ((1 / ((6 * exp(-3.5 * abs(x[adj] + 0.05))) * 0.278))^2) + 0.35 * (1 / ((6 * exp(-3.5 * abs(x[adj] + 0.05))) * 0.278)) * (abs(x[adj])*100)))*N) }
    } else {
      cost_function <- function(x){ 1 / ((1.5 * W + 2.0 * (W + L) * (L / W)^2 + N * (W + L) * (1.5 * (V^2) + 0.35 * V * (abs(x[adj])*100)))*N) }
    }

    #set the labels to be used within the returned plot
    main.title <- "Accumulated cost isolines around origin"
    legend.cost <- "energy expenditure cost (Megawatts)"
    if (V==0) {
      sub.title <- paste0("Cost based on the Pandolf et al.'s metabolic energy expenditure cost function \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V is based on the Tobler on-path hiking function")
      sub.title.lcp.plot <- paste0("LCP(s) and cost distance(s) based on the Pandolf et al.'s metabolic energy expenditure cost function \n cost in Megawatts; parameters: W: ", W, "; L: ", L, "; N: ", N, "; V is based on the Tobler on-path hiking function \nblack dot=start location\n red dot(s)=destination location(s)")
    } else {
      sub.title <- paste0("Cost based on the Pandolf et al.'s metabolic energy expenditure cost function \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V: ", V)
      sub.title.lcp.plot <- paste0("LCP(s) and cost distance(s) based on the Pandolf et al.'s metabolic energy expenditure cost function \n cost in Megawatts; parameters: W: ", W, "; L: ", L, "; N: ", N, "; V: ", V, "\nblack dot=start location\n red dot(s)=destination location(s)")
    }
  }

  if(funct=="pcf") {
    #Pandolf et al.'s metabolic energy expenditure cost function with CORRECTION FACTOR
    #note: V is velocity in m/s; the slope is expressed in percent
    #if V is set to 0 by the user, it is worked out from the DTM using the on-path Tobler hiking function
    #which is expressed as its reciprocal and is multiplied by 0.278 to turn kmh to m/s
    if (V==0) {
      cost_function <- function(x){ ifelse(abs(x[adj])*100 > 0 , 1 / (1.5 * W + 2.0 * (W+L) * (L/W)^2 + N * (W+L) * (1.5 * (1 / ((6 * exp(-3.5 * abs(x[adj] + 0.05))) * 0.278))^2 + 0.35 * (1 / ((6 * exp(-3.5 * abs(x[adj] + 0.05))) * 0.278)) * (abs(x[adj])*100))),
                                           1 / ((1.5 * W + 2.0 * (W+L) * (L/W)^2 + N * (W+L) * (1.5 * (1 / ((6 * exp(-3.5 * abs(x[adj] + 0.05))) * 0.278))^2 + 0.35 * (1 / ((6 * exp(-3.5 * abs(x[adj] + 0.05))) * 0.278)) * (abs(x[adj])*100))) - (N * ((abs(x[adj])*100) * (W+L) * (1 / ((6 * exp(-3.5 * abs(x[adj] + 0.05))) * 0.278))/3.5) - ((W+L) * ((abs(x[adj])*100)+6)^2/W) + (25-(1 / ((6 * exp(-3.5 * abs(x[adj] + 0.05))) * 0.278))^2)))) }
    } else {
      cost_function <- function(x){ ifelse(abs(x[adj])*100 > 0 , 1 / (1.5 * W + 2.0 * (W+L) * (L/W)^2 + N * (W+L) * (1.5 * V^2 + 0.35 * V * (abs(x[adj])*100))),
                                           1 / ((1.5 * W + 2.0 * (W+L) * (L/W)^2 + N * (W+L) * (1.5 * V^2 + 0.35 * V * (abs(x[adj])*100))) - (N * ((abs(x[adj])*100) * (W+L) * V/3.5) - ((W+L) * ((abs(x[adj])*100)+6)^2/W) + (25-V^2)))) }
    }

    #set the labels to be used within the returned plot
    main.title <- "Accumulated cost isolines around origin"
    legend.cost <- "energy expenditure cost (Megawatts)"
    if (V==0) {
      sub.title <- paste0("Cost based on the Pandolf et al.'s metabolic energy expenditure cost function with correction factor \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V is based on the Tobler on-path hiking function")
      sub.title.lcp.plot <- paste0("LCP(s) and cost distance(s) based on the Pandolf et al.'s metabolic energy expenditure cost function (with correction factor) \n cost in Megawatts; parameters: W: ", W, "; L: ", L, "; N: ", N, "; V is based on the Tobler on-path hiking function \nblack dot=start location\n red dot(s)=destination location(s)")
    } else {
      sub.title <- paste0("Cost based on the Pandolf et al.'s metabolic energy expenditure cost function with correction factor \nparameters: W: ", W, "; L: ", L, "; N: ", N, "; V: ", V)
      sub.title.lcp.plot <- paste0("LCP(s) and cost distance(s) based on the Pandolf et al.'s metabolic energy expenditure cost function with correction factor \n cost in Megawatts; parameters: W: ", W, "; L: ", L, "; N: ", N, "; V: ", V, "\nblack dot=start location\n red dot(s)=destination location(s)")
    }
  }

  if(funct=="ls") {
    #Llobera-Sluckin's metabolic energy expenditure cost function (in KJ/m)
    cost_function <- function(x){ 1 / ((2.635 + (17.37 * abs(x[adj])) + (42.37 * abs(x[adj])^2) - (21.43 * abs(x[adj])^3) + (14.93 * abs(x[adj])^4)) * N) }

    #set the labels to be used within the returned plot
    main.title <- "Accumulated cost isolines around origin"
    sub.title <- paste0("Cost based on the Llobera-Sluckin's metabolic energy expenditure cost function \n terrain factor N=", N)
    legend.cost <- "energy expenditure cost (KJ/m)"
    sub.title.lcp.plot <- paste0("LCP(s) and cost distance(s) based on the Llobera-Sluckin's metabolic energy expenditure cost function \n cost in KJ/m; terrain factor N=", N, "\nblack dot=start location\n red dot(s)=destination location(s)")
  }

  if(funct=="b") {
    #Bellavia 2002's cost function; it measures abstract cost;
    #slope is originally in degrees; (atan(abs(x[adj]))*180/pi) turns rise/run into degrees
    cost_function <- function(x){ 1 / (((atan(abs(x[adj]))*180/pi)+1) * N) }

    #set the labels to be used within the returned plot
    main.title <- "Accumulated cost isolines around origin"
    sub.title <- paste0("Cost based on the Bellavia's cost function \n terrain factor N=", N)
    legend.cost <- "cost"
    sub.title.lcp.plot <- paste0("LCP(s) and cost distance(s) based on the Bellavia's cost function \n terrain factor N=", N, "\nblack dot=start location\n red dot(s)=destination location(s)")
  }

  if(funct=="m") {
    #Minetti et al.'s cost function
    #slope as gradient (i.e., rise/run); cost in J/(kg*m)
    cost_function <- function(x){ 1 / (((280.5 * abs(x[adj])^5) - (58.7 * abs(x[adj])^4) - (76.8 * abs(x[adj])^3) + (51.9 * abs(x[adj])^2) + (19.6 * abs(x[adj])) + 2.5) * N) }

    #set the labels to be used within the returned plot
    main.title <- "Accumulated cost isolines around origin"
    sub.title <- paste0("Cost based on the Minetti et al.'s metabolic cost function \n cost in J / (Kg*m) \n terrain factor N=", N)
    legend.cost <- "metabolic cost J / (Kg*m)"
    sub.title.lcp.plot <- paste0("LCP(s) and cost distance(s) based on the Minetti et al.'s metabolic cost function \ncost in J / (Kg*m); terrain factor N=", N, "\nblack dot=start location\n red dot(s)=destination location(s)")
  }

  if(funct=="ma") {
    #Marin Arroyo's slope-dependant walking-time cost function;
    #the cost function is originally expressed in seconds; here it is the reciprocal of time (1/T) that is used in order to eventually get
    #T/1. Slope is in percent.
    #Note: unlike the original formula, here d (distance travelled in meter) is not taken into account since 'gdistance' takes care of the cells' dimension
    #when calculating accumulated costs.
    cost_function <- function(x){ ifelse((abs(x[adj])*100) < 0, 1 / ((0.6 * ((abs(x[adj])*100)/23+1))*N), 1 / ((0.6 * ((abs(x[adj])*100)/11+1))*N)) }

    #set the labels to be used within the returned plot
    main.title <- paste0("Walking-time isochrones (in ", time, ") around origin")
    sub.title <- paste0("Walking-time based on the Marin Arroyo's hiking function \n terrain factor N=", N)
    legend.cost <- paste0("walking-time (", time,")")
    sub.title.lcp.plot <- paste0("LCP(s) and walking-time distance(s) based on the Marin Arroyo's hiking function \n terrain factor N=", N, "\nblack dot=start location\n red dot(s)=destination location(s)")
  }

  if(funct=="a") {
    #Ardigo et al's energy cost function; slope as rise/run (gradient)
    #cost in J/(kg*m)
    #if V is set to 0 by the user, it is worked out from the DTM using the on-path Tobler hiking function
    #which is expressed as its reciprocal and is multiplied by 0.278 to turn kmh to m/s
    if (V==0) {
      cost_function <- function(x){ 1 / ((1.866 * exp(4.911*abs(x[adj])) * (1 / ((6 * exp(-3.5 * abs(x[adj] + 0.05))) * 0.278))^2 - 3.773 * exp(3.416*abs(x[adj])) * (1 / ((6 * exp(-3.5 * abs(x[adj] + 0.05))) * 0.278)) + (45.71*abs(x[adj]^2) + 18.90*abs(x[adj])) + 4.456) * N) }
    } else {
      cost_function <- function(x){ 1 / ((1.866 * exp(4.911*abs(x[adj])) * V^2 - 3.773 * exp(3.416*abs(x[adj])) * V + (45.71*abs(x[adj]^2) + 18.90*abs(x[adj])) + 4.456) * N) }
    }

    #set the labels to be used within the returned plot
    main.title <- "Accumulated cost isolines around origin"
    legend.cost <- "metabolic cost J / (Kg*m)"
    if (V==0) {
      sub.title <- paste0("Cost based on the Ardigo et al.'s metabolic cost function \n cost in J / (Kg*m) \nparameters: N=", N, "; V is based on the Tobler on-path hiking function")
      sub.title.lcp.plot <- paste0("LCP(s) and cost distance(s) based on the Ardigo et al.'s metabolic cost function \ncost in J / (Kg*m); parameters: N=", N, "; V is based on the Tobler on-path hiking function \nblack dot=start location\n red dot(s)=destination location(s)")
    } else {
      sub.title <- paste0("Cost based on the Ardigo et al.'s metabolic cost function \n cost in J / (Kg*m) \nparameters: N=", N, "; V=", V)
      sub.title.lcp.plot <- paste0("LCP(s) and cost distance(s) based on the Ardigo et al.'s metabolic cost function \ncost in J / (Kg*m); parameters: N=", N, "; V=", V, "\nblack dot=start location\n red dot(s)=destination location(s)")
    }
  }

  if(funct=="h") {
    #energetic expenditure in calories per km, following Hare 2004;
    #walking speed is worked out from the DTM using the on-path Tobler hiking function,
    #which is expressed as its reciprocal; speed in km/h as per original Tobler's function and as requested by Hare's function
    cost_function <- function(x){ 1 / ((48+30/(1/(6 * exp(-3.5 * abs(x[adj] + 0.05)))))*N) }

    #set the labels to be used within the returned plot
    main.title <- "Accumulated cost isolines around origin"
    sub.title <- paste0("Cost based on the Hare's metabolic cost function \n cost in cal/km \n terrain factor N=", N)
    legend.cost <- "metabolic cost cal/km"
    sub.title.lcp.plot <- paste0("LCP(s) and cost distance(s) based on the Hare's metabolic cost function \ncost in cal/km; terrain factor N=", N, "\nblack dot=start location\n red dot(s)=destination location(s)")
  }

  if(funct=="e") {
    #Eastman's cost function; it measures abstract cost;
    #slope is originally in degrees; (atan(abs(x[adj]))*180/pi) turns rise/run into degrees
    cost_function <- function(x){ 1 / ((0.031 * (atan(abs(x[adj]))*180/pi)^2 - 0.025 * (atan(abs(x[adj]))*180/pi) + 1)*N) }

    #set the labels to be used within the returned plot
    main.title <- "Accumulated cost isolines around origin"
    sub.title <- paste0("Cost based on the Eastman's cost function \n terrain factor N=", N)
    legend.cost <- "cost"
    sub.title.lcp.plot <- paste0("LCP(s) and cost distance(s) based on the Eastman's cost function \n terrain factor N=", N, "\nblack dot=start location\n red dot(s)=destination location(s)")
  }

  #cost calculation for walking-speed-based cost functions
  if (funct=="t" | funct=="tofp" | funct=="mp" | funct=="icmonp" | funct=="icmoffp" | funct=="icfonp" | funct=="icfoffp" | funct=="alb" | funct=="gkrs" | funct=="r" | funct=="ks") {

    #restrict the speed calculation to adjacent cells by creating an index for adjacent cells (adj) with the function 'adjacent'
    adj <- raster::adjacent(dtm, cells=1:raster::ncell(dtm), pairs=TRUE, directions=move)

    speed <- slope

    #apply the cost function to the adjacent cells of the speed dataset, which is equal to the slope dataset as per previous step
    speed[adj] <- cost_function(slope)

    #turn the walking speed from kmh to ms (0.278=1000/3600)
    speed <- speed * 0.278

    #correct the speed values taking into account the distance between cell centers
    Conductance <- gdistance::geoCorrection(speed)

    #if the 'barrier' parameter in not NULL...
    if (is.null(barrier)==FALSE) {
      #...create a barrier transitional layers via the 'create_barrier_cs' function
      barrier.TL <- leastcostpath::create_barrier_cs(dtm, barrier, neighbours=move, field=field)
      #use map algebra to "incorporate" the barrier to the conductance layer
      Conductance <- Conductance*barrier.TL
    }

    #cost-surface raster to be exported (the division turns ms back to kmh)
    cost.surface.to.export <- raster::raster(speed) / 0.278
  }

  #cost calculation for OTHER TYPES of cost functions;
  #note the Uriarte Gonzalez's and Marin Arroyo's slope-dependant walking-time cost function are in this group since (unlike the above functions)
  #they expresses cost as time NOT speed
  if (funct=="ree" | funct=="hrz" | funct=="wcs" | funct=="vl" | funct=="p" | funct=="ug" | funct=="ls" | funct=="b" | funct=="m" | funct=="ma" | funct=="a" | funct=="pcf" | funct=="h" | funct=="e") {

    #restrict the cost calculation to adjacent cells by creating an index for adjacent cells (adj) with the function 'adjacent'
    adj <- raster::adjacent(dtm, cells=1:raster::ncell(dtm), pairs=TRUE, directions=move)

    cost <- slope

    #apply the cost function to the adjacent cells of the cost dataset, which is equal to the slope dataset as per previous step
    cost[adj] <- cost_function(slope)

    #correct the cost values taking into account the distance between cell centers
    Conductance <- gdistance::geoCorrection(cost)

    #if the 'barrier' parameter in not NULL...
    if (is.null(barrier)==FALSE) {
      #...create a barrier transitional layers via the 'create_barrier_cs' function
      barrier.TL <- leastcostpath::create_barrier_cs(dtm, barrier, neighbours=move, field=field)
      #use map algebra to "incorporate" the barrier to the conductance layer
      Conductance <- Conductance*barrier.TL
    }

    #cost-surface raster to be exported
    cost.surface.to.export <- raster::raster(cost)
  }

  #accumulate the cost outwards from the origin
  accum_final <- gdistance::accCost(Conductance, sp::coordinates(origin))

  #if user select the Tobler's, the modified Tobler's, the Irmischer-Clarke's,
  #the Uriarte Gonzalez's, the Alberti's function, or other speed-related functions, turn seconds into the user-defined time-scale
  if (funct=="t" | funct=="tofp" | funct=="mp" | funct=="icmonp" | funct=="icmoffp" | funct=="icfonp" | funct=="icfoffp" | funct=="ug" | funct=="alb" | funct=="gkrs" | funct=="r" | funct=="ks" | funct=="ma"){
    if (time=="h") {
      #turn seconds into hours
      accum_final <- accum_final / 3600
    } else {
      #turn seconds into minutes
      accum_final <- accum_final / 60
    }
  }

  #if user select the Val Leusen's, the Pandolf et al.'s, or the Pandolf et al.'s with correction factor function,
  #turn the cost from Watts to Megawatts
  if (funct=="vl" | funct=="p" | funct=="pcf") {
    accum_final <- accum_final / 1000000
  }

  #if no break value is entered, set the breaks to one tenth of the range of the values of the final accumulated cost surface
  if(is.null(breaks)==TRUE){
    #exclude the inf values from the calculation
    breaks <- round((max(accum_final[][is.finite(accum_final[])]) - min(accum_final[][is.finite(accum_final[])])) / 10,2)
  }

  #crop the final accumulated dataset to the extent of the input dtm so that NA cell (e.g., cells corresponding to the sea)
  #can be excluded
  accum_final <- raster::mask(accum_final, dtm)

  #set the break values for the isolines, again excluding inf values
  levels <- seq(min(accum_final[][is.finite(accum_final[])]), max(accum_final[][is.finite(accum_final[])]), breaks)

  if (graph.out==TRUE) {

    #produce the ingredients for the hillshade raster
    #to be used in both the rendered plots
    slope <- raster::terrain(dtm, opt = "slope")
    aspect <- raster::terrain(dtm, opt = "aspect")
    hill <- raster::hillShade(slope, aspect, angle = 45, direction = 0)

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
                   col = topo.colors(255))

      #add the hillshade
      raster::plot(hill,
                   col = grey(0:100/100),
                   legend = FALSE,
                   alpha=transp,
                   add=TRUE)

      #add the contours
      raster::contour(accum_final,
                      add=TRUE,
                      levels=levels,
                      labcex=cex.breaks,
                      drawlabels = cont.lab)

      #add the origin
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

      #add the origin
      raster::plot(origin,
                   pch=20,
                   add=TRUE)
    }

  }

  #calculate and store the contours as a SpatialLinesDataFrame
  isolines <- raster::rasterToContour(accum_final, levels=levels)


  #if 'destin' is NOT NULL, calculate the least-cost path(s) from the origin to the destination(s);
  #the 'Conductance' transitional layer is used
  if(is.null(destin)==FALSE){
    #calculate the least-cost path(s)
    sPath <- gdistance::shortestPath(Conductance, sp::coordinates(origin), sp::coordinates(destin), output="SpatialLines")

    sPath.back <- NULL

    #if 'return.base' is set to TRUE
    if(return.base==TRUE){
      #create an empty list to store the LCPs back to the origin
      sPath.back <- list()
      #loop through the destinations, for each calculate the LCPs to the origin, and store those in the
      #previously created list
      for(i in 1:length(destin)) {
        sPath.back[[i]]<- gdistance::shortestPath(Conductance, sp::coordinates(destin[i,]), sp::coordinates(origin), output="SpatialLines")
      }
      #merge the individual LCPs
      sPath.back <- base::do.call(rbind, sPath.back)
    }


    if (graph.out==TRUE) {
      #plot the dtm
      raster::plot(dtm, main="Digital Terrain Model with Least-cost Path(s)",
                   sub=sub.title.lcp.plot,
                   cex.main=0.90,
                   cex.sub=0.7,
                   legend.lab="Elevation (masl)")

      #add the hillshade
      raster::plot(hill,
                   col = grey(0:100/100),
                   legend = FALSE,
                   alpha=transp,
                   add=TRUE)

      #add the origin
      raster::plot(origin, add=TRUE, pch=20)

      #add the LCPs
      raster::plot(sPath, add=TRUE)

      #if 'return.base' is set to TRUE, plot the LCPs back to the origin
      if(return.base==TRUE){
        for(i in 1:length(destin)) {
          raster::plot(sPath.back, lty=rb.lty, add=TRUE)
        }
      }

      #add the destination(s)
      raster::plot(destin,
                   add=TRUE,
                   pch=20,
                   col="red")
    }

    #calculate the length of the least-cost paths and store the values by appending them to a new variable of the sPath object
    sPath$length <- rgeos::gLength(sPath, byid=TRUE)

    #same as above if the LCPs back to the origin had been calculated
    if(return.base==TRUE){
      sPath.back$length <- rgeos::gLength(sPath.back, byid=TRUE)
    }

    #extract the cost from the accum_final to the destination location(s), appending the data to a new column
    destin$cost <- raster::extract(accum_final, destin)

    #if user select the Tobler's, the modified Tobler's, the Irmischer-Clarke's,
    #the Uriarte Gonzalez's, the Alberti's function, or other functions producing time cost...
    if (funct=="t" | funct=="tofp" | funct=="mp" | funct=="icmonp" | funct=="icmoffp" | funct=="icfonp" | funct=="icfoffp" | funct=="ug" | funct=="alb" | funct=="gkrs" | funct=="r" | funct=="ks" | funct=="ma"){
      #create a new columns in the 'destin' layer to store decimal hours turned into sessagesimal format
      if (time=="h") {
        destin$cost_hms <- as.character(chron::times(destin$cost / 24))
      } else {
        destin$cost_hms <- as.character(chron::times(destin$cost / (24*60)))
      }
    }

    if (graph.out==TRUE) {
      #if destin.lab is TRUE...
      if(destin.lab==TRUE){
        #if the dataset 'destin' contains the column "cost_hms" (as per the above step),
        #use that column to get the points' labels (this means that cost-functions that produced
        #walking time distance have been used)
        if(is.null(destin$cost_hms)==FALSE){
          raster::text(sp::coordinates(destin),
                       labels=destin$cost_hms,
                       pos = 4,
                       cex=cex.lcp.lab)
        } else {
          #otherwise (i.e., in case energy cost function has been used), use the "cost" column;
          #notice that in case of cost function returning time, the cost is still stored in the
          #"cost" column, but is the time in sessagesimal that will be used for plotting, not in decimal;
          #decimal time has been nonetheless stored
          raster::text(sp::coordinates(destin),
                       labels=round(destin$cost,2),
                       pos = 4,
                       cex=cex.lcp.lab)
        }
      }
    }

    #if export is TRUE, export the LPCs and the destinatin locations with cost as a shapefile
    if(export==TRUE){
      rgdal::writeOGR(sPath, ".", paste0("LCPs_", funct), driver="ESRI Shapefile")
      rgdal::writeOGR(destin, ".", paste0("dest_loc_w_cost_", funct), driver="ESRI Shapefile")
      #same as above if the LCPs back to the origin had been calculated
      if(return.base==TRUE){
        rgdal::writeOGR(sPath.back, ".", paste0("LCPs_back_", funct), driver="ESRI Shapefile")
      }
    }

  } else {
    sPath=NULL
    sPath.back=NULL
    dest.loc.w.cost=NULL
  }

  #if export is TRUE, export the cost-surface and the accumulated cost surface as a raster file,
  #the isolines as a shapefile
  if(export==TRUE){
    raster::writeRaster(accum_final, paste0("accum_cost_surf_", funct), format="GTiff")
    raster::writeRaster(cost.surface.to.export, paste0("cost_surf_", funct), format="GTiff")
    rgdal::writeOGR(isolines, ".", paste0("isolines_", funct), driver="ESRI Shapefile")
  }

  #if no DTM was provided (i.e., if 'studyplot' is not NULL), export the downloaded DTM as a raster file
  if(export==TRUE & is.null(studyplot)==FALSE){
    raster::writeRaster(dtm, "dtm", format="GTiff")
  }

  if (graph.out==TRUE) {
    #restore the original graphical device's settings if previously modified
    if(is.null(destin)==FALSE & oneplot==TRUE){
      par(mfrow = c(1,1))
    }
  }

  #restore the advice for error messages on the R console, which has been deactivated
  #at the beginning of function
  options(warn = 1)

  #store the results in a list
  results <- list("dtm"=dtm,
                  "cost.surface"=cost.surface.to.export,
                  "accumulated.cost.raster"=accum_final,
                  "isolines" = isolines,
                  "LCPs"=sPath,
                  "LCPs.back"=sPath.back,
                  "dest.loc.w.cost"=destin,
                  "conductance"=Conductance)
}
