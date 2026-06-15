#' The movecost cost-function registry
#'
#' \code{mc_cost_functions()} returns a data frame describing every cost function
#' implemented in the package. The registry is the single internal source from
#' which all analysis functions draw their cost models: each function is defined
#' once, here, and reused everywhere (in movecost <= 2.x the same definitions were
#' repeated across functions; the registry eliminates that duplication).\cr
#'
#' Each cost function maps the terrain slope between two adjacent cells
#' (expressed as rise over run, signed: positive uphill, negative downhill, and
#' computed in the direction of movement) to either a walking speed or the
#' reciprocal of a per-metre cost. In what follows, \strong{x} stands for slope as
#' rise/run calculated between adjacent cells in the direction of travel.\cr
#'
#' \strong{Functions expressing cost as walking time (speed in km/h)}:\cr
#'
#' \strong{Tobler's hiking function (on-path)} (\code{"t"}):\cr
#' \eqn{(6 * exp(-3.5 * abs(x + 0.05))) * (1/N)}\cr
#' Tobler W. (1993), Three Presentations on Geographical Analysis and Modeling,
#' NCGIA Technical Report 93-1.\cr
#'
#' \strong{Tobler's hiking function (off-path)} (\code{"tofp"}):\cr
#' \eqn{(6 * exp(-3.5 * abs(x + 0.05))) * 0.6}\cr
#' as per Tobler's indication, the off-path walking speed is reduced by 0.6.\cr
#'
#' \strong{Marquez-Perez et al.'s modified Tobler hiking function} (\code{"mp"}):\cr
#' \eqn{(4.8 * exp(-5.3 * abs((x * 0.7) + 0.03))) * (1/N)}\cr
#' Marquez-Perez J., Vallejo-Villalta I., Alvarez-Francoso J.I. (2017),
#' "Estimated travel time for walking trails in natural areas",
#' Geografisk Tidsskrift-Danish Journal of Geography, 117:1, 53-62,
#' \doi{10.1080/00167223.2017.1316212}.\cr
#'
#' \strong{Irmischer-Clarke's hiking function (male, on-path)} (\code{"icmonp"}):\cr
#' \eqn{((0.11 + exp(-(abs(x)*100 + 5)^2 / (2 * 30^2))) * 3.6) * (1/N)}\cr
#' Irmischer I.J., Clarke K.C. (2018), Measuring and modeling the speed of human
#' navigation, Cartography and Geographic Information Science, 45(2), 177-186,
#' \doi{10.1080/15230406.2017.1292150}. The original functions express speed in
#' m/s; they are multiplied by 3.6 to obtain km/h for consistency with the other
#' Tobler-related cost functions; slope is in percent.\cr
#'
#' \strong{Irmischer-Clarke's hiking function (male, off-path)} (\code{"icmoffp"}):\cr
#' \eqn{(0.11 + 0.67 * exp(-(abs(x)*100 + 2)^2 / (2 * 30^2))) * 3.6}\cr
#' the Gaussian denominator is \eqn{2*30^2} (=1800), the 2*sigma^2 form of
#' equation 4 of Irmischer & Clarke (2018), shared by the on-path (equation 2)
#' and female (equation 5) variants; unchanged from movecost <= 2.x.\cr
#'
#' \strong{Irmischer-Clarke's hiking function (female, on-path)} (\code{"icfonp"}):\cr
#' \eqn{((0.95 * (0.11 + exp(-(abs(x) * 100 + 5)^2/(2 * 30^2)))) * 3.6) * (1/N)}\cr
#'
#' \strong{Irmischer-Clarke's hiking function (female, off-path)} (\code{"icfoffp"}):\cr
#' \eqn{(0.95 * (0.11 + 0.67 * exp(-(abs(x) * 100 + 2)^2/(2 * 30^2)))) * 3.6}\cr
#'
#' \strong{Uriarte Gonzalez's walking-time cost function} (\code{"ug"}):\cr
#' \eqn{1 / ((0.0277 * (abs(x)*100) + 0.6115) * N)}\cr
#' Chapa Brunet T., Garcia J., Mayoral Herrera V., Uriarte Gonzalez A. (2008),
#' GIS landscape models for the study of preindustrial settlement patterns in
#' Mediterranean areas, in Geoinformation Technologies for Geo-Cultural
#' Landscapes (pp. 255-273), CRC Press, \doi{10.1201/9780203881613.ch12}.
#' The function originally expresses walking time in seconds; its reciprocal is
#' used internally so that the accumulated value corresponds to time. Unlike in
#' the original formulation, the pixel resolution is not embedded in the formula
#' because the graph machinery accounts for the actual distance between cell
#' centres.\cr
#'
#' \strong{Marin Arroyo's walking-time cost function} (\code{"ma"}):\cr
#' \eqn{ifelse(x < 0, 1 / ((0.6 * ((abs(x)*100)/23+1))*N), 1 / ((0.6 * ((abs(x)*100)/11+1))*N))}\cr
#' Marin Arroyo A.B. (2009), The use of optimal foraging theory to estimate Late
#' Glacial site catchment areas from a central place: the case of eastern
#' Cantabria, Spain, Journal of Anthropological Archaeology 28, 27-36.
#' \strong{Note}: in movecost <= 2.x the downhill branch could never be selected
#' because the test was applied to the absolute slope; version 3.0 evaluates the
#' signed slope, so downhill movement now correctly uses the gentler divisor
#' (23), as in the original publication. Results for downhill movement therefore
#' differ (correctly) from movecost <= 2.x.\cr
#'
#' \strong{Alberti's Tobler hiking function modified for pastoral foraging
#' excursions} (\code{"alb"}):\cr
#' \eqn{(6 * exp(-3.5 * abs(x + 0.05))) * 0.25}\cr
#' Alberti G. (2019), Locating potential pastoral foraging routes in Malta
#' through the use of a Geographic Information System. The Tobler function is
#' rescaled by 0.25, i.e. the ratio between the average flock speed (1.5 km/h)
#' and the maximum human walking speed (about 6.0 km/h) on favourable slopes.\cr
#'
#' \strong{Garmy, Kaddouri, Rozenblat, and Schneider's hiking function}
#' (\code{"gkrs"}):\cr
#' \eqn{(4 * exp(-0.008 * ((atan(abs(x))*180/pi)^2))) * (1/N)}\cr
#' slope in degrees; see Herzog I. (2020), Spatial Analysis Based on Cost
#' Functions, in Gillings M., Haciguzeller P., Lock G. (eds), "Archaeological
#' Spatial Analysis. A Methodological Guide", Routledge: New York, 333-358.\cr
#'
#' \strong{Rees' hiking function} (\code{"r"}):\cr
#' \eqn{((1 / (0.75 + 0.09 * abs(x) + 14.6 * (abs(x))^2)) * 3.6) * (1/N)}\cr
#' Rees W.G. (2004), Least-cost paths in mountainous terrain,
#' Computers & Geosciences, 30(3), 203-209. Speed transposed from m/s to km/h.\cr
#'
#' \strong{Kondo-Seino's modified Tobler hiking function} (\code{"ks"}):\cr
#' \eqn{ifelse(x >= -0.07, (5.1 * exp(-2.25 * abs(x + 0.07))) * (1/N), (5.1 * exp(-1.5 * abs(x + 0.07))) * (1/N))}\cr
#' Kondo Y., Seino Y. (2010), GPS-aided Walking Experiments and Data-driven
#' Travel Cost Modeling on the Historical Road of Nakasendo-Kisoji (Central
#' Highland Japan), in Frischer B., Webb Crawford J., Koller D. (eds), Making
#' History Interactive, CAA Proceedings of the 37th International Conference
#' (BAR International Series S2079), Archaeopress, Oxford, 158-165.
#' \strong{Note}: in movecost <= 2.x the steep-downhill branch could never be
#' selected because the test was applied to the absolute slope; version 3.0
#' evaluates the signed slope as intended by Kondo and Seino, so walking speed on
#' downhill gradients steeper than -7\% now correctly uses the gentler decay
#' coefficient (1.5). Results on such gradients therefore differ (correctly)
#' from movecost <= 2.x.\cr
#'
#' \strong{Tripcevich's hiking function} (\code{"trp"}):\cr
#' \eqn{((4.028*46^2)/(((atan(abs(x))*180/pi)+4.127)^2+46^2))*(1/N)}\cr
#' Tripcevich N. (2008), Estimating Llama caravan travel speeds:
#' ethno-archaeological fieldwork with a Peruvian salt caravan.\cr
#'
#' \strong{Function for wheeled vehicles}:\cr
#'
#' \strong{Wheeled-vehicle critical slope cost function} (\code{"wcs"}):\cr
#' \eqn{1 / ((1 + ((abs(x)*100) / sl.crit)^2) * N)}\cr
#' where \code{sl.crit} (critical slope, in percent) is "the transition where
#' switchbacks become more effective than direct uphill or downhill paths",
#' typically in the range 8-16; Herzog I. (2016), Potential and Limits of
#' Optimal Path Analysis, in Bevan A., Lake M. (eds), Computational Approaches
#' to Archaeological Spaces (pp. 179-211), Routledge.\cr
#'
#' \strong{Functions expressing abstract cost}:\cr
#'
#' \strong{Relative energetic expenditure} (\code{"ree"}):\cr
#' \eqn{1 / ((tan((atan(abs(x))*180/pi)*pi/180) / tan(1*pi/180)) * N)}\cr
#' Conolly J., Lake M. (2006), Geographic Information Systems in Archaeology,
#' Cambridge University Press, p. 220.\cr
#'
#' \strong{Bellavia's cost function} (\code{"b"}):\cr
#' \eqn{1 / (((atan(abs(x))*180/pi)+1) * N)}\cr
#' see Herzog I. (2020), cited above.\cr
#'
#' \strong{Eastman's cost function} (\code{"e"}):\cr
#' \eqn{1 / ((0.031*(atan(abs(x))*180/pi)^2 - 0.025*(atan(abs(x))*180/pi) + 1) * N)}\cr
#' Vaissie E. (2021), Mobility of Paleolithic Populations: Biomechanical
#' Considerations and Spatiotemporal Modelling, PaleoAnthropology 2021(1),
#' 120-144 (with reference to Eastman 1999).\cr
#'
#' \strong{Functions expressing cost as metabolic energy expenditure}:\cr
#'
#' \strong{Pandolf et al.'s cost function (in Watts)} (\code{"p"}):\cr
#' \eqn{1 / (1.5 * W + 2.0 * (W + L) * (L / W)^2 + N * (W + L) * (1.5 * (V^2) + 0.35 * V * (abs(x)*100)))}\cr
#' where W is the walker's body weight (kg), L the carried load (kg), V the
#' velocity in m/s, N the terrain coefficient (the \eqn{\eta} of the original
#' publication, multiplying the velocity-dependent term once). If V is set to
#' 0, it is worked out internally from the Tobler on-path function and varies
#' with the slope. \strong{Note}: movecost <= 2.x additionally multiplied the
#' whole expression by N, applying the terrain factor twice; fixed in version
#' 3.0 (see NEWS.md). With the default N = 1 the results coincide.
#' Pandolf K.B., Givoni B., Goldman R.F. (1977), Predicting energy expenditure
#' with loads while standing or walking very slowly, Journal of Applied
#' Physiology, 43(4), 577-581, \doi{10.1152/jappl.1977.43.4.577}.\cr
#'
#' \strong{Pandolf et al.'s cost function with downhill correction (in Watts)}
#' (\code{"pcf"}):\cr
#' \eqn{ifelse(x >= 0, 1 / (M), 1 / (M - CF))} where M is the Pandolf
#' et al. metabolic rate above and \eqn{CF = N * (G*(W+L)*V/3.5) - ((W+L)*(G+6)^2/W) + (25-V^2)}
#' with G the slope magnitude in percent; Yokota M., Berglund L.G., Santee
#' W.R., Buller M.J., Hoyt R.W. (2004), U.S. Army Research Institute of
#' Environmental Medicine Technical Report T04-09.
#' \strong{Note}: in movecost <= 2.x the test selecting the downhill branch was
#' applied to the absolute slope, so the downhill correction could never be
#' selected; version 3.0 evaluates the signed slope, applying the Yokota et
#' al. correction to downhill movement as intended. Results for downhill
#' movement therefore differ (correctly) from movecost <= 2.x.\cr
#'
#' \strong{Minetti et al.'s metabolic cost function (J/(kg*m))} (\code{"m"}):\cr
#' \eqn{1 / (((280.5 * abs(x)^5) - (58.7 * abs(x)^4) - (76.8 * abs(x)^3) + (51.9 * abs(x)^2) + (19.6 * abs(x)) + 2.5) * N)}\cr
#' Minetti A.E., Moia C., Roi G.S., Susta D., Ferretti G. (2002), Energy cost of
#' walking and running at extreme uphill and downhill slopes, Journal of Applied
#' Physiology 93, 1039-1046. Valid for slopes in the range -0.5/0.5; outside
#' this range its output becomes counterintuitive (see Herzog 2013), in which
#' case the Herzog polynomial approximation (\code{"hrz"}) is preferable.\cr
#'
#' \strong{Herzog's metabolic cost function (J/(kg*m))} (\code{"hrz"}):\cr
#' \eqn{1 / (((1337.8 * abs(x)^6) + (278.19 * abs(x)^5) - (517.39 * abs(x)^4) - (78.199 * abs(x)^3) + (93.419 * abs(x)^2) + (19.825 * abs(x)) + 1.64) * N)}\cr
#' Herzog I. (2016), cited above; 6th-degree polynomial approximation of
#' Minetti et al.'s function.\cr
#'
#' \strong{Van Leusen's metabolic cost function (in Watts)} (\code{"vl"}):\cr
#' \eqn{1 / (1.5 * W + 2.0 * (W + L) * (L / W)^2 + N * (W + L) * (1.5 * (V^2) + 0.35 * V * ((abs(x)*100) + 10)))}\cr
#' (terrain factor applied once, as in the Pandolf et al. equation it
#' modifies; movecost <= 2.x applied it twice, see the note under the Pandolf
#' entry above.)
#' Van Leusen P.M. (2002), Pattern to process; with the amendments reported by
#' Herzog I. (2013), Least-cost Paths - Some Methodological Issues,
#' Internet Archaeology 36.\cr
#'
#' \strong{Llobera-Sluckin's metabolic cost function (kJ/m)} (\code{"ls"}):\cr
#' \eqn{1 / ((2.635 + (17.37 * abs(x)) + (42.37 * abs(x)^2) - (21.43 * abs(x)^3) + (14.93 * abs(x)^4)) * N)}\cr
#' Llobera M., Sluckin T.J. (2007), Zigzagging: Theoretical insights on climbing
#' strategies, Journal of Theoretical Biology 249, 206-217.\cr
#'
#' \strong{Ardigo et al.'s metabolic cost function (J/(kg*m))} (\code{"a"}):\cr
#' \eqn{1 / ((1.866 * exp(4.911*abs(x)) * V^2 - 3.773 * exp(3.416*abs(x)) * V + (45.71*abs(x)^2 + 18.90*abs(x)) + 4.456) * N)}\cr
#' Ardigo L.P., Saibene F., Minetti A.E. (2003), The optimal locomotion on
#' gradients: walking, running or cycling?, Eur J Appl Physiol 90, 365-371.
#' If V is set to 0, it is worked out internally from the Tobler on-path
#' function.\cr
#'
#' \strong{Hare's metabolic cost function (cal/km)} (\code{"h"}):\cr
#' \eqn{1 / ((48 + 30 / (6 * exp(-3.5 * abs(x + 0.05)))) * N)}\cr
#' Hare T.S. (2004), Using Measures of Cost Distance in the Estimation of Polity
#' Boundaries in the Post Classic Yautepec valley, Mexico, Journal of
#' Archaeological Science 31.\cr
#'
#' \strong{Terrain factor (N)}:\cr
#' virtually all the implemented cost functions (with few exceptions, namely the
#' off-path functions and Alberti's, which natively embed one) can take into
#' account a terrain factor (parameter \code{N}; 1 by default) representing the
#' easiness/difficulty of moving on different terrain types.\cr
#'
#' The following reference list of terrain factors is based on the data
#' collected in Herzog, I. (2020), Spatial Analysis Based on Cost Functions, in
#' Gillings M., Haciguzeller P., Lock G. (eds), "Archaeological Spatial
#' Analysis. A Methodological Guide", Routledge: New York, 340 (with previous
#' references). It is divided into two sections: section (a) reports the terrain
#' factors to be used for cost functions measuring time; section (b) for
#' functions measuring cost other than time.\cr
#'
#' (a) for time-based cost functions:
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
#' (b) for cost functions measuring cost other than time:
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
#' \strong{Why speed-type and reciprocal-type functions are treated
#' differently}:\cr
#' internally the package works with conductance (the reciprocal of cost), so
#' walking-speed functions are used as they are (the accumulated value being the
#' reciprocal of speed, i.e. pace, integrated over distance, which yields time),
#' while the other functions are reciprocated so that the accumulated value
#' corresponds to the original cost (see Nakoinz-Knitter (2016), "Modelling
#' Human Behaviour in Landscapes", Springer, p. 183).
#'
#' @return A data frame with one row per cost function and columns:
#' \itemize{
#'   \item \code{code}: the identifier to pass as the \code{funct} parameter;
#'   \item \code{name}: extended name of the cost function;
#'   \item \code{type}: \code{"walking time"}, \code{"energy"}, \code{"abstract"},
#'     or \code{"wheeled vehicle"};
#'   \item \code{cost_unit}: unit of the accumulated cost
#'     (before any user-selected time conversion);
#'   \item \code{uses_N}, \code{uses_WLV}, \code{uses_slcrit}: whether the
#'     function responds to the terrain factor, to the W/L/V parameters, and to
#'     the critical-slope parameter respectively.
#' }
#'
#' @examples
#' # list all the implemented cost functions
#' mc_cost_functions()
#'
#' @seealso \code{\link{mc_surface}}
#' @export
mc_cost_functions <- function() {
  reg <- .cf_registry()
  out <- do.call(rbind, lapply(reg, function(e) {
    data.frame(code = e$code, name = e$name, type = e$type,
               cost_unit = e$unit, uses_N = e$uses_N,
               uses_WLV = e$uses_WLV, uses_slcrit = e$uses_slcrit,
               stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out
}

#' Internal: full registry of cost-function definitions
#'
#' Each entry provides a factory ('make') that, given the user parameters,
#' returns a vectorised function of the SIGNED slope (rise over run, in the
#' direction of travel). 'group' is either "speed" (the function returns a
#' speed in km/h, converted internally to m/s) or "recip" (the function
#' returns the reciprocal of a per-metre cost).
#'
#' @keywords internal
#' @noRd
.cf_registry <- function() {
  # Tobler on-path speed (km/h) used to derive V internally when V = 0
  tobler_ms <- function(x) (6 * exp(-3.5 * abs(x + 0.05))) * 0.278

  list(
    t = list(
      code = "t", name = "Tobler's hiking function (on-path)",
      type = "walking time", unit = "hours (or minutes)", group = "speed",
      time_output = TRUE, megawatt = FALSE,
      uses_N = TRUE, uses_WLV = FALSE, uses_slcrit = FALSE,
      label = function(p) paste0("Tobler's on-path hiking function\nterrain factor N=", p$N),
      make = function(p) function(x) (6 * exp(-3.5 * abs(x + 0.05))) * (1 / p$N)
    ),
    tofp = list(
      code = "tofp", name = "Tobler's hiking function (off-path)",
      type = "walking time", unit = "hours (or minutes)", group = "speed",
      time_output = TRUE, megawatt = FALSE,
      uses_N = FALSE, uses_WLV = FALSE, uses_slcrit = FALSE,
      label = function(p) "Tobler's off-path hiking function",
      make = function(p) function(x) (6 * exp(-3.5 * abs(x + 0.05))) * 0.6
    ),
    mp = list(
      code = "mp", name = "Marquez-Perez et al.'s modified Tobler hiking function",
      type = "walking time", unit = "hours (or minutes)", group = "speed",
      time_output = TRUE, megawatt = FALSE,
      uses_N = TRUE, uses_WLV = FALSE, uses_slcrit = FALSE,
      label = function(p) paste0("Marquez-Perez et al.'s modified Tobler hiking function\nterrain factor N=", p$N),
      make = function(p) function(x) (4.8 * exp(-5.3 * abs((x * 0.7) + 0.03))) * (1 / p$N)
    ),
    icmonp = list(
      code = "icmonp", name = "Irmischer-Clarke's hiking function (male, on-path)",
      type = "walking time", unit = "hours (or minutes)", group = "speed",
      time_output = TRUE, megawatt = FALSE,
      uses_N = TRUE, uses_WLV = FALSE, uses_slcrit = FALSE,
      label = function(p) paste0("(male, on-path) Irmischer-Clarke's hiking function\nterrain factor N=", p$N),
      make = function(p) function(x) ((0.11 + exp(-(abs(x) * 100 + 5)^2 / (2 * 30^2))) * 3.6) * (1 / p$N)
    ),
    icmoffp = list(
      code = "icmoffp", name = "Irmischer-Clarke's hiking function (male, off-path)",
      type = "walking time", unit = "hours (or minutes)", group = "speed",
      time_output = TRUE, megawatt = FALSE,
      uses_N = FALSE, uses_WLV = FALSE, uses_slcrit = FALSE,
      label = function(p) "(male, off-path) Irmischer-Clarke's hiking function",
      # Gaussian denominator 2*30^2 (=1800), the 2*sigma^2 form of eq. 4 of
      # Irmischer & Clarke (2018); identical to movecost <= 2.x.
      make = function(p) function(x) (0.11 + 0.67 * exp(-(abs(x) * 100 + 2)^2 / (2 * 30^2))) * 3.6
    ),
    icfonp = list(
      code = "icfonp", name = "Irmischer-Clarke's hiking function (female, on-path)",
      type = "walking time", unit = "hours (or minutes)", group = "speed",
      time_output = TRUE, megawatt = FALSE,
      uses_N = TRUE, uses_WLV = FALSE, uses_slcrit = FALSE,
      label = function(p) paste0("(female, on-path) Irmischer-Clarke's hiking function\nterrain factor N=", p$N),
      make = function(p) function(x) ((0.95 * (0.11 + exp(-(abs(x) * 100 + 5)^2 / (2 * 30^2)))) * 3.6) * (1 / p$N)
    ),
    icfoffp = list(
      code = "icfoffp", name = "Irmischer-Clarke's hiking function (female, off-path)",
      type = "walking time", unit = "hours (or minutes)", group = "speed",
      time_output = TRUE, megawatt = FALSE,
      uses_N = FALSE, uses_WLV = FALSE, uses_slcrit = FALSE,
      label = function(p) "(female, off-path) Irmischer-Clarke's hiking function",
      make = function(p) function(x) (0.95 * (0.11 + 0.67 * exp(-(abs(x) * 100 + 2)^2 / (2 * 30^2)))) * 3.6
    ),
    ug = list(
      code = "ug", name = "Uriarte Gonzalez's walking-time cost function",
      type = "walking time", unit = "hours (or minutes)", group = "recip",
      time_output = TRUE, megawatt = FALSE,
      uses_N = TRUE, uses_WLV = FALSE, uses_slcrit = FALSE,
      label = function(p) paste0("Uriarte Gonzalez's walking-time cost function\nterrain factor N=", p$N),
      make = function(p) function(x) 1 / ((0.0277 * (abs(x) * 100) + 0.6115) * p$N)
    ),
    ma = list(
      code = "ma", name = "Marin Arroyo's walking-time cost function",
      type = "walking time", unit = "hours (or minutes)", group = "recip",
      time_output = TRUE, megawatt = FALSE,
      uses_N = TRUE, uses_WLV = FALSE, uses_slcrit = FALSE,
      label = function(p) paste0("Marin Arroyo's walking-time cost function\nterrain factor N=", p$N),
      make = function(p) function(x) {
        ifelse(x < 0,
               1 / ((0.6 * ((abs(x) * 100) / 23 + 1)) * p$N),
               1 / ((0.6 * ((abs(x) * 100) / 11 + 1)) * p$N))
      }
    ),
    alb = list(
      code = "alb", name = "Alberti's Tobler hiking function modified for pastoral foraging excursions",
      type = "walking time", unit = "hours (or minutes)", group = "speed",
      time_output = TRUE, megawatt = FALSE,
      uses_N = FALSE, uses_WLV = FALSE, uses_slcrit = FALSE,
      label = function(p) "Alberti's modified Tobler hiking function",
      make = function(p) function(x) (6 * exp(-3.5 * abs(x + 0.05))) * 0.25
    ),
    gkrs = list(
      code = "gkrs", name = "Garmy, Kaddouri, Rozenblat, and Schneider's hiking function",
      type = "walking time", unit = "hours (or minutes)", group = "speed",
      time_output = TRUE, megawatt = FALSE,
      uses_N = TRUE, uses_WLV = FALSE, uses_slcrit = FALSE,
      label = function(p) paste0("Garmy et al.'s hiking function\nterrain factor N=", p$N),
      make = function(p) function(x) (4 * exp(-0.008 * ((atan(abs(x)) * 180 / pi)^2))) * (1 / p$N)
    ),
    r = list(
      code = "r", name = "Rees' hiking function",
      type = "walking time", unit = "hours (or minutes)", group = "speed",
      time_output = TRUE, megawatt = FALSE,
      uses_N = TRUE, uses_WLV = FALSE, uses_slcrit = FALSE,
      label = function(p) paste0("Rees' hiking function\nterrain factor N=", p$N),
      make = function(p) function(x) ((1 / (0.75 + 0.09 * abs(x) + 14.6 * (abs(x))^2)) * 3.6) * (1 / p$N)
    ),
    ks = list(
      code = "ks", name = "Kondo-Seino's modified Tobler hiking function",
      type = "walking time", unit = "hours (or minutes)", group = "speed",
      time_output = TRUE, megawatt = FALSE,
      uses_N = TRUE, uses_WLV = FALSE, uses_slcrit = FALSE,
      label = function(p) paste0("Kondo-Seino's hiking function\nterrain factor N=", p$N),
      make = function(p) function(x) {
        ifelse(x >= -0.07,
               (5.1 * exp(-2.25 * abs(x + 0.07))) * (1 / p$N),
               (5.1 * exp(-1.5  * abs(x + 0.07))) * (1 / p$N))
      }
    ),
    trp = list(
      code = "trp", name = "Tripcevich's hiking function",
      type = "walking time", unit = "hours (or minutes)", group = "speed",
      time_output = TRUE, megawatt = FALSE,
      uses_N = TRUE, uses_WLV = FALSE, uses_slcrit = FALSE,
      label = function(p) paste0("Tripcevich's hiking function\nterrain factor N=", p$N),
      make = function(p) function(x) ((4.028 * 46^2) / (((atan(abs(x)) * 180 / pi) + 4.127)^2 + 46^2)) * (1 / p$N)
    ),
    wcs = list(
      code = "wcs", name = "Wheeled-vehicle critical slope cost function",
      type = "wheeled vehicle", unit = "abstract cost", group = "recip",
      time_output = FALSE, megawatt = FALSE,
      uses_N = TRUE, uses_WLV = FALSE, uses_slcrit = TRUE,
      label = function(p) paste0("wheeled-vehicle critical slope cost function\ncritical slope set to ", p$sl.crit, " percent; terrain factor N=", p$N),
      make = function(p) function(x) 1 / ((1 + ((abs(x) * 100) / p$sl.crit)^2) * p$N)
    ),
    ree = list(
      code = "ree", name = "Relative energetic expenditure cost function",
      type = "abstract", unit = "abstract cost", group = "recip",
      time_output = FALSE, megawatt = FALSE,
      uses_N = TRUE, uses_WLV = FALSE, uses_slcrit = FALSE,
      label = function(p) paste0("relative energetic expenditure cost function\nterrain factor N=", p$N),
      make = function(p) function(x) 1 / ((tan((atan(abs(x)) * 180 / pi) * pi / 180) / tan(1 * pi / 180)) * p$N)
    ),
    b = list(
      code = "b", name = "Bellavia's cost function",
      type = "abstract", unit = "abstract cost", group = "recip",
      time_output = FALSE, megawatt = FALSE,
      uses_N = TRUE, uses_WLV = FALSE, uses_slcrit = FALSE,
      label = function(p) paste0("Bellavia's cost function\nterrain factor N=", p$N),
      make = function(p) function(x) 1 / (((atan(abs(x)) * 180 / pi) + 1) * p$N)
    ),
    e = list(
      code = "e", name = "Eastman's cost function",
      type = "abstract", unit = "abstract cost", group = "recip",
      time_output = FALSE, megawatt = FALSE,
      uses_N = TRUE, uses_WLV = FALSE, uses_slcrit = FALSE,
      label = function(p) paste0("Eastman's cost function\nterrain factor N=", p$N),
      make = function(p) function(x) 1 / ((0.031 * (atan(abs(x)) * 180 / pi)^2 - 0.025 * (atan(abs(x)) * 180 / pi) + 1) * p$N)
    ),
    p = list(
      code = "p", name = "Pandolf et al.'s metabolic energy expenditure cost function",
      type = "energy", unit = "Megawatts", group = "recip",
      time_output = FALSE, megawatt = TRUE,
      uses_N = TRUE, uses_WLV = TRUE, uses_slcrit = FALSE,
      label = function(p) paste0("Pandolf et al.'s metabolic energy expenditure cost function\ncost in Megawatts; parameters: W: ", p$W, "; L: ", p$L, "; N: ", p$N, "; V: ", if (p$V == 0) "Tobler-derived" else p$V),
      # NOTE: the terrain factor is applied ONCE, as the coefficient eta of
      # Pandolf et al. (1977), multiplying the velocity-dependent term.
      # movecost <= 2.x additionally multiplied the whole expression by N,
      # counting terrain twice; fixed in 3.0 (see NEWS.md). With the default
      # N = 1 the two formulations coincide.
      make = function(p) {
        if (p$V == 0) {
          function(x) {
            V_int <- tobler_ms(x)
            1 / (1.5 * p$W + 2.0 * (p$W + p$L) * (p$L / p$W)^2 +
                   p$N * (p$W + p$L) * (1.5 * V_int^2 + 0.35 * V_int * (abs(x) * 100)))
          }
        } else {
          function(x) 1 / (1.5 * p$W + 2.0 * (p$W + p$L) * (p$L / p$W)^2 +
                             p$N * (p$W + p$L) * (1.5 * (p$V^2) + 0.35 * p$V * (abs(x) * 100)))
        }
      }
    ),
    pcf = list(
      code = "pcf", name = "Pandolf et al.'s cost function with downhill correction factor",
      type = "energy", unit = "Megawatts", group = "recip",
      time_output = FALSE, megawatt = TRUE,
      uses_N = TRUE, uses_WLV = TRUE, uses_slcrit = FALSE,
      label = function(p) paste0("Pandolf et al.'s cost function with downhill correction\ncost in Megawatts; parameters: W: ", p$W, "; L: ", p$L, "; N: ", p$N, "; V: ", if (p$V == 0) "Tobler-derived" else p$V),
      make = function(p) {
        function(x) {
          V_use <- if (p$V == 0) tobler_ms(x) else p$V
          G <- abs(x) * 100
          M_up <- 1.5 * p$W + 2.0 * (p$W + p$L) * (p$L / p$W)^2 +
            p$N * (p$W + p$L) * (1.5 * V_use^2 + 0.35 * V_use * G)
          M_dn <- M_up - (p$N * (G * (p$W + p$L) * V_use / 3.5) -
                            ((p$W + p$L) * (G + 6)^2 / p$W) + (25 - V_use^2))
          # signed-slope test: the Yokota et al. correction applies to
          # DOWNHILL movement only (in movecost <= 2.x the test was applied
          # to the absolute slope, so the correction was never selected on
          # downhill gradients; see NEWS.md)
          ifelse(x >= 0, 1 / M_up, 1 / M_dn)
        }
      }
    ),
    m = list(
      code = "m", name = "Minetti et al.'s metabolic energy expenditure cost function",
      type = "energy", unit = "J/(kg*m)", group = "recip",
      time_output = FALSE, megawatt = FALSE,
      uses_N = TRUE, uses_WLV = FALSE, uses_slcrit = FALSE,
      label = function(p) paste0("Minetti et al.'s metabolic cost function\ncost in J/(kg*m); terrain factor N=", p$N),
      make = function(p) function(x) {
        a <- abs(x)
        1 / (((280.5 * a^5) - (58.7 * a^4) - (76.8 * a^3) + (51.9 * a^2) + (19.6 * a) + 2.5) * p$N)
      }
    ),
    hrz = list(
      code = "hrz", name = "Herzog's metabolic cost function",
      type = "energy", unit = "J/(kg*m)", group = "recip",
      time_output = FALSE, megawatt = FALSE,
      uses_N = TRUE, uses_WLV = FALSE, uses_slcrit = FALSE,
      label = function(p) paste0("Herzog's metabolic cost function\ncost in J/(kg*m); terrain factor N=", p$N),
      make = function(p) function(x) {
        a <- abs(x)
        1 / (((1337.8 * a^6) + (278.19 * a^5) - (517.39 * a^4) -
                (78.199 * a^3) + (93.419 * a^2) + (19.825 * a) + 1.64) * p$N)
      }
    ),
    vl = list(
      code = "vl", name = "Van Leusen's metabolic energy expenditure cost function",
      type = "energy", unit = "Megawatts", group = "recip",
      time_output = FALSE, megawatt = TRUE,
      uses_N = TRUE, uses_WLV = TRUE, uses_slcrit = FALSE,
      label = function(p) paste0("Van Leusen's metabolic energy expenditure cost function\ncost in Megawatts; parameters: W: ", p$W, "; L: ", p$L, "; N: ", p$N, "; V: ", if (p$V == 0) "Tobler-derived" else p$V),
      # NOTE: terrain factor applied once (Pandolf's eta), as in the
      # publication; movecost <= 2.x multiplied the whole expression by N a
      # second time (fixed in 3.0, see NEWS.md).
      make = function(p) {
        if (p$V == 0) {
          function(x) {
            V_int <- tobler_ms(x)
            1 / (1.5 * p$W + 2.0 * (p$W + p$L) * (p$L / p$W)^2 +
                   p$N * (p$W + p$L) * (1.5 * V_int^2 + 0.35 * V_int * ((abs(x) * 100) + 10)))
          }
        } else {
          function(x) 1 / (1.5 * p$W + 2.0 * (p$W + p$L) * (p$L / p$W)^2 +
                             p$N * (p$W + p$L) * (1.5 * (p$V^2) + 0.35 * p$V * ((abs(x) * 100) + 10)))
        }
      }
    ),
    ls = list(
      code = "ls", name = "Llobera-Sluckin's metabolic energy expenditure cost function",
      type = "energy", unit = "kJ/m", group = "recip",
      time_output = FALSE, megawatt = FALSE,
      uses_N = TRUE, uses_WLV = FALSE, uses_slcrit = FALSE,
      label = function(p) paste0("Llobera-Sluckin's metabolic cost function\ncost in kJ/m; terrain factor N=", p$N),
      make = function(p) function(x) {
        a <- abs(x)
        1 / ((2.635 + (17.37 * a) + (42.37 * a^2) - (21.43 * a^3) + (14.93 * a^4)) * p$N)
      }
    ),
    a = list(
      code = "a", name = "Ardigo et al.'s metabolic energy expenditure cost function",
      type = "energy", unit = "J/(kg*m)", group = "recip",
      time_output = FALSE, megawatt = FALSE,
      uses_N = TRUE, uses_WLV = TRUE, uses_slcrit = FALSE,
      label = function(p) paste0("Ardigo et al.'s metabolic cost function\ncost in J/(kg*m); parameters: N: ", p$N, "; V: ", if (p$V == 0) "Tobler-derived" else p$V),
      make = function(p) {
        if (p$V == 0) {
          function(x) {
            V_int <- tobler_ms(x)
            1 / ((1.866 * exp(4.911 * abs(x)) * V_int^2 -
                    3.773 * exp(3.416 * abs(x)) * V_int +
                    (45.71 * abs(x)^2 + 18.90 * abs(x)) + 4.456) * p$N)
          }
        } else {
          function(x) 1 / ((1.866 * exp(4.911 * abs(x)) * p$V^2 -
                              3.773 * exp(3.416 * abs(x)) * p$V +
                              (45.71 * abs(x)^2 + 18.90 * abs(x)) + 4.456) * p$N)
        }
      }
    ),
    h = list(
      code = "h", name = "Hare's metabolic energy expenditure cost function",
      type = "energy", unit = "cal/km", group = "recip",
      time_output = FALSE, megawatt = FALSE,
      uses_N = TRUE, uses_WLV = FALSE, uses_slcrit = FALSE,
      label = function(p) paste0("Hare's metabolic cost function\ncost in cal/km; terrain factor N=", p$N),
      make = function(p) function(x) 1 / ((48 + 30 / (6 * exp(-3.5 * abs(x + 0.05)))) * p$N)
    )
  )
}

#' Internal: fetch one registry entry, with validation
#' @keywords internal
#' @noRd
.cf_get <- function(funct) {
  reg <- .cf_registry()
  if (!funct %in% names(reg)) {
    stop("Unknown cost function '", funct, "'.\nAvailable codes: ",
         paste(names(reg), collapse = ", "),
         "\nSee mc_cost_functions() for details.", call. = FALSE)
  }
  reg[[funct]]
}
