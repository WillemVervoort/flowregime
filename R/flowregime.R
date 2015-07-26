#' Environmental Flow Regime Components
#'
#' A package for computing environmental flow regime components
#' from hydrologic time series. See the vignette to get started.
#' @name flowregime-package
#' @docType package
#' @import zoo xts
#' @importFrom EcoHydRology BaseflowSeparation
NULL

#' Missouri River Flows at Sioux City, IA
#' 
#' Missouri River daily discharge recorded at USGS gauge 06486000 in 
#'   Sioux City, IA from 2004-01-01 to 2015-05-25.
#' @docType data
#' @keywords datasets
#' @name siouxcity
#' @usage data(siouxcity)
#' @format An \code{xts} object.
NULL

#' Separate Flow Data
#'
#' Separate flow data into baseflow and quickflow. 
#' 
#' @details This function is basically a \code{xts} wrapper for
#'   \code{EcoHydRology::BaseflowSeparation}.
#' @seealso \code{\link[EcoHydRology]{BaseflowSeparation}}. 
#' @param ts A time series of class \code{xts}.
#' @param filter_parameter Filter parameter for flow separation algorithm. 
#'   Default value is that recommended by Nathan and McMahon (1990).
#' @param passes The number of times to pass the filter over the data.
#' @return An \code{xts} object containing the columns "baseflow" and 
#'   "quickflow".
#'
#' @references Nathan, R. J., and T. A. McMahon. "Evaluation of automated 
#'   techniques for base flow and recession analyses." Water Resources 
#'   Research 26.7 (1990): 1465-1473. 
#'   \url{http://dx.doi.org/10.1029/WR026i007p01465}.
#'
#' @export
separate_flow = function(ts, filter_parameter = 0.925, passes = 3){
  f = BaseflowSeparation(coredata(ts), filter_parameter = filter_parameter, 
    passes = passes)
  names(f) = c("baseflow", "quickflow")
  return(xts(f, order.by = index(ts)))
}

#' Time To Rise
#'
#' Compute the time for flow to rise up to or above a given threshold.
#'
#' @param ts A time series of class \code{xts}.
#' @param ut The upper flow threshold for which the time to rise is 
#'   computed. If omitted, the maximum flow of the series is used.
#' @param lt The lower flow threshold from which to start computing 
#'   the time to rise. If omitted, the minimum flow of the series is 
#'   used.
#' @return The time to rise, in a format dependent on the value of 
#'   \code{index(ts)}.
#'
#' @examples
#' data(siouxcity)
#' lower = quantile(siouxcity, 0.25)
#' upper= quantile(siouxcity, 0.75)
#' time_to_rise(siouxcity, upper, lower)
#'
#' @export
time_to_rise = function(ts, ut, lt){
  if(missing(lt))
    lt = min(ts)
  if(missing(ut))
    ut = max(ts)
  u = which(ts >= ut)
  if(length(u) < 1)
    return(NA)
  else
    high = min(u)
  l = which(ts[1:high] <= lt)
  if(length(l) < 1)
    return(NA)
  else
    low = max(l)
  return(index(ts)[high] - index(ts)[low])
}

#' Time To Recede
#'
#' Compute the time for flow to fall down to or below a given threshold.
#'
#' @param ts A time series of class \code{xts}.
#' @param lt The lower flow threshold for which the time to recede is 
#'   computed. If omitted, the minimum flow of the series is used.
#' @param ut The upper flow threshold from which to start computing 
#'   the time to recede. If omitted, the maximum flow of the series is 
#'   used.
#' @return The time to recede, in a format dependent on the value of 
#'   \code{index(ts)}.
#'
#' @examples
#' data(siouxcity)
#' lower = quantile(siouxcity, 0.25)
#' upper= quantile(siouxcity, 0.75)
#' time_to_recede(siouxcity, lower, upper)
#'
#' @export
time_to_recede = function(ts, lt, ut){
  if(missing(ut))
    ut = max(ts)
  u = which(ts >= ut)
  if(length(u) < 1)
    return(NA)
  else
    high = max(u)
  if(missing(lt))
    lt = min(ts[high:length(index(ts))])
  l = which(ts[high:length(index(ts))] <= lt)
  if(length(l) < 1)
    return(NA)
  else
    low = min(l)
  return(index(ts)[low] - index(ts)[high])
}

#' High Flow Duration
#'
#' Compute the longest continuous period during which flow is 
#'   at or above a given threshold.
#'
#' @param ts A time series of class \code{xts}.
#' @param ut The upper flow threshold used to compute duration.
#' @return The duration of the longest period where flow is at
#'   or above the threshold, in a format dependent on the value 
#'   of \code{index(ts)}.
#'
#' @examples
#' data(siouxcity)
#' upper = quantile(siouxcity, 0.75)
#' high_flow_duration(siouxcity, upper)
#'
#' @export
high_flow_duration = function(ts, ut){
  pd = ts
  pd[ts < ut] = NA
  m = na.contiguous(pd)
  if(length(m) < 1)
    return(NA)
  else
    return(tail(index(m), 1) - head(index(m), 1))
}

#' Low Flow Duration
#'
#' Compute the longest continuous period during which flow is 
#'   at or below a given threshold.
#'
#' @param ts A time series of class \code{xts}.
#' @param lt The lower flow threshold used to compute duration.
#' @return The duration of the longest period where flow is at
#'   or below the threshold, in a format dependent on the value 
#'   of \code{index(ts)}.
#'
#' @examples
#' data(siouxcity)
#' lower = quantile(siouxcity, 0.25)
#' low_flow_duration(siouxcity, lower)
#'
#' @export
low_flow_duration = function(ts, lt){
  pd = ts
  pd[ts > lt] = NA
  m = na.contiguous(pd)
  if(length(m) < 1)
    return(NA)
  else
    return(tail(index(m), 1) - head(index(m), 1))
}

#' Total Time Above Threshold
#'
#' Compute the total amount of time that flow is 
#'   at or above a given threshold.
#'
#' @param ts A time series of class \code{xts}.
#' @param ut The upper flow threshold.
#' @return The total amount of time that flow is at or 
#'   above the threshold, in a format dependent on the 
#'   value of \code{index(ts)}.
#'
#' @examples
#' data(siouxcity)
#' upper = quantile(siouxcity, 0.75)
#' total_time_above_threshold(siouxcity, upper)
#'
#' @export
total_time_above_threshold = function(ts, ut){
  pd = ts
  d = which(pd >= ut)
  if(length(d) < 1)
    return(NA)
  dt = diff(index(pd))
  if(length(unique(dt)) > 1)
    warning("Time series is irregular. Total time calculation ",
      "may be erroneous.")
  dt = as.difftime(c(dt, dt[length(dt)]), units = units(dt))
  sum(dt[d])
}

#' Total Time Below Threshold
#'
#' Compute the total amount of time that flow is 
#'   at or below a given threshold.
#'
#' @param ts A time series of class \code{xts}.
#' @param lt The lower flow threshold.
#' @return The total amount of time that flow is at or 
#'   below the threshold, in a format dependent on the 
#'   value of \code{index(ts)}.
#'
#' @examples
#' data(siouxcity)
#' lower = quantile(siouxcity, 0.25)
#' total_time_below_threshold(siouxcity, lower)
#'
#' @export
total_time_below_threshold = function(ts, lt){
  pd = ts
  d = which(pd <= lt)
  if(length(d) < 1)
    return(NA)
  dt = diff(index(pd))
  if(length(unique(dt)) > 1)
    warning("Time series is irregular. Total time calculation ",
      "may be erroneous.")
  dt = as.difftime(c(dt, dt[length(dt)]), units = units(dt))
  sum(dt[d])
}

#' Number of High Flow Pulses
#'
#' Compute the number of high flow pulses (peaks) above a threshold.
#' 
#' @param ts A time series of class \code{xts}.
#' @param width The window size within which to detect peaks. Must be odd.
#' @param ut The upper flow threshold above which to identify peaks.
#' @param ... Other arguments passed to method \code{rollapply} (except for 
#'   \code{align}; see 'details' for more information).
#' @return The number of peaks above the threshold.
#'
#' @details The method \code{rollapply} identifies peaks by testing whether the 
#'   value at the center of the rolling window is the maximum value within the 
#'   window. The function therefore sets the \code{rollapply} argument 
#'   \code{align = "center"}. Larger windows can be used to essentially reduce 
#'   the tolerance for what is considered a 'peak' in the presence of noise.
#'
#' @examples
#' data(siouxcity)
#' number_of_pulses(siouxcity["2005-01-01::2005-01-23"])
#' number_of_pulses(siouxcity["2005-01-01::2005-01-23"], width = 5)
#' number_of_pulses(siouxcity["2005-01-01::2005-01-23"], width = 5, ut = 15000)
#'
#' @export
number_of_pulses = function(ts, width = 3, ut = 0, ...){
  if(as.integer(width) != width){
    warning("Rounding argument 'width' to nearest integer.")
    width = round(width)
  }
  if(width < 1 | width %% 2 == 0)
    stop("Argument 'width' must be an odd positive integer.")    
  mid = (width - 1) %/% 2 + 1
  peaks = coredata(rollapply(ts, width, function(x) which.max(x) == mid, 
    align = "center", ...))
  peaks[ts < ut] = FALSE
  sum(peaks, na.rm = TRUE)
}


