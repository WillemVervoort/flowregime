#' Environmental Flow Regime Components
#'
#' A package for computing environmental flow regime components
#' from hydrologic time series. See the vignette to get started.
#' @name flowregime-package
#' @aliases flowregime
#' @docType package
#' @import zoo
#' @import xts
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
#' @param which Logical: If \code{TRUE}, return the index locations of the 
#'   rising limb instead of the count.
#' @return The time to rise, in a format dependent on the value of 
#'   \code{index(ts)}.
#'
#' @examples
#' data(siouxcity)
#' time_to_rise(siouxcity['2011'])
#' time_to_rise(siouxcity['2011'], which = TRUE)
#' time_to_rise(siouxcity['2011'], lt = 70000)
#'
#' @export
time_to_rise = function(ts, ut, lt, which = FALSE){
  if(missing(lt))
    lt = min(ts)
  if(missing(ut))
    ut = max(ts)
  u = which(ts >= ut)
  if(length(u) < 1)
    if(which)
      return(u)
    else
      return(NA)
  else
    high = min(u)
  l = which(ts[1:high] <= lt)
  if(length(l) < 1)
    if(which)
      return(l)
    else
      return(NA)
  low = max(l)
  if(which)
    index(ts)[seq(from = low, to = high, by = 1)]
  else
    index(ts)[high] - index(ts)[low]
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
#' @param which Logical: If \code{TRUE}, return the index locations of the 
#'   recession limb instead of the count.
#' @return The time to recede, in a format dependent on the value of 
#'   \code{index(ts)}.
#'
#' @examples
#' data(siouxcity)
#' time_to_recede(siouxcity['2011'])
#' time_to_recede(siouxcity['2011'], which = TRUE)
#' time_to_recede(siouxcity['2011'], lt = 70000, which = TRUE)
#'
#' @export
time_to_recede = function(ts, lt, ut, which = FALSE){
  if(missing(ut))
    ut = max(ts)
  u = which(ts >= ut)
  if(length(u) < 1)
    if(which)
      return(u)
    else
      return(NA)
  else
    high = max(u)
  if(missing(lt))
    lt = min(ts[high:length(index(ts))])
  l = which(ts[high:length(index(ts))] <= lt)
  if(length(l) < 1)
    if(which)
      return(l)
    else
      return(NA)
  low = min(l) + high - 1
  if(which)
    index(ts)[seq(from = high, to = low, by = 1)]
  else
    index(ts)[low] - index(ts)[high]
}

#' High Flow Duration
#'
#' Compute the longest continuous period during which flow is 
#'   at or above a given threshold.
#'
#' @param ts A time series of class \code{xts}.
#' @param ut The upper flow threshold used to compute duration.
#' @param which Logical: If \code{TRUE}, return the index locations of the 
#'   duration instead of the count.
#' @return The duration of the longest period where flow is at
#'   or above the threshold, in a format dependent on the value 
#'   of \code{index(ts)}.
#'
#' @examples
#' data(siouxcity)
#' high_flow_duration(siouxcity['2011'], 70000)
#' high_flow_duration(siouxcity['2011'], 70000, which = TRUE)
#'
#' @export
high_flow_duration = function(ts, ut, which = FALSE){
  pd = ts
  pd[ts < ut] = NA
  if(all(is.na(pd)))
    if(which)
      return(integer(0))
    else
      return(NA)
  m = na.contiguous(pd)
  if(which)
    index(ts)[seq(from = which(index(ts) == head(index(m), 1)), 
      to = which(index(ts) == tail(index(m), 1)), by = 1)]
  else
    tail(index(m), 1) - head(index(m), 1)
}

#' Low Flow Duration
#'
#' Compute the longest continuous period during which flow is 
#'   at or below a given threshold.
#'
#' @param ts A time series of class \code{xts}.
#' @param lt The lower flow threshold used to compute duration.
#' @param which Logical: If \code{TRUE}, return the index locations of the 
#'   duration instead of the count.
#' @return The duration of the longest period where flow is at
#'   or below the threshold, in a format dependent on the value 
#'   of \code{index(ts)}.
#'
#' @examples
#' data(siouxcity)
#' low_flow_duration(siouxcity['2006-06/2007-06'], 18000)
#' low_flow_duration(siouxcity['2006-06/2007-06'], 18000, which = TRUE)
#'
#' @export
low_flow_duration = function(ts, lt, which = FALSE){
  pd = ts
  pd[ts > lt] = NA
  if(all(is.na(pd)))
    if(which)
      return(integer(0))
    else
      return(NA)
  m = na.contiguous(pd)
  if(which)
    index(ts)[seq(from = which(index(ts) == head(index(m), 1)), 
      to = which(index(ts) == tail(index(m), 1)), by = 1)]
  else
    tail(index(m), 1) - head(index(m), 1)
}

#' Total Time Above Threshold
#'
#' Compute the total amount of time that flow is 
#'   at or above a given threshold.
#'
#' @param ts A time series of class \code{xts}.
#' @param ut The upper flow threshold.
#' @param which Logical: If \code{TRUE}, return the index locations of the 
#'   high flow records.
#' @return The total amount of time that flow is at or 
#'   above the threshold, in a format dependent on the 
#'   value of \code{index(ts)}.
#'
#' @examples
#' data(siouxcity)
#' total_time_above_threshold(siouxcity['2010'], 60000)
#' total_time_above_threshold(siouxcity['2010'], 60000, which = TRUE)
#'
#' @export
total_time_above_threshold = function(ts, ut, which = FALSE){
  pd = ts
  d = which(pd >= ut)
  if(which)
    return(index(ts)[d])
  if(length(d) < 1)
    return(0)
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
#' @param which Logical: If \code{TRUE}, return the index locations of the 
#'   low flow records.
#' @return The total amount of time that flow is at or 
#'   below the threshold, in a format dependent on the 
#'   value of \code{index(ts)}.
#'
#' @examples
#' data(siouxcity)
#' total_time_below_threshold(siouxcity['2007-06/2008-06'], 15000)
#' total_time_below_threshold(siouxcity['2007-06/2008-06'], 15000, which = TRUE)
#'
#' @export
total_time_below_threshold = function(ts, lt, which = FALSE){
  pd = ts
  d = which(pd <= lt)
  if(which)
    return(index(ts)[d])
  if(length(d) < 1)
    return(0)
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
#' @param ut The upper flow threshold above which to identify peaks.
#' @param ws The window size within which to detect peaks. Must be odd.
#' @param which Logical: If \code{TRUE}, return the index locations of the 
#'   peaks instead of the count.
#' @return The number of peaks above the threshold.
#'
#' @details The method \code{rollapply} identifies peaks by testing whether the 
#'   value at the center of the rolling window is the maximum value within the 
#'   window. Larger windows can be used to essentially reduce the tolerance for 
#'   what is considered a 'peak' in the presence of noise.
#'
#' @examples
#' data(siouxcity)
#' number_of_pulses(siouxcity['2009'])
#' number_of_pulses(siouxcity['2009'], which = TRUE)
#' number_of_pulses(siouxcity['2009'], ws = 7, which = TRUE)
#' number_of_pulses(siouxcity['2009'], ws = 7, ut = 32000)
#'
#' @export
number_of_pulses = function(ts, ut = 0, ws = 3, which = FALSE){
  if(as.integer(ws) != ws){
    warning("Rounding argument 'ws' to nearest integer.")
    ws = round(ws)
  }
  if(ws < 1 | ws %% 2 == 0)
    stop("Argument 'ws' must be an odd positive integer.")    
  mid = (ws - 1) %/% 2 + 1
  peaks = as.logical(coredata(rollapply(ts, ws, function(x) 
    which.max(x) == mid, align = "center", fill = FALSE, by.column = TRUE)))
  peaks[ts < ut] = FALSE
  if(which)
    index(ts)[peaks]
  else
    sum(peaks, na.rm = TRUE)
}


