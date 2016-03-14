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

#' Roanoke River Flows at Roanoke Rapids, NC
#' 
#' Roanoke River daily discharge recorded at USGS gauge 02080500 in 
#'   Roanoke Rapids, NC from 1913-01-01 to 1991-12-31.
#' @docType data
#' @keywords datasets
#' @name roanokerapids
#' @usage data(roanokerapids)
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

#' Longest High Flow Duration
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
#' longest_high_flow_duration(siouxcity['2011'], 70000)
#' longest_high_flow_duration(siouxcity['2011'], 70000, which = TRUE)
#'
#' @export
longest_high_flow_duration = function(ts, ut, which = FALSE){
  if(any(is.na(ts)))
    warning("NA values detected. Result may be erroneous.")
  pd = ts
  pd[ts < ut] = NA
  if(all(is.na(pd)))
    if(which)
      return(integer(0))
    else
      return(0)
  m = na.contiguous(pd)
  if(which)
    index(ts)[seq(from = which(index(ts) == head(index(m), 1)), 
      to = which(index(ts) == tail(index(m), 1)), by = 1)]
  else
    tail(index(m), 1) - head(index(m), 1)
}

#' Longest Low Flow Duration
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
#' longest_low_flow_duration(siouxcity['2006-06/2007-06'], 18000)
#' longest_low_flow_duration(siouxcity['2006-06/2007-06'], 18000, which = TRUE)
#'
#' @export
longest_low_flow_duration = function(ts, lt, which = FALSE){
  if(any(is.na(ts)))
    warning("NA values detected. Result may be erroneous.")
  pd = ts
  pd[ts > lt] = NA
  if(all(is.na(pd)))
    if(which)
      return(integer(0))
    else
      return(0)
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
#' @param min.dur The minimum duration required for a pulse to be counted.
#' @param which Logical: If \code{TRUE}, return the index locations of the 
#'   peaks instead of the count.
#' @return The number of peaks above the threshold.
#'
#' @examples
#' data(siouxcity)
#' number_of_high_pulses(siouxcity['2009'], ut = 32000)
#' number_of_high_pulses(siouxcity['2009'], ut = 32000, which = TRUE)
#'
#' @export
number_of_high_pulses = function(ts, ut, min.dur = NA, which = FALSE){
  if(any(is.na(ts)))
    warning("NA values detected. Result may be erroneous.")
  pd = ts
  pd[ts < ut] = NA
  if(all(is.na(pd)))
    return(0)
  idx <- 1 + cumsum(is.na(coredata(pd)))
  not.na <- !is.na(coredata(pd))
  pulses = split(index(pd)[not.na], idx[not.na])
  if(!is.na(min.dur))
    pulses = pulses[sapply(pulses, function(x) max(x) - min(x) >= min.dur)]
  pulselocs = sapply(pulses, which.max)
  if(which)
    index(ts)[pulselocs]
  else
    length(pulselocs)
}

#' Number of Low Flow Pulses
#'
#' Compute the number of low flow pulses below a threshold.
#' 
#' @param ts A time series of class \code{xts}.
#' @param lt The lower flow threshold below which to identify pulses.
#' @param min.dur The minimum duration required for a pulse to be counted.
#' @param which Logical: If \code{TRUE}, return the index locations of the 
#'   pulses instead of the count.
#' @return The number of pulses below the threshold.
#'
#' @examples
#' data(siouxcity)
#' number_of_low_pulses(siouxcity['2009'], lt = 12000)
#' number_of_low_pulses(siouxcity['2009'], lt = 12000, which = TRUE)
#'
#' @export
number_of_low_pulses = function(ts, lt, min.dur = NA, which = FALSE){
  if(any(is.na(ts)))
    warning("NA values detected. Result may be erroneous.")
  pd = ts
  pd[ts >= lt] = NA
  if(all(is.na(pd)))
    return(0)
  idx <- 1 + cumsum(is.na(coredata(pd)))
  not.na <- !is.na(coredata(pd))
  pulses = split(index(pd)[not.na], idx[not.na])
  if(!is.na(min.dur))
    pulses = pulses[sapply(pulses, function(x) max(x) - min(x) >= min.dur)]
  if(which){
    pulselocs = sapply(pulses, which.min)
    index(ts)[pulselocs]
  } else
    length(pulses)
}

#' Mean Duration of High Flow Pulses
#'
#' Compute the average duration of high flow pulses above a threshold.
#'
#' @param ts A time series of class \code{xts}.
#' @param ut The upper flow threshold above which to identify pulses.
#' @return The mean duration of high flow pulses.
#'
#' @examples
#' data(siouxcity)
#' mean_high_pulse_duration(siouxcity['2009'], ut = 32000)
#'
#' @export
mean_high_pulse_duration = function(ts, ut){
  if(any(is.na(ts)))
    warning("NA values detected. Result may be erroneous.")
  interval = diff(index(ts))
  if(length(unique(interval)) > 1){
    warning("time series is irregular. Results may be erroneous.")
  }
  interval = interval[[1]]
  pd = ts
  pd[ts < ut] = NA
  if(all(is.na(pd)))
    return(0)
  idx <- 1 + cumsum(is.na(coredata(pd)))
  not.na <- !is.na(coredata(pd))
  pulses = split(index(pd)[not.na], idx[not.na])
  pulselengths = sapply(pulses, function(x) interval + max(x) - min(x))
  mean(pulselengths)
}

#' Mean Duration of Low Flow Pulses
#'
#' Compute the average duration of low flow pulses below a threshold.
#'
#' @param ts A time series of class \code{xts}.
#' @param lt The lower flow threshold below which to identify pulses.
#' @return The mean duration of high flow pulses.
#'
#' @examples
#' data(siouxcity)
#' mean_low_pulse_duration(siouxcity['2009'], lt = 12000)
#'
#' @export
mean_low_pulse_duration = function(ts, lt){
  if(any(is.na(ts)))
    warning("NA values detected. Result may be erroneous.")
  interval = diff(index(ts))
  if(length(unique(interval)) > 1)
    warning("time series is irregular. Results may be erroneous.")
  interval = interval[[1]]
  pd = ts
  pd[ts > lt] = NA
  if(all(is.na(pd)))
    return(0)
  idx <- 1 + cumsum(is.na(coredata(pd)))
  not.na <- !is.na(coredata(pd))
  pulses = split(index(pd)[not.na], idx[not.na])
  pulselengths = sapply(pulses, function(x) interval + max(x) - min(x))
  mean(pulselengths)
}

#' Number of No-Flow Days
#'
#' Compute the number of days for which flow is zero (to some tolerance).
#'
#' @param ts A time series of class \code{xts}.
#' @param tol A tolerance factor. Flows less than \code{tol} are considered to 
#'   be essentially zero. 
#' @return The number of no-flow days.
#'
#' @examples
#' data(siouxcity)
#' number_of_no_flow_days(siouxcity['2009'])
#' number_of_no_flow_days(siouxcity['2009'], tol = 0.1)
#'
#' @export
number_of_no_flow_days = function(ts, tol = 0){
  pd = ts
  noflow = pd[coredata(pd) <= tol]
  length(unique(format(index(noflow), "%Y-%m-%d")))
}

#' Average Minimum Flow
#'
#' Compute the n-average minimum flow.
#'
#' @param ts A time series of class \code{xts}. Assumes a regular timeseries.
#' @param n The moving-window size.
#' @param which Logical: If \code{TRUE}, return the index location of the 
#'   n-average minimum flow instead of the flow magnitude.
#' @param return The n-average minimum flow.
#' 
#' @examples
#' @data(siouxcity)
#' average_minimum_flow(siouxcity['2009'], 1)
#' average_minimum_flow(siouxcity['2009'], 7)
#' average_minimum_flow(siouxcity['2009'], 90, which = TRUE)
#'
#' @export
average_minimum_flow = function(ts, n = 1, which = FALSE){
  if(n < 1 || n != round(n))
    stop("argument 'n' must be a positive integer")
  if(length(index(ts)) < n){
    warning("window size 'n' exceeds length of 'ts'. Returning NA")
    return(NA)
  }
  durflow = rollapply(ts, n, mean, align = "center", fill = NULL)
  if(which)
    index(durflow)[which.min(durflow)]
  else
    coredata(durflow)[which.min(durflow)]
}

#' Average Maximum Flow
#'
#' Compute the n-average maximum flow.
#'
#' @param ts A time series of class \code{xts}. Assumes a regular timeseries.
#' @param n The moving-window size.
#' @param which Logical: If \code{TRUE}, return the index location of the 
#'   n-average maximum flow instead of the flow magnitude.
#' @param return The n-average maximum flow.
#' 
#' @examples
#' @data(siouxcity)
#' average_maximum_flow(siouxcity['2009'], 1)
#' average_maximum_flow(siouxcity['2009'], 7)
#' average_maximum_flow(siouxcity['2009'], 90, which = TRUE)
#'
#' @export
average_maximum_flow = function(ts, n = 1, which = FALSE){
  if(n < 1 || n != round(n))
    stop("argument 'n' must be a positive integer")
  if(length(index(ts)) < n){
    warning("window size 'n' exceeds length of 'ts'. Returning NA")
    return(NA)
  }
  durflow = rollapply(ts, n, mean, align = "center", fill = NULL)
  if(which)
    index(durflow)[which.max(durflow)]
  else
    coredata(durflow)[which.max(durflow)]
}

#' Number Of Flow Reversals
#'
#' Compute the number of times that the flow rate of change reverses.
#'
#' @param ts A time series of class \code{xts}. Assumes a regular timeseries.
#' @param which Logical: If \code{TRUE}, return the index location of the 
#'   flow reversals instead of the total number of reversals.
#' @param return The total number of flow reversals.
#' 
#' @examples
#' @data(siouxcity)
#' number_of_reversals(siouxcity['2009-01'])
#' number_of_reversals(siouxcity['2009-01'], which = TRUE)
#'
#' @export
number_of_reversals = function(ts, which = FALSE){
  if(length(index(ts)) < 3){
    warning("Argument 'ts' has less than 3 data points. Returning NA")
    return(NA)
  }
  revfun = function(x) 
    ifelse(all(c(x[[1]], x[[3]]) < x[[2]]) | all(c(x[[1]], x[[3]]) > x[[2]]), 
      TRUE, FALSE)
  revflow = as.logical(coredata(rollapply(ts, 3, revfun, align = 'right', 
    fill = FALSE)))
  if(which)
    index(ts)[revflow]
  else
    sum(revflow)  
}

