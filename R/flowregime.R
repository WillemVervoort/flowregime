#' Environmental Flow Regime Components
#'
#' A package for computing environmental flow regime components
#' from hydrologic time series. See the vignette to get started.
#' @name flowregime-package
#' @docType package
#' @import xts
NULL

#' Time To Rise
#'
#' Compute the time for flow to rise up to or above a given threshold.
#'
#' @param ts A time series, e.g. of class \code{zoo}.
#' @param ut The upper flow threshold for which the time to rise is 
#'   computed. If omitted, the maximum flow of the series is used.
#' @param lt The lower flow threshold from which to start computing 
#'   the time to rise. If omitted, the minimum flow of the series is 
#'   used.
#' @return the time to rise, in a format dependent on the value of 
#'   \code{index(ts)}.
#'
#' @examples
#' testts = xts(c(sort(rlnorm(1000)), rev(sort(rlnorm(489)))), 
#'   order.by = seq(as.POSIXct("2010-01-01 00:00:00"), 
#'   as.POSIXct("2010-02-01 00:00:00"), length.out = 1489), 
#'   frequency=2*24*365)
#' lower = quantile(testts, 0.25)
#' upper= quantile(testts, 0.75)
#' time_to_rise(testts, upper, lower)
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
#' @param ts A time series, e.g. of class \code{zoo}.
#' @param lt The lower flow threshold for which the time to recede is 
#'   computed. If omitted, the minimum flow of the series is used.
#' @param ut The upper flow threshold from which to start computing 
#'   the time to recede. If omitted, the maximum flow of the series is 
#'   used.
#' @return the time to recede, in a format dependent on the value of 
#'   \code{index(ts)}.
#'
#' @examples
#' testts = zoo(c(sort(rlnorm(1000)), rev(sort(rlnorm(489)))), 
#'   order.by = seq(as.POSIXct("2010-01-01 00:00:00"), 
#'   as.POSIXct("2010-02-01 00:00:00"), length.out = 1489), 
#'   frequency=2*24*365)
#' lower = quantile(testts, 0.25)
#' upper= quantile(testts, 0.75)
#' time_to_recede(testts, lower, upper)
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
#' @param ts A time series, e.g. of class \code{zoo}.
#' @param ut The upper flow threshold used to compute duration.
#' @return The duration of the longest period where flow is at
#'   or above the threshold, in a format dependent on the value 
#'   of \code{index(ts)}.
#'
#' @examples
#' testts = zoo(c(sort(rlnorm(1000)), rev(sort(rlnorm(489)))), 
#'   order.by = seq(as.POSIXct("2010-01-01 00:00:00"), 
#'   as.POSIXct("2010-02-01 00:00:00"), length.out = 1489), 
#'   frequency=2*24*365)
#' upper = quantile(testts, 0.75)
#' high_flow_duration(testts, upper)
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
#' @param ts A time series, e.g. of class \code{zoo}.
#' @param lt The lower flow threshold used to compute duration.
#' @return The duration of the longest period where flow is at
#'   or below the threshold, in a format dependent on the value 
#'   of \code{index(ts)}.
#'
#' @examples
#' testts = zoo(c(sort(rlnorm(1000)), rev(sort(rlnorm(489)))), 
#'   order.by = seq(as.POSIXct("2010-01-01 00:00:00"), 
#'   as.POSIXct("2010-02-01 00:00:00"), length.out = 1489), 
#'   frequency=2*24*365)
#' lower = quantile(testts, 0.25)
#' low_flow_duration(testts, lower)
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
#' @param ts A time series, e.g. of class \code{zoo}.
#' @param ut The upper flow threshold.
#' @return The total amount of time that flow is at or 
#'   above the threshold, in a format dependent on the 
#'   value of \code{index(ts)}.
#'
#' @examples
#' testts = zooreg(c(sort(rlnorm(1000)), rev(sort(rlnorm(489)))), 
#'   order.by = seq(as.POSIXct("2010-01-01 00:00:00"), 
#'   as.POSIXct("2010-02-01 00:00:00"), length.out = 1489), 
#'   frequency=2*24*365)
#' upper = quantile(testts, 0.75)
#' total_time_above_threshold(testts, upper)
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
#' @param ts A time series, e.g. of class \code{zoo}.
#' @param lt The lower flow threshold.
#' @return The total amount of time that flow is at or 
#'   below the threshold, in a format dependent on the 
#'   value of \code{index(ts)}.
#'
#' @examples
#' testts = zooreg(c(sort(rlnorm(1000)), rev(sort(rlnorm(489)))), 
#'   order.by = seq(as.POSIXct("2010-01-01 00:00:00"), 
#'   as.POSIXct("2010-02-01 00:00:00"), length.out = 1489), 
#'   frequency=2*24*365)
#' lower = quantile(testts, 0.25)
#' total_time_below_threshold(testts, lower)
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


