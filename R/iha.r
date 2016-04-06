#' Indicators of Hydrologic Alteration
#'
#' Compute the Indicators of Hydrologic Alteration for a daily time series.
#'
#' @param ts A time series of class \code{xts}.
#' @param yearstart A character vector specifying the month and day signifying
#'   the start of the water year, in the format "mm-dd". Default is "01-01".
#' @param yearend A character vector specifying the month and day signifying
#'   the end of the water year, in the format "mm-dd". Default is "12-31".
#' @param groups The parameter groups to compute. Computes all groups (1-5) by
#'   default.
#' @param ut The upper flow threshold for identifying high flow pulses. IHA 
#'   recommends the 75th flow percentile for pre-impact conditions.
#' @param lt The lower flow threshold for identifying low flow pulses. IHA 
#'   recommends the 25th flow percentile for pre-impact conditions.
#' @param parametric Logical: perform parametric (mean) or non-parametric
#'   (median) analysis.
#' @param keep.raw Logical: return the hydrologic attributes for each year 
#'   of record as attribute 'raw' of the output dataframe. 
#' @return A 3-column dataframe. If \code{stat = TRUE}, the dataframe contains 
#'   the parameter names, the central tendencies and the measures of 
#'   dispersion. If \code{stat = FALSE}, the dataframe contains the parameter 
#'   names, the value of the attribute and the year of record (YoR) for which 
#'   the attribute was computed.
#'
#' @details The IHA method requires a regular time series with no missing 
#'   values. 
#'
#' @section Notes:
#'   \itemize{
#'     \item Low- and high-flow pulses occurring at the start of a year of 
#'       record are ignored when calculating number of pulses and mean duration
#'       of pulses.
#'     \item IHA v7 rounds averages of flow values to the nearest whole number for
#'       individual water years, and this rounding propagates through the 
#'       scorecard. Numerical results may differ slightly from IHAv7.
#'       
#'
#' } 
#'
#' @references Richter, B. D., Baumgartner, J. V., Powell, J. and Braun, D. P.
#'   (1996), A Method for Assessing Hydrologic Alteration within Ecosystems.
#'   Conservation Biology, 10: 1163-1174. doi: 10.1046/j.1523-1739.1996.10041163.x
#'
#' @examples
#' data(siouxcity)
#' IHA(siouxcity['2009/2011'], ut = 32000, lt = 12000)
#' IHA(siouxcity['2009/2011'], ut = 32000, lt = 12000, keep.raw = FALSE)
#' IHA(siouxcity['2009-10-01/2011-09-30'], yearstart = "10-01", 
#'   yearend = "09-30", ut = 32000, lt = 12000)
#'
#' @export
IHA = function(ts, yearstart = "01-01", yearend = "12-31", groups = 1:5, 
  ut, lt, parametric = TRUE, keep.raw = TRUE){
  # argument checking
  if(yearstart == yearend)
    stop("Arguments 'yearstart' and 'yearend' must be different")
  if(!all(groups %in% 1:5))
    stop("Value of argument 'groups' is invalid.")
  if(4 %in% groups){
    if(missing(ut) | missing(lt))
      stop("Arguments 'ut' and 'lt' are required to compute group 4")
  }
  # check that time series is regular
  pd = ts
  # define water years
  startidx = index(pd)[which(format(index(pd), "%m-%d") == yearstart)]
  endidx = index(pd)[which(format(index(pd), "%m-%d") == yearend)]
  if(length(startidx) != length(endidx))
    stop("Flow record is incomplete")
  sets = paste(startidx, endidx, sep = "/")
  # helper function for computing groups
  looper = function(record, periods, fun, ...){
    res = setNames(vector("list", length(periods)), periods)
    for(p in periods){
      r = record[p]
      res[[p]] = fun(r, ...)
      res[[p]]["YoR"] = p
    }
    ret = do.call(rbind.data.frame, res)
    rownames(ret) = NULL
    ret
  }
  # compute groups
  gd = vector("list", 5)
  if(1 %in% groups)
    gd[[1]] = looper(pd, sets, group1, parametric)
  if(2 %in% groups)
    gd[[2]] = looper(pd, sets, group2)
  if(3 %in% groups)
    gd[[3]] = looper(pd, sets, group3)
  if(4 %in% groups)
    gd[[4]] = group4(pd, ut, lt, sets, parametric)
  if(5 %in% groups)
    gd[[5]] = looper(pd, sets, group5, parametric)  
  res = do.call(rbind.data.frame, gd)  
# compute stats for each group
  pnames = unique(res$parameter)
  pcentral = setNames(vector("numeric", length(pnames)), pnames)
  pdispersion = pcentral
  for(p in pnames){
    pcentral[[p]] = mean(res[res$parameter == p, "value"])
    pdispersion[[p]] = sd(res[res$parameter == p, "value"])/pcentral[[p]]
  }
  # fix for circular stats
  for(p in c("minima JD", "maxima JD")){
    pcentral[[p]] = circ_stat(res[res$parameter == p, "value"], mean)
    pdispersion[[p]] = circ_stat(res[res$parameter == p, "value"], sd)/pcentral[[p]]
}
  res2 = data.frame(parameter = pnames, central.tendency = pcentral, 
    dispersion = pdispersion, row.names = NULL)  
  if(keep.raw)
    structure(res2, raw = res)
  else
    res2
}

group1 = function(r, para){
  if(para){
    fun = mean
    tag = "mean"
  } else {
    fun = median
    tag = "median"
  }
  mnths = unique(format(index(r), "%b"))
  res = sapply(mnths, function(x) 
    fun(coredata(r)[which(format(index(r), "%b") == x)]))
  data.frame(parameter = paste0(tag, " (", mnths, ")"), value = res, 
    row.names = NULL)
}

group2 = function(r){
  rollnums = c(1, 3, 7, 30, 90)
  minima = vector("numeric", length = length(rollnums))
  maxima = vector("numeric", length = length(rollnums))
  names(minima) = paste0(rollnums, "-day minima")
  names(maxima) = paste0(rollnums, "-day maxima")
  for(i in seq_along(rollnums)){
    rollr = rollapply(r, rollnums[i], mean, align = "center", fill = NULL)
    minima[[i]] = min(coredata(rollr))
    maxima[[i]] = max(coredata(rollr))
  }
  res = c(minima, maxima)
  res[["no. zero-flow days"]] = number_of_no_flow_days(r, 0)
  res[["baseflow index"]] = baseflow_index(r)

  data.frame(parameter = names(res), value = res, row.names = NULL)
}

group3 = function(r){
  max1 = get_jd(as.POSIXlt(index(r)[which.max(coredata(r))]))
  min1 = get_jd(as.POSIXlt(index(r)[which.min(coredata(r))]))  
  data.frame(parameter = c("minima JD", "maxima JD"), value = c(min1, max1))
}

group4 = function(rec, ut, lt, periods, para){
    if(para){
    fun = mean
    tag = "mean"
  } else {
    fun = median
    tag = "median"
  }
  interval = diff(index(rec))
  if(length(unique(interval)) > 1){
    warning("Time series is irregular. Results may be erroneous.")
  }
  interval = interval[[1]]
  # get pulses and durations
  hpf = number_of_high_pulses(rec, ut, which = TRUE, return.which = "first", 
    ignore.first = TRUE)
  hpl = number_of_high_pulses(rec, ut, which = TRUE, return.which = "last", 
    ignore.first = TRUE)
  hpdur = interval + hpl - hpf
  lpf = number_of_low_pulses(rec, lt, which = TRUE, return.which = "first", 
    ignore.first = TRUE, min.dur = as.difftime(1, units = "days"))
  lpl = number_of_low_pulses(rec, lt, which = TRUE, return.which = "last", 
    ignore.first = TRUE, min.dur = as.difftime(1, units = "days"))
  lpdur = interval + lpl - lpf
  # split by period
  hpidx = lapply(periods, function(x) hpf[which(hpf %in% index(rec[x]))])  
  hduridx = lapply(periods, function(x) hpdur[which(hpf %in% index(rec[x]))])
  lpidx = lapply(periods, function(x) lpf[which(lpf %in% index(rec[x]))])  
  lduridx = lapply(periods, function(x) lpdur[which(lpf %in% index(rec[x]))])
  res = vector("list", length(periods))
  nm = 
  for(i in seq_along(res)){
    resv = setNames(c(length(lpidx[[i]]), fun(lduridx[[i]]), 
      length(hpidx[[i]]), fun(hduridx[[i]])), c("no. low pulses", 
      paste(tag, "low pulse duration"), "no. high pulses", 
      paste(tag, "high pulse duration")))    
    res[[i]] = data.frame(parameter = names(resv), value = resv, 
      YoR = periods[[i]], row.names = NULL)
  }
  do.call(rbind.data.frame, res)
}

group5 = function(r, para){
  if(para){
    fun = mean
    tag = "mean"
  } else {
    fun = median
    tag = "median"
  }
  diffv = tail(diff(coredata(r)), -1)
  difft = tail(diff(index(r)), -1)
  diffs = diffv/as.numeric(difft)
  whichpos = which(diffs > 0)
  whichneg = which(diffs < 0)
  meanpos = fun(diffs[whichpos])
  meanneg = fun(diffs[whichneg])
  reversals = number_of_reversals(r, which = FALSE)
  res = setNames(c(meanpos, meanneg, reversals), c(paste(tag, "rise rate"), 
  paste(tag, "fall rate"), "no. reversals"))
  data.frame(parameter = names(res), value = res, row.names = NULL)
}

#' Compare IHA Results
#'
#' Compare the results of two IHA analyses, e.g. pre- and post-impact periods.
#'
#' @param pre The pre-impact IHA statistics, i.e. the output of 
#'   \code{IHA(..., keep.raw = TRUE)}.
#' @param post The post-impact IHA statistics.
#' @param as.percent If \code{TRUE}, return the differences as a relative 
#'   percent difference (\code{(post - pre)/pre}). Otherwise, return the 
#'   magnitude of difference.
#' @param cl Logical: If \code{TRUE}, compute confidence limits for each 
#'   parameter.
#' @return A dataframe containing either the magnitude of difference or 
#'   relative percent difference of each parameter, and (optionally) the upper
#'   and lower confidence intervals.
#'
#' @examples
#' data(siouxcity)
#' pre = IHA(siouxcity['2004/2009'], ut = 32000, lt = 12000)
#' post = IHA(siouxcity['2010/2014'], ut = 32000, lt = 12000)
#' compareIHA(pre, post)
#' compareIHA(pre, post, as.percent = TRUE, cl = TRUE)
#'
#' @export
compareIHA = function(pre, post, as.percent = FALSE, cl = FALSE){
  if(!identical(sort(pre$parameter), sort(post$parameter)))
    stop("'pre' and 'post' analyses do not match.")
  post = post[match(pre$parameter, post$parameter),]
  ctendiff = post$central.tendency - pre$central.tendency
  dispdiff = post$dispersion - pre$dispersion
  if(as.percent){
    ctendiff = 100*ctendiff/pre$central.tendency
    dispdiff = 100*dispdiff/pre$dispersion
  }
  res = data.frame(pre$parameter, ctendiff, dispdiff)
  if(cl){
    if(is.null(attr(pre, "raw")) || is.null(attr(post, "raw")))
      warning("Cannot compute confidence levels for outputs of ",
        "IHA(..., keep.raw = FALSE).")
    else{
      #res["cl.upper"] = 
      #res["cl.lower"] = 
    }
  }
  res
}

# circular statistics for IHA
circ_stat = function(v, statfun){
  q1bin = v < 92
  q2bin = (v > 91) & (v < 184)
  q3bin = (v > 183) & (v < 276)
  q4bin = v > 275
  # identify largest bin
  q1size = sum(q1bin)
  q2size = sum(q2bin)
  q3size = sum(q3bin)
  q4size = sum(q4bin)
  whichbin = which.max(c(q1size, q2size, q3size, q4size))
  if(whichbin %in% c(2, 3)){
    res = statfun(v)
  } else if(whichbin == 1){
    v[q4bin] = v[q4bin] - 366
    res = statfun(v)
    res = ifelse(res < 0, res + 366, res)
  } else {
    v[q1bin] = v[q1bin] + 366
    res = statfun(v)
    res = ifelse(res > 366, res - 366, res)
  }
  # warn about data scatter
  n = length(v)
  if(((whichbin == 1) && (q3size > 0.1*n)) ||
     ((whichbin == 2) && (q4size > 0.1*n)) ||
     ((whichbin == 3) && (q1size > 0.1*n)) ||
     ((whichbin == 4) && (q2size > 0.1*n)))
    warning("Dates of extreme flows are widely distributed throughout the ",
      "year. Use statistics with caution.")     
  res
}

#' IHA Julian Day
#'
#' Compute Julian Day as per IHA v7.
#'
#' @param d Date or POSIX*t object
#' @return Integer Julian date. 
#'
#' @details IHA defines Jan 1 as JD 1 and Dec 31 as JD 366 for both leap years
#'   and non-leap years. In non-leap years, JD 60 (Feb 29) is skipped, so that
#'   e.g. Feb 28 2003 is JD 59 and Mar 01 2003 is JD 61.
#'
get_jd = function(d){
  year = d$year + 1900
  month = d$mon + 1
  jd = d$yday + 1
  if(!((year %% 4 == 0) && ((year %% 100 != 0) || (year %% 400 == 0))) &&
    month > 2)
    jd + 1
  else
    jd
}  
