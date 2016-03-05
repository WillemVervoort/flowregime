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
  ut, lt, keep.raw = TRUE){
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
  # define groups to compute
  groupfuns = list(group1, group2, group3, function(x) group4(x, ut, lt), 
    group5)[groups]
  # calculate groups for each year in record
  records = setNames(vector('list', length(sets)), sets)
  for(n in names(records)){
    r = pd[n]
    res = lapply(groupfuns, function(fun) fun(r))
    records[[n]] = do.call(rbind.data.frame, res)
    records[[n]]["YoR"] = n
  }
  res = do.call(rbind.data.frame, records)
  rownames(res) = NULL
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

group1 = function(r){
  mnths = unique(format(index(r), "%b"))
  res = sapply(mnths, function(x) 
    mean(coredata(r)[which(format(index(r), "%b") == x)]))
  data.frame(parameter = paste0("mean (", mnths, ")"), value = res, 
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
  data.frame(parameter = names(res), value = res, row.names = NULL)
}

group3 = function(r){
  max1 = format(index(r)[which.max(coredata(r))], "%m-%d")
  min1 = format(index(r)[which.min(coredata(r))], "%m-%d")
  jds = setNames(seq(1:366), format(seq(as.Date("2004-01-01"), 
    as.Date("2004-12-31"), by = "1 day"), "%m-%d"))
  max1 = jds[[which(max1 == names(jds))]]
  min1 = jds[[which(min1 == names(jds))]]
  data.frame(parameter = c("minima JD", "maxima JD"), value = c(min1, max1))
}

group4 = function(r, ut, lt){
  hp = number_of_high_pulses(r, ut, which = FALSE)
  lp = number_of_low_pulses(r, lt, which = FALSE)
  hpd = mean_high_pulse_duration(r, ut)
  lpd = mean_low_pulse_duration(r, lt)
  res = c("no. low pulses" = lp, "mean low pulse duration" = lpd, 
    "no. high pulses" = hp, "mean high pulse duration" = hpd)
  data.frame(parameter = names(res), value = res, row.names = NULL)
  
}

group5 = function(r){
  diffs = tail(diff(r), -1)
  whichpos = which(diffs > 0)
  whichneg = which(diffs < 0)
  meanpos = mean(diffs[whichpos])
  meanneg = mean(diffs[whichneg])
  reversalfun = function(x)
    ifelse(((x[2] > x[1]) && (x[3] < x[2])) || (x[2] < x[1]) && (x[3] > x[2]), 
      TRUE, FALSE)
  reversals = sum(rollapply(r, 3, reversalfun, align = "right", fill = NULL))
  res = c("mean rise rate" = meanpos, "mean fall rate" = meanneg, 
    "no. reversals" = reversals)
  data.frame(parameter = names(res), value = res, row.names = NULL)
}

#' Compare IHA Results
#'
#' Compare the results of two IHA analyses, e.g. pre- and post-impact periods.
#'
#' @param pre The pre-impact IHA statistics, i.e. the output of 
#'   \code{IHA(..., stat = TRUE)}.
#' @param post The post-impact IHA statistics.
#' @param as.percent If \code{TRUE}, return the differences as a relative 
#'   percent difference (\code{(post - pre)/pre}). Otherwise, return the 
#'   magnitude of difference.
#' @param cl Logical: If \code{TRUE}, compute confidence limits for each 
#'   parameter.
#' @return A dataframe containing either the magnitude of difference or 
#'   relative percent difference of each parameter, and (optionally) the upper
#'   and lowe confidence intervals.
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
  if(!identical(sort(pre$parameters), sort(post$parameters)))
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


date_stat = function(v, statfun, yearstart, yearend){
  starti = as.integer(substring(yearstart, 1, 2))
  endi = as.integer(substring(yearend, 1, 2)) 
  vmon = as.integer(substring(v, 1,2))
  if(starti > endi){
    # water year
    vdates = c(
      as.Date(paste0("0003-", v[vmon >= starti])),
      as.Date(paste0("0004-", v[vmon <= endi]))
    )
  } else {
    vdates = as.Date(paste("0000-", v))
  }
  statfun(vdates)
}


circ_stat = function(v, statfun){
  q1bin = v < 92
  q2bin = (v > 91) && (v < 184)
  q3bin = (v > 183) && (v < 276)
  q4bin = v > 275
  whichbin = which.max(c(q1bin, q2bin, q3bin, q4bin))
  if(whichbin == 1 & (q3bin >= 0.1*length(v) || q4bin >= 0.1*length(v))){
    warning("spread")
  }
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
  res
}
