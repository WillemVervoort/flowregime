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
#' @param stats Logical: If \code{TRUE}, return the inter-annual statistics 
#'   of each parameter. Otherwise, return the hydrologic attributes for each 
#'   year of record. 
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
#' IHA(siouxcity['2009/2011'], ut = 32000, lt = 12000, stats = FALSE)
#' IHA(siouxcity['2009-10-01/2011-09-30'], yearstart = "10-01", 
#'   yearend = "09-30", ut = 32000, lt = 12000)
#'
#' @export
IHA = function(ts, yearstart = "01-01", yearend = "12-31", groups = 1:5, 
  ut, lt, stats = TRUE){
  # argument checking
  if(yearstart == yearend)
    stop("Arguments 'yearstart' and 'yearend' must be different")
  if(!all(groups %in% 1:5))
    stop("Value of argument 'groups' is invalid.")
  if(4 %in% groups){
    if(missing(ut) | missing(lt))
      stop("Arguments 'ut' and 'lt' are required to compute group 4")
  }
  ## check that time series is regular
  pd = ts
  # define water years
  startidx = index(pd)[which(format(index(pd), "%m-%d") == yearstart)]
  endidx = index(pd)[which(format(index(pd), "%m-%d") == yearend)]
  if(length(startidx) != length(endidx))
    stop("Flow record is incomplete.")
  sets = paste(startidx, endidx, sep = "/")
  # calculate groups for each year in record
  records = setNames(vector('list', length(sets)), sets)
  for(n in names(records)){
    r = pd[n]
    res = list(group1(r), group2(r), group3(r), group4(r, ut, lt), group5(r))
    records[[n]] = do.call(rbind.data.frame, res[groups])
    records[[n]]["YoR"] = n
  }
  res = do.call(rbind.data.frame, records)
  rownames(res) = NULL
  if(!stats)
    return(res)
  # compute stats for each group
  pnames = unique(res$parameter)
  pcentral = setNames(vector("numeric", length(pnames)), pnames)
  pdispersion = pcentral
  for(p in pnames){
    pcentral[[p]] = mean(res[res$parameter == p, "value"])
    pdispersion[[p]] = sd(res[res$parameter == p, "value"])/pcentral[[p]]
  }
  data.frame(parameter = pnames, central.tendency = pcentral, 
    dispersion = pdispersion, row.names = NULL)  
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
#' @return a dataframe containing either the magnitude of difference or 
#'   relative percent difference of each parameter, and (optionally) the upper
#'   and lowe confidence intervals.
#'
#' @examples
#' data(siouxcity)
#' pre = iha(siouxcity['2004/2009'], ut = 32000, lt = 12000)
#' post = iha(siouxcity['2010/2014'], ut = 32000, lt = 12000)
#' compareIHA(pre, post)
#' compareIHA(pre, post, as.percent = TRUE, cl = TRUE)
#'
#' @export
compareIHA = function(pre, post, as.percent = FALSE, cl = FALSE){
  if(!identical(sort(pre$paramters, post$parameters)))
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
    #res["cl.upper"] = 
    #res["cl.lower"] = 
  }
  res
}

group1 = function(r){
  mnths = unique(format(index(r), "%b"))
  res = sapply(mnths, function(x) 
    mean(coredata(r)[which(format(index(r), "%b") == x)]))
  data.frame(parameter = paste0("mean (", mnths, ")"), value = res, 
    row.names = NULL)
}

group2 = function(r, rollnums = c(1, 3, 7, 30, 90)){
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
  max1 = as.POSIXlt(format(index(r)[which.max(coredata(r))], "%Y-%m-%d"))$yday
  min1 = as.POSIXlt(format(index(r)[which.min(coredata(r))], "%Y-%m-%d"))$yday
  data.frame(parameter = c("maxima JD", "minima JD"), value = c(max1,  min1))
}

group4 = function(r, ut, lt){
  hp = number_of_high_pulses(r, ut, which = FALSE)
  lp = number_of_low_pulses(r, lt, which = FALSE)
  hpd = mean_high_pulse_duration(r, ut)
  lpd = mean_low_pulse_duration(r, lt)
  res = c("no. high pulses" = hp, "no. low pulses" = lp, 
    "mean high pulse duration" = hpd, "mean low pulse duration" = lpd)
  data.frame(parameter = names(res), value = res, row.names = NULL)
  
}

group5 = function(r){
  diffs = tail(diff(r), -1)
  whichpos = which(diffs > 0)
  whichneg = which(diffs < 0)
  meanpos = mean(diffs[whichpos])
  meanneg = mean(diffs[whichneg])
  numpos = length(whichpos)
  numneg = length(whichneg)
  res = c("no. of rises" = numpos, "no. of falls" = numneg, 
    "mean rise rate" = meanpos, "mean fall rate" = meanneg)
  data.frame(parameter = names(res), value = res, row.names = NULL)
}