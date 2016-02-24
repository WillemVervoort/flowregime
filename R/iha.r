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
#' @return A two-column dataframe containing the IHA parameter names 
#'   (\code{parameters}) and computed values (\code{values}).
#'
#' @details The IHA method requires a regular time series with no missing 
#'   values. 
#' @references Richter, B. D., Baumgartner, J. V., Powell, J. and Braun, D. P.
#'   (1996), A Method for Assessing Hydrologic Alteration within Ecosystems.
#'   Conservation Biology, 10: 1163–1174. doi: 10.1046/j.1523-1739.1996.10041163.x
#'
#' @examples
#' data(siouxcity)
#' iha(siouxcity['2009/2011'], ut = 32000, lt = 12000)
#' iha(siouxcity['2009-10-01/2011-09-30'], yearstart = "10-01", 
#'   yearend = "09-30", ut = 32000, lt = 12000)
#'
#' @export
iha = function(ts, yearstart = "01-01", yearend = "12-31", groups = 1:5, 
  ut, lt){
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
    records[[n]]["period"] = n
  }
  res = do.call(rbind.data.frame, records)
  rownames(res) = NULL
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