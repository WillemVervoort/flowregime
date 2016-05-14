#' IHA Range of Variability Approach
#' 
#' Analyse IHA results using the Range of Variability Approach.
#'
#' @param pre The pre-impact IHA statistics, i.e. the output of 
#'   \code{IHA(...)}.
#' @param post The post-impact IHA statistics.
#'@param rvacat The RVA categories, as defined by e.g. RVA_categories(pre)
#' @return A 5-column dataframe containing the Hydrologic Alteration Factor for 
#'   each parameter and RVA category, and the upper and lower bounds defining 
#'   the RVA categories for each parameter.
#'
#' @details The \code{boundaries} argument defines cutoffs relative to the mean
#'   (if \code{parametric = TRUE}) or median (if \code{parametric = FALSE}), 
#'   e.g. \code{boundaries = c(-17, 17)} represents cutoffs at the 33rd and 
#'   67th percentile when \code{parametric = FALSE}, and 
#'   \code{boundaries = c(-0.5, 0.5)} represents cutoffs at +/- 0.5 standard
#'   deviations from the mean when \code{parametric = TRUE}.
#'
#' @importFrom stats setNames
#' @export
RVA = function(pre, post, rvacat){
  if(!identical(sort(pre$parameter), sort(post$parameter)))
    stop("'pre' and 'post' analyses do not match.")
  check_RVA_categories(rvacat, pre)
  # define RVA boundaries
  # get pre and post counts for each RVA category
  pre.freq = rva_freq(pre, rvacat)
  post.freq = rva_freq(post, rvacat)
  # calculate HA Factor
  ratio = length(unique(post$YoR))/length(unique(pre$YoR))
  get_HAF = function(o, e) 
    (o/ratio - e)/e
  params = names(pre.freq)
  HAF = setNames(vector("list", length(params)), params)
  for(i in params){
    HAF[[i]] = data.frame(parameter = i, 
      boundary.lower = rvacat[rvacat$parameter == i, 'lower'],
      boundary.upper = rvacat[rvacat$parameter == i, 'upper'],
      HAF.lower = get_HAF(post.freq[[i]]$low, pre.freq[[i]]$low),
      HAF.middle = get_HAF(post.freq[[i]]$mid, pre.freq[[i]]$mid),
      HAF.upper = get_HAF(post.freq[[i]]$high, pre.freq[[i]]$high),
      row.names = NULL, stringsAsFactors = FALSE)
  }
  res = do.call(rbind.data.frame, HAF)
  rownames(res) = NULL
  res
}

# get category counts
rva_freq = function(iha, cats){
  params = unique(iha$parameter)
  res = setNames(vector("list", length(params)), params)
  for(i in params){
    vals = iha[iha$parameter == i, 'value']
    low = vals < cats[cats$parameter == i, 'lower']
    high = vals > cats[cats$parameter == i, 'upper']
    mid = !(low | high)
    res[[i]] = data.frame(low = sum(low), mid = sum(mid), high = sum(high))
  }
  res
}

#' Define RVA Category Boundaries
#'
#' Compute the cutoff values for specifying the three categories used in the
#' RVA approach.
#'
#' @param pre The pre-impact IHA statistics, i.e. the output of 
#'   \code{IHA(...)}.
#' @param boundaries a two-element numeric vector of cutoffs that define the
#'   three RVA categories.
#' @param parametric A logical vector with 1 or 2 elements: Compute lower and
#'   upper boundaries for RVA categories based on parametric (mean) or 
#'   non-parametric (median) analysis.
#' @param na.rm Logical: Remove NA values when computing means/medians.
#'
#' @return A 4-column dataframe containing the upper and lower bounds of the 
#'   RVA categories for each parameter.
#'
#' @importFrom stats quantile
#' @importFrom stats setNames
#' @export
build_RVA_categories = function(pre, boundaries, parametric = c(FALSE, FALSE), 
  na.rm = FALSE){
  if(!all(is.logical(parametric)))
    stop("Value of argument 'parametric' not recognized.")
  if(length(parametric) < 2)
    parametric = rep(parametric, 2)
  if(missing(boundaries)){
    boundaries = c(
      if(parametric[[1]]) -1 else -17,
      if(parametric[[2]]) 1 else 17
    )
    warning("Argument 'boundaries' not specified. Using default values ", 
      paste(boundaries, collapse = " and "), ".")
  }
  bnd = boundaries[order(boundaries)]
  para = parametric[order(boundaries)]
  pfun = function(v, i) 
    mean(v, na.rm = na.rm) + i*sd(v)
  npfun = function(v, i) 
    quantile(v, 0.01*(50 + i), na.rm = na.rm)[[1]]
  lfun = if(para[[1]]) pfun else npfun
  ufun = if(para[[2]]) pfun else npfun  
  params = unique(pre$parameter)
  parambounds = setNames(vector("list", length(params)), params)
  for(i in params){
    vals = pre[pre$parameter == i, "value"]
    parambounds[[i]] = data.frame(parameter = i, lower = lfun(vals, bnd[[1]]), 
      upper = ufun(vals, bnd[[2]]), row.names = NULL, stringsAsFactors = FALSE)
  }
  res = do.call(rbind.data.frame, parambounds)
  rownames(res) = NULL
  res
}

#' Check RVA Categories
#'
#' Check the suitability of RVA categories. If category cutoffs are defined 
#' outside the range of values in the pre-impact dataset for any parameter, 
#' a warning will be issued.
#'
#' @param rvacat The RVA categories, as defined by e.g. 
#'   \code{build_RVA_categories(...)}.
#' @param pre The pre-impact IHA statistics used to define the RVA categories, 
#'   i.e. the output of \code{IHA(...)}.
#' @return A dataframe of 3 columns:
#'
#' @export
check_RVA_categories = function(rvacat, pre){
  params = unique(pre$parameter)
  if(!all(sort(params) == sort(rvacat$parameter)))
    stop("Parameters in argument 'pre' do not match those in argument 'rvacat'.")
  lower.ok = setNames(vector("logical", length(params)), params)
  upper.ok = lower.ok
  for(i in params){
    vals = pre[pre$parameter == i, 'value']
    lower.ok[[i]] = if(rvacat[rvacat$parameter == i, 'lower'] > min(vals)) TRUE else FALSE
    upper.ok[[i]] = if(rvacat[rvacat$parameter == i, 'upper'] < max(vals)) TRUE else FALSE
  }
  if(any(!upper.ok | !lower.ok))
    warning("The following RVA categories definitions should be checked: '",
      paste(params[!upper.ok | !lower.ok], collapse = "', '"), "'.")
  data.frame(parameter = params, lower.ok = lower.ok, upper.ok = upper.ok,
    row.names = NULL, stringsAsFactors = FALSE)
}


#' Replace RVA Categories
#'
#' Replace faulty RVA category definitions.
#'
#' @param rvacat The RVA category specification containing faulty values, i.e.
#'   the output of \code{build_RVA_categories(...)}.
#' @param faulty Identifies elements of \code{rvacat} that are faulty, e.g. the
#'   output of \code{check_RVA_categories(rvacat, ...)}.
#' @param newcat The new RVA category specifications, i.e. the output of 
#'   \code{build_RVA_categories(...)}.
#' @return A combination of \code{rvacat} and \code{newcat}, where the faulty
#'   elements of \code{rvacat} as defined by \code{mask} are replaced with
#'   values from \code{newcat}
#'
#' @export
replace_RVA_categories = function(rvacat, faulty, newcat){
  mask = !faulty[order(faulty, rvacat), 2:3]
  newcat = newcat[order(newcat, rvacat),]
  rvacat[, 2:3][mask] = newcat[, 2:3][mask]
  rvacat
}
