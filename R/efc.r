#' Default thresholds for EFC Analysis
#'
#' @param ts A time series of class \code{xts}.
#' @param method The method used to compute environmental flow components:
#'   may be either 'standard' or 'advanced'.
#' @return a list containing the following thresholds:
#'
#'
#' @seealso \link{EFC}
#'
#' 
EFC_default_thresholds = function(ts, method = c("standard", "advanced")){
  method = match.arg(tolower(method), c("standard", "advanced"))
}

#' IHA Environmental Flow Components
#'
#' Compute the IHA Environmental Flow Components (EFCs).
#'
#' @param ts A time series of class \code{xts}.
#' @param yearstart A character vector specifying the month and day signifying
#'   the start of the water year, in the format "mm-dd". Default is "01-01".
#' @param yearend A character vector specifying the month and day signifying
#'   the end of the water year, in the format "mm-dd". Default is "12-31".
#' @param method The method used to compute environmental flow components:
#'   may be either 'standard' or 'advanced'.
#' @param thresholds A list containing all or some of the following elements:
#' @return An \code{xts} time series of EFC factors.
#'
#'
#'
EFC = function(ts, yearstart = "01-01", yearend = "12-31", 
  method = c("standard", "advanced"), thresholds = list()){
  method = match.arg(tolower(method), c("standard", "advanced"))
  if(length(thresholds) < 1)
    thresholds = EFC_default_thresholds(ts, method)
  if(method == "standard")
    if(!all(c("high flow") %in% names(thresholds)))
      stop("Missing elements in argument 'thresholds'.")
  else
    if(!all(c("high flow", "low flow", "high flow start rate", 
      "high flow end rate") %in% names(thresholds)))
      stop("Missing elements in argument 'thresholds'.")
  highclasses = names(thresholds)[c("high flow",  "small flood minimum", 
    "large flood minimum") %in% names(thresholds)]
  lowclasses = names(thresholds)[c("low flow",  "extreme low flow") %in% 
    names(thresholds)]
  message("High flow classes specified: ", paste(highclasses, sep = ", "))
  message("Low flow classes specified: ", paste(lowclasses, sep = ", "))
  
  
  if(method == "standard")
    pass1 = ifelse(coredata(ts) > thresholds[["high flow"]], "high flow pulse", 
      "low flow")
  else{
    pass1 = rep(NA, length(ts))
    pass[1] = ifelse(coredata(ts)[1] > thresholds[["low flow"]], 
      ifelse(coredata(ts)[1] > thresholds[["high flow"]], 
        "ascending limb", "descending limb"), 
      "low flow")
  }
  
}
