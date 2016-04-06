#' Default thresholds for EFC Analysis
#'
EFC_default_thresholds = function(ts){

}

#' IHA Environmental Flow Components
#'
#' Compute the IHA Environmental Flow Components (EFCs).
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
    pass[1] = ifelse(cordata(ts)[1] > thresholds[["low flow"]], 
      ifelse(coredata(ts)[1] > thresholds[["high flow"]], 
        "ascending limb", "descending limb"), 
      "low flow")
  }
  
}
