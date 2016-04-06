#' IHA Range of Variability Approach
#' 
#' Compute the IHA Range of Variability.
#'
#' @param pre The pre-impact IHA statistics, i.e. the output of 
#'   \code{IHA(..., keep.raw = TRUE)}.
#' @param post The post-impact IHA statistics.
#' @param boundaries a two-element numeric vector of cutoffs that define the
#'   three RVA categories.
#' @param parametric Logical: perform parametric (mean) or non-parametric
#'   (median) analysis.
#' @return A dataframe.
#'
#' @details The \code{boundaries} argument defines cutoffs relative to the mean
#'   (if \code{parametric = TRUE}) or median (if \code{parametric = FALSE}), 
#'   e.g. \code{boundaries = c(-17, 17)} represents cutoffs at the 33rd and 
#'   67th percentile when \code{parametric = FALSE}, and 
#'   \code{boundaries = c(-0.5, 0.5)} represents cutoffs at +/- 0.5 standard
#'   deviations from the mean when \code{parametric = TRUE}.
#'
RVA = function(pre, post, boundaries, parametric = TRUE){
  boundaries = sort(unique(boundaries))
  if(length(boundaries) != 2)
    stop("Argument 'boundaries' must contain two values.")
  if(!identical(sort(pre$parameters), sort(post$parameters)))
    stop("'pre' and 'post' analyses do not match.")
  post = post[match(pre$parameter, post$parameter),]
  pre.raw = attr(pre, "raw")
  post.raw = attr(post, "raw")
  if(is.null(pre.raw) || is.null(post.raw))
    stop("Cannot compute RVA for outputs of IHA(..., keep.raw = FALSE).")
  for(i in unique(pre.raw$parameter)){
    if(parametric)
      bound = mean
  }  

}