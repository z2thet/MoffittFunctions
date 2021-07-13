#' Calculate percent change from baseline
#' testing adding a function to MoffittFunctions package
#' @param baseline Baseline or starting value
#' @param nextone Next value. This could be a vector of values.
#'
#' @export
percentchange <- function(baseline,nextone){
  .check_numeric_input(baseline)
  .check_numeric_input(nextone)
  return((nextone-baseline)/baseline*100) }
