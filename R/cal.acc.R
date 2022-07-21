#' Title
#'
#' @param predictions
#' @param truelabels
#'
#' @return
#' @export
#'
#' @examples
cal.acc <- function(predictions, truelabels) {

  return(length(which(predictions == truelabels)) / length(truelabels))

}