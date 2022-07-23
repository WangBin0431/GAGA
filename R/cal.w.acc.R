#' Title
#'
#' @param predictions
#' @param truelabels
#'
#' @return
#' @export cal.w.acc
#'
#' @examples
cal.w.acc <- function(predictions, truelabels) {

  tmp = table(truelabels)
  # weights = tmp/sum(tmp)
  # print(weights)
  lvls <- levels(factor(truelabels))
  accs <- lapply(lvls, function(x) {
    idx <- which(truelabels == x)

    return(cal.acc(predictions[idx], truelabels[idx]))

  })
  #browser()
  acc <- mean(unlist(accs))

  return(acc)

}
