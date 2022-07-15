#Ref http://tecdat.cn/

calculate.accuracy <- function(predictions, truelabels) {
  
  return(length(which(predictions == truelabels)) / length(truelabels))
  
}

calculate.w.accuracy <- function(predictions, truelabels) {
  
  tmp = table(truelabels)
  # weights = tmp/sum(tmp)
  # print(weights)
  lvls <- levels(factor(truelabels))
  accs <- lapply(lvls, function(x) {
    idx <- which(truelabels == x)
    
    return(calculate.accuracy(predictions[idx], truelabels[idx]))
    
  })
  #browser()
  acc <- mean(unlist(accs))
  
  return(acc)
  
} 
