#' Title
#'
#' @param fit
#' @param newx
#'
#' @return
#' @export
#'
#' @examples
predict.GAGA = function(fit,newx){
  if(inherits(fit, "GAGA", which = FALSE)==FALSE)stop("predict.GAGA need a GAGA object as its input")
  family = class(fit)[2]
  Ey=switch(family,
            "gaussian"=predict.LM_GAGA(fit,newx),
            "poisson"=predict.poisson_GAGA(fit,newx),
            "binomial"=predict.logistic_GAGA(fit,newx),
            "multinomial"=predict.multinomial_GAGA(fit,newx),
            "cox"=predict.cox_GAGA(fit,newx)
  )
  return(Ey)
}

predict.LM_GAGA = function(fit,newx){
  Eb = fit$beta
  Ey = newx%*%Eb
  return(Ey)
}

predict.poisson_GAGA = function(fit,newx){
  Eb = fit$beta
  Ey = exp(newx%*%Eb)
  return(Ey)
}

predict.logistic_GAGA = function(fit,newx){
  test_size = nrow(newx)
  classnames = fit$classnames
  Eb = fit$beta
  t = 1/(1+exp(-newx%*%Eb))
  Ey = rep(0,test_size)
  Ey[t>0.5] = 1
  Ey=classnames[Ey+1]
  return(Ey)
}

predict.multinomial_GAGA = function(fit,newx){
  test_size = nrow(newx)
  classnames = fit$classnames
  Eb = fit$beta

  Ey = rep(0,test_size)
  z = newx%*%Eb
  t = exp(z)/(1+rowSums(exp(z)))
  t = cbind(t,1-rowSums(t))
  for(jj in 1:test_size){
    Ey[jj] = which.max(t[jj,])
  }
  Ey = classnames[Ey]
  return(Ey)
}

predict.cox_GAGA = function(fit,newx){
  Eb = fit$beta
  Ey = newx%*%Eb
  return(Ey)
}
