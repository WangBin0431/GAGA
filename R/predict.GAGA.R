#' Get predictions from a GAGA fit object
#'
#' Gives fitted values from a fitted GAGA object.
#' @param fit Fitted "GAGA" object.
#' @param newx Matrix of new values for x at which predictions are to be made. Must be a
#' matrix.  If the intercept term needs to be considered in the estimation process, then
#' the first column of \code{X} must be all 1s.
#'
#' @return Predictions
#' @export predict.GAGA
#'
#' @examples
#' set.seed(2022)
#' p_size = 30
#' sample_size=300
#' R1 = 3
#' R2 = 2
#' rate = 0.5 #Proportion of value zero in beta
#' # Set true beta
#' zeroNum = round(rate*p_size)
#' ind = sample(1:p_size,zeroNum)
#' beta_true = runif(p_size,0,R2)
#' beta_true[ind] = 0
#' X = R1*matrix(rnorm(sample_size * p_size), ncol = p_size)
#' y=X%*%beta_true + rnorm(sample_size,mean=0,sd=2)
#' # Estimation
#' fit = GAGA(X,y,alpha = 3,family="gaussian")
#' Eb = fit$beta
#' #Create testing data
#' X_t = R1*matrix(rnorm(sample_size * p_size), ncol = p_size)
#' y_t=X_t%*%beta_true + rnorm(sample_size,mean=0,sd=2)
#' #Prediction
#' Ey = predict.GAGA(fit,newx=X_t)
#'
#' cat("\n err:", norm(Eb-beta_true,type="2")/norm(beta_true,type="2"))
#' cat("\n acc:", cal.w.acc(as.character(Eb!=0),as.character(beta_true!=0)))
#' cat("\n perr:", norm(Ey-y_t,type="2")/sqrt(sample_size))
#'
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
