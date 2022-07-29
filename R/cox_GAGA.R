#' Fit a Cox model via the GAGA algorithm.
#'
#' Fit a Cox model via the Global Adaptive Generative Adjustment algorithm.
#' Part of this function refers to the coxphfit function in MATLAB 2016b.
#' @param X Input matrix, of dimension nobs*nvars; each row is an observation.
#' If the intercept term needs to be considered in the estimation process, then the first column of \code{X} must be all 1s.
#' @param t A n*2 matrix, one column should be named "time", indicating the survival time;
#' the other column must be named "status", and consists of 0 and 1, 0 indicates that the row of data is censored, 1 is opposite.
#' @param alpha Hyperparameter. The suggested value for alpha is 2 or 3.
#' @param itrNum The number of iteration steps. In general, 20 steps are enough.
#' If the condition number of \code{X} is large, it is recommended to greatly increase the
#' number of iteration steps.
#' @param thresh Convergence threshold for beta Change, if \code{max(abs(beta-beta_old))<threshold}, return.
#' @param flag It identifies whether to make model selection. The default is \code{TRUE}.
#' @param lamda_0 The initial value of the regularization parameter for ridge regression.
#' The running result of the algorithm is not sensitive to this value.
#' @param fdiag It identifies whether to use diag Approximation to speed up the algorithm.
#'
#' @return Coefficient vector.
#' @export cox_GAGA
#'
#' @examples
#' set.seed(2022)
#' library(mvtnorm)
#' p_size = 50
#' sample_size = 500
#' test_size = 1000
#' R1 = 3
#' R2 = 1
#' rate = 0.5 #Proportion of value zero in beta
#' censoringRate = 0.25 #Proportion of censoring data in observation data
#' # Set true beta
#' zeroNum = round(rate*p_size)
#' ind = sample(1:p_size,zeroNum)#'
#' beta_true = runif(p_size,-R2,R2)
#' beta_true[ind] = 0
#' # Generate training samples
#' cov_mat=matrix(1:p_size*p_size,p_size,p_size) ##covariance matrix
#' for(i in 1:p_size){for(j in 1:p_size) {if(i!=j){cov_mat[i,j]=0.0}else{cov_mat[i,j]=1}}}
#' X = R1*rmvnorm(n=sample_size, mean=rep(0, nrow(cov_mat)), sigma=cov_mat)
#' z = X%*%beta_true
#' u = runif(sample_size,0,1)
#' t = ((-log(1-u)/(3*exp(z)))*100)^(0.1)
#' cs = rep(0,sample_size)
#' csNum = round(censoringRate*sample_size)
#' ind = sample(1:sample_size,csNum)#'
#' cs[ind] = 1
#' t[ind] = runif(csNum,0,0.8)*t[ind]
#' y = cbind(t,1 - cs)
#' colnames(y) = c("time", "status")
#' #Estimation
#' fit = GAGA(X,y,alpha=2,family="cox")
#' Eb = fit$beta
#'
#' #Generate testing samples
#' X_t = R1*rmvnorm(n=test_size, mean=rep(0, nrow(cov_mat)), sigma=cov_mat)
#' z = X_t%*%beta_true
#' u = runif(sample_size,0,1)
#' t = ((-log(1-u)/(3*exp(z)))*100)^(0.1)
#' cs = rep(0,test_size)
#' csNum = round(censoringRate*test_size)
#' ind = sample(1:test_size,csNum)#'
#' cs[ind] = 1
#' t[ind] = runif(csNum,0,0.8)*t[ind]
#' y_t = cbind(t,1 - cs)
#' colnames(y_t) = c("time", "status")
#' #Prediction
#' pred = predict.GAGA(fit,newx=X_t)
#'
#' cat("\n err:", norm(Eb-beta_true,type="2")/norm(beta_true,type="2"))
#' cat("\n acc:", cal.w.acc(as.character(Eb!=0),as.character(beta_true!=0)))
#' cat("\n Cindex:", cal.cindex(pred,y_t))
#'
cox_GAGA = function(X,t,alpha=2,itrNum=20,thresh=1.e-3,flag=TRUE,lamda_0=0.5,fdiag=TRUE){

  exitflag = FALSE
  vnames=colnames(X)
  if(is.null(vnames))vnames=paste("V",seq(ncol(X)),sep="")
  fit = list()
  class(fit) = c("GAGA","cox")

  eps = 1.e-19
  n = nrow(X)
  p = ncol(X)
  y = t[,"time"]
  tmp = sort(y,index.return = T)
  sorty = tmp$x
  idx = tmp$ix
  X = X[idx,]
  cens = 1-t[idx,"status"]

  #In the future, "freq" can be used to expand the situation that each observation frequency is not equal
  freq = rep(1,n)
  sumf = max(1, sum(freq))
  baseX = (t(freq)%*%X) / sumf
  X = X - t(matrix(rep(baseX,n),ncol = n))

  atrisk = mrank(sorty)

  b = rep(1,p)*lamda_0
  b_old = b

  for (index in 1:itrNum){
    if(index == itrNum || exitflag){
      db = (b-b_old)
      b = b/alpha
    }
    b = as.vector(b)
    if(index==1){
      stdX = t(sqrt(t(freq)%*%(X*X)/sumf))
      beta = rep(0,p)
      t = stdX!=0
      beta[t] = 0.01/stdX[t]
    }

    bfgsItr = 20
    tmpout  = getEb.cox(beta,X,sorty,freq,cens,atrisk,b,bfgsItr,fdiag)
    beta = tmpout$beta
    cov_beta = tmpout$gg
    E_pow_beta = diag(cov_beta) + beta*beta

    inv_b = E_pow_beta/alpha

    if(flag&&(index==itrNum || exitflag)){
      tmpQ = sum(db<=100)
      if(tmpQ == 0){
        beta = rep(0,Q)
        fit$beta = as.vector(beta)
        names(fit$beta) = vnames
        fit$alpha = alpha
        fit$itrNum = index
        fit$fdiag = fdiag
        return(fit)
      }
      cov0 = getDDfu.cox(beta,X,sorty,freq,cens,atrisk,rep(0,p),fdiag)
      #browser()
      beta[E_pow_beta<diag(cov0)] = 0
      #beta[db>20] = 0;
      break

    }else{
      b_old = b
      b = 1/inv_b
    }

    if(index==1){
      beta_old = beta
    }else{
      if(max(abs(beta-beta_old))<thresh)exitflag = TRUE
      beta_old = beta
    }

  }#for (index in 1:itrNum)
  fit$beta = as.vector(beta)
  names(fit$beta) = vnames
  fit$alpha = alpha
  fit$itrNum = index
  fit$fdiag = fdiag
  return(fit)
}

################################################################################
mrank = function(sorty){
  atrisk = rep(0,length(sorty))
  atrisk[1] = 1
  if(length(sorty)==1){
    return(atrisk)
  }

  for(ii in 2:length(sorty)){
    if(sorty[ii]==sorty[ii-1]){
      atrisk[ii] = atrisk[ii-1]
    }else{
      atrisk[ii] = ii
    }
  }
  return(atrisk)
}


#################################################################################
getEb.cox  = function(u0,X,sorty,freq,cens,atrisk,b,bfgsItr,fdiag){
  n = nrow(X)
  p = ncol(X)
  tmpout = negloglike.cox(u0,X,sorty,freq,cens,atrisk,b,fdiag)
  v = tmpout$v
  g = tmpout$g
  D = tmpout$gg
  u = u0
  for(ii in 1:bfgsItr){
    d = -D%*%g
    LL = optimize(Func_lambda.cox, c(0, 2), tol = 1.e-19, u=u, X=X, sorty=sorty, freq=freq, cens=cens, atrisk=atrisk, b=b,d=d)
    u = u + LL$minimum*d
    tmpout = negloglike.cox(u,X,sorty,freq,cens,atrisk,b,fdiag)
    v = tmpout$v
    g = tmpout$g
    D = tmpout$gg

    if(norm(g,type="2")/p<1.e-16){
      break;
    }
  }
  tmpout$beta = u
  return(tmpout)
}

Func_u.cox = function(u,X,sorty,freq,cens,atrisk,b){
  obsfreq = freq*(1-cens)
  Xu = X%*%u
  r = exp(Xu)
  risksum = spatstat.utils::revcumsum(freq*r) #rev(cumsum(rev(freq*r)))
  risksum = risksum[atrisk]
  v = -t(obsfreq)%*%(Xu - log(risksum)) + 0.5*sum(b*(u^2))
  return(v)
}

Func_lambda.cox = function(lambda,u,X,sorty,freq,cens,atrisk,b,d){
  return(Func_u.cox(u+lambda*d,X,sorty,freq,cens,atrisk,b))
}

getDDfu.cox = function(u,X,sorty,freq,cens,atrisk,b,fdiag){
  obsfreq = freq*(1-cens)
  Xu = X%*%u
  r = exp(Xu)
  tmp1 = as.numeric(r)*freq
  risksum = spatstat.utils::revcumsum(tmp1) #rev(cumsum(rev(freq*r)))
  risksum = risksum[atrisk]
  n = nrow(X)
  p = ncol(X)
  Xr = X*tmp1
  Xrsum = apply(Xr,2,spatstat.utils::revcumsum) #apply(apply(apply(Xr,2,rev),2,cumsum),2,rev)
  Xrsum = Xrsum[atrisk,]
  A = Xrsum/risksum

  if(fdiag == FALSE){
    t1 = rep(1:p,p)
    t2 = sort(t1)
    XXr = X[,t1]*X[,t2]*(tmp1)
    XXrsum = apply(XXr,2,spatstat.utils::revcumsum)  #apply(apply(apply(XXr,2,rev),2,cumsum),2,rev)
    XXrsum = XXrsum[atrisk,]/risksum
    gg = t(obsfreq)%*%XXrsum
    dim(gg) = c(p,p)
    gg = MASS::ginv(gg - t(A)%*%(A*obsfreq) + diag(b))
    return(gg)
  }else{
    XXr = X*X*tmp1
    XXrsum = apply(XXr,2,spatstat.utils::revcumsum)
    XXrsum = XXrsum[atrisk,]/risksum
    gg = t(obsfreq)%*%XXrsum
    gg = diag(as.numeric(1/(gg - apply(A*A*obsfreq,2,sum)+b)))
    return(gg)
  }

}

negloglike.cox = function(u,X,sorty,freq,cens,atrisk,b,fdiag){
  tmpout = list()
  obsfreq = freq*(1-cens)
  Xu = X%*%u
  r = exp(Xu)
  risksum = spatstat.utils::revcumsum(freq*r) #rev(cumsum(rev(freq*r)))
  risksum = risksum[atrisk]
  tmpout$v = -t(obsfreq)%*%(Xu - log(risksum)) + 0.5*sum(b*(u^2))

  n = nrow(X)
  p = ncol(X)

  tmp1 = as.numeric(r)*freq
  Xr = X*tmp1
  Xrsum = apply(Xr,2,spatstat.utils::revcumsum) #apply(apply(apply(Xr,2,rev),2,cumsum),2,rev)
  Xrsum = Xrsum[atrisk,]
  A = Xrsum/risksum
  tmpout$g = t(-t(obsfreq)%*%(X-A)) + (u*b)


  if(fdiag == FALSE){
    t1 = rep(1:p,p)
    t2 = sort(t1)
    XXr = X[,t1]*X[,t2]*(tmp1)
    XXrsum = apply(XXr,2,spatstat.utils::revcumsum)  #apply(apply(apply(XXr,2,rev),2,cumsum),2,rev)
    XXrsum = XXrsum[atrisk,]/risksum
    gg = t(obsfreq)%*%XXrsum
    dim(gg) = c(p,p)
    gg = MASS::ginv(gg - t(A)%*%(A*obsfreq) + diag(b))
  }else{
    XXr = X*X*tmp1
    XXrsum = apply(XXr,2,spatstat.utils::revcumsum)
    XXrsum = XXrsum[atrisk,]/risksum
    gg = t(obsfreq)%*%XXrsum
    gg = diag(as.numeric(1/(gg - apply(A*A*obsfreq,2,sum)+b)))
  }
  tmpout$gg = gg
  return(tmpout)
}
