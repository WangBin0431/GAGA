#Part of this function refers to the coxphfit function in MATLAB 2016b
#' Title
#'
#' @param X
#' @param t
#' @param alpha
#' @param itrNum
#' @param flag
#' @param lamda_0
#' @param fdiag
#'
#' @return
#' @export
#'
#' @examples
cox_GAGA = function(X,t,alpha=2,itrNum=20,flag=TRUE,lamda_0=0.5,fdiag=TRUE){
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
    if(index == itrNum){
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

    if(flag&&(index==itrNum)){
      tmpQ = sum(db<=100)
      if(tmpQ == 0){
        beta = rep(0,Q)
        return(as.vector(beta))
      }
      cov0 = getDDfu.cox(beta,X,sorty,freq,cens,atrisk,rep(0,p),fdiag)
      #browser()
      beta[E_pow_beta<diag(cov0)] = 0
      #beta[db>20] = 0;

    }else{
      b_old = b
      b = 1/inv_b
    }

  }#for (index in 1:itrNum)
  return(as.vector(beta))
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
