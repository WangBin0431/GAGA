
#' Title
#'
#' @param X
#' @param y
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
logistic_GAGA = function(X,y,alpha=1,itrNum=30,flag=TRUE,lamda_0=0.001,fdiag=TRUE){
  eps = 1.e-19
  n = nrow(X)
  p = ncol(X)
  b = rep(1,p)*lamda_0
  b_old = b

  for (index in 1:itrNum){
    if(index == itrNum){
      db = (b-b_old)
      b = b/alpha
    }

    if(index==1){
      beta = rep(0,p)
      cov_beta = getDDfu.logistic(beta,X,y,b,fdiag)
      D0 = cov_beta
    }
    b = as.vector(b)
    bfgsItr = 20
    #browser()
    beta  = getEb.logistic(X,y,b,beta,D0,bfgsItr,fdiag)
    cov_beta = getDDfu.logistic(beta,X,y,b,fdiag)

    D0 = cov_beta

    E_pow_beta = diag(cov_beta) + beta*beta

    inv_b = E_pow_beta/alpha

    if(flag&&(index==itrNum)){
      cov0 = getDDfu.logistic(beta,X,y,rep(0,p),fdiag)
      beta[E_pow_beta<diag(cov0)] = 0
    }else{
      b_old = b
      b = 1/inv_b
    }
  }#for (index in 1:itrNum)

  return(as.vector(beta))
}
##############################################################################

#################################################################################
getEb.logistic  = function(X,y,b,beta,D0,bfgsItr,fdiag){
  n = nrow(X)
  p = ncol(X)
  #bfgs
  D = D0
  u = beta
  g = Dfu.logistic(u,X,y,b)
  for(ii in 1:bfgsItr){
    d = -D%*%g
    LL = optimize(Func_lambda.logistic, c(0, 2), tol = 1.e-19, u=u, X=X, y=y, b=b,d=d)
    #browser()
    u = u + LL$minimum*d
    D = getDDfu.logistic(u,X,y,b,fdiag)
    g = Dfu.logistic(u,X,y,b)
    if(norm(g,type="2")/p<1.e-16){
      break;
    }
  }
  return(u)
}



sigmod = function(x){1/(1+exp(-1*x))}

Func_u.logistic = function(u,X,y,b){

  eps = 2.2204e-16
  Xu = X%*%u
  tmp1 = y*log(sigmod(Xu)+eps)
  tmp2 = (1-y)*log(1-sigmod(Xu)+eps)
  v = -1*sum(tmp1 + tmp2) + 0.5*sum(b*(u^2))
  return(v)
}

Func_lambda.logistic = function(lambda,u,X,y,b,d){
  return(Func_u.logistic(u+lambda*d,X,y,b))
}

Dfu.logistic = function(u,X,y,b){
  Xu = X%*%u
  tmp = y - sigmod(Xu)
  g = -1*(t(X)%*%tmp)+(u*b)
  return(g)
}

getDDfu.logistic = function(u,X,y,b,fdiag){
  if(fdiag==FALSE){
    h = sigmod(X%*%u)
    #gg = solve(t(X)%*%diag(as.vector(h*(1-h)))%*%X+diag(b),tol = 0)
    #gg = ginv(t(X)%*%diag(as.vector(h*(1-h)))%*%X+diag(b))
    gg = MASS::ginv(t(X)%*%(X*as.vector(h*(1-h)))+diag(b))
    return(gg)
  }else{
    h = sigmod(X%*%u)
    gg = diag(1/diag(t(X)%*%(X*as.vector(h*(1-h)))+diag(b)))
    return(gg)
  }
}



