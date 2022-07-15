
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
poisson_GAGA = function(X,y,alpha=1,itrNum=30,flag=TRUE,lamda_0=0.5,fdiag=TRUE){
  # print("This is poisson_GAGA")
  # return(2)
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
    b = as.vector(b)
    if(index==1){
      tmpy = y
      tmpy[tmpy<=0] = 0.1;
      beta = solve(t(X)%*%X + diag(b),tol = 0)%*%(t(X)%*%log(tmpy))
      cov_beta = getDDfu.poisson(beta,X,y,b,fdiag)
      D0 = cov_beta
    }

    bfgsItr = 20
    beta  = getEb.poisson(X,y,b,beta,D0,bfgsItr,fdiag)
    cov_beta = getDDfu.poisson(beta,X,y,b,fdiag)

    D0 = cov_beta

    E_pow_beta = diag(cov_beta) + beta*beta

    inv_b = E_pow_beta/alpha

    if(flag&&(index==itrNum)){
      tmpQ = sum(db<=100)
      if(tmpQ == 0){
        beta = rep(0,p)
        return(as.vector(beta))
      }
      cov0 = getDDfu.poisson(beta,X,y,rep(0,p),fdiag)
      beta[E_pow_beta<diag(cov0)] = 0
      beta[db>20] = 0;

    }else{
      b_old = b
      b = 1/inv_b
    }

  }#for (index in 1:itrNum)
  return(as.vector(beta))
}


#################################################################################
getEb.poisson  = function(X,y,b,beta,D0,bfgsItr,fdiag){
  n = nrow(X)
  p = ncol(X)
  D = D0
  u = beta
  g = Dfu.poisson(u,X,y,b)
  for(ii in 1:bfgsItr){
    d = -D%*%g
    LL = optimize(Func_lambda.poisson, c(0, 2), tol = 1.e-19, u=u, X=X, y=y, b=b,d=d)
    #browser()
    u = u + LL$minimum*d
    D = getDDfu.poisson(u,X,y,b,fdiag)
    g = Dfu.poisson(u,X,y,b)
    if(norm(g,type="2")/p<1.e-16){
      break;
    }
  }
  return(u)
}





Func_u.poisson = function(u,X,y,b){

  eps = 2.2204e-16
  Xu = X%*%u
  v = -t(y)%*%Xu + sum(exp(Xu)) + 0.5*sum(b*(u^2))
  return(v)
}

Func_lambda.poisson = function(lambda,u,X,y,b,d){
  return(Func_u.poisson(u+lambda*d,X,y,b))
}

Dfu.poisson = function(u,X,y,b){
  Xu = X%*%u
  tmp = exp(Xu) - y
  g = (t(X)%*%tmp)+(u*b)
  return(g)
}

getDDfu.poisson = function(u,X,y,b,fdiag){
  if(fdiag == FALSE){
    gg = MASS::ginv(t(X)%*%(X*as.vector(exp(X%*%u)))+diag(b))
    return(gg)
  }else{
    gg = diag(1/diag(t(X)%*%(X*as.vector(exp(X%*%u)))+diag(b)))
    return(gg)
  }
}
