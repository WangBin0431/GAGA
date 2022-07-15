
#' Title
#'
#' @param X design matrix
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
multinomial_GAGA = function(X,y,alpha=2,itrNum=50,flag=TRUE,lamda_0=0.001,fdiag=TRUE){
  eps = 1.e-19
  N = nrow(X)
  P = ncol(X)
  K = ncol(y)
  C = K-1

  b = rep(1,P*C)*lamda_0
  b_old = b

  for (index in 1:itrNum){
    if(index == itrNum){
      db = (b-b_old)
      b = b/alpha
    }

    if(index==1){
      beta = rep(0,P*C)
      dim(beta) = c(P,C)
      cov_beta = getDDfu.multinomial(beta,X,y,matrix(b,c(P,C)),fdiag)
      D0 = cov_beta
    }
    b = as.vector(b)
    bfgsItr = 20
    beta  = getEb.multinomial(X,y,matrix(b,c(P,C)),beta,D0,bfgsItr,fdiag)
    cov_beta = getDDfu.multinomial(beta,X,y,matrix(b,c(P,C)),fdiag)

    D0 = cov_beta

    E_pow_beta = diag(cov_beta) + as.vector(beta*beta)

    inv_b = E_pow_beta/alpha

    if(flag&&(index==itrNum)){
      cov0 = getDDfu.multinomial(beta,X,y,matrix(rep(0,P*C),c(P,C)),fdiag)

      beta[E_pow_beta<diag(cov0)] = 0
      # beta[db>20] = 0;

    }else{
      b_old = b
      b = 1/inv_b
      sigm2 = 1
    }

  }#for (index in 1:itrNum)
  return(beta)
}


#################################################################################
getEb.multinomial  = function(X,y,b,beta,D0,bfgsItr,fdiag){
  N = nrow(X)
  P = ncol(X)
  D = D0
  u = beta
  g = Dfu.multinomial(u,X,y,b)
  for(ii in 1:bfgsItr){
    d = -D%*%as.vector(g)
    dim(d) = dim(beta)
    LL = optimize(Func_lambda.multinomial, c(0, 2), tol = 1.e-19, u=u, X=X, y=y, b=b,d=d)
    u = u + LL$minimum*d
    D = getDDfu.multinomial(u,X,y,b,fdiag)
    g = Dfu.multinomial(u,X,y,b)
    if(norm(g,type="2")/P<1.e-16){
      break;
    }
  }
  return(u)
}


sigmod = function(x){1/(1+exp(-1*x))}

Func_u.multinomial = function(u,X,y,b){
  eps = 2.2204e-16
  N = nrow(y)
  P = nrow(u)
  C = ncol(u)
  Xu = X%*%u
  v = -sum(y[1:N,1:C]*Xu) + sum(log(1+rowSums(exp(Xu)))) + 0.5*sum(b*u*u)
  return(v)
}

Func_lambda.multinomial = function(lambda,u,X,y,b,d){
  return(Func_u.multinomial(u+lambda*d,X,y,b))
}

Dfu.multinomial = function(u,X,y,b){
  Xu = X%*%u
  P = nrow(u)
  C = ncol(u)
  Xt = t(X)
  tmp1 = exp(Xu)
  tmp2 = 1/(1+rowSums(tmp1))
  g = Xt%*%(-y[,1:C] + (tmp2*tmp1)) + u*b
  return(g)
}

getDDfu.multinomial = function(u,X,y,b,fdiag){
  #browser()
  Xu = X%*%u
  P = nrow(u)
  C = ncol(u)
  Xt = t(X)
  invgg = rep(0,P*C*P*C)
  dim(invgg) = c(P*C,P*C)
  tmp2 = 1+rowSums(exp(Xu))
  tmp22 = tmp2^2
  for(ii in 1:C){
    indexii = (ii-1)*P+1
    for(jj in 1:C){
      indexjj = (jj-1)*P+1
      if(ii==jj){
        tmp1 = exp(Xu[,jj])
        tmp4 = tmp1*(tmp2-tmp1)/tmp22
        invgg[indexjj:(indexjj+P-1),indexii:(indexii+P-1)] = Xt%*%(tmp4*X) + diag(b[,jj])
      }else{
        if(fdiag == FALSE){
          #browser()
          tmp1 = exp(Xu[,jj]+Xu[,ii])
          tmp3 = Xt%*%(X*(tmp1/tmp22))
          invgg[indexjj:(indexjj+P-1),indexii:(indexii+P-1)] = tmp3
          invgg[indexii:(indexii+P-1),indexjj:(indexjj+P-1)] = tmp3
        }
      }
    }#for(jj in 1:C)
  }#for(ii in 1:C)
  if(fdiag==TRUE){
    gg = diag(1/diag(invgg))
  }else{
    #gg = ginv(invgg)
    gg = solve(invgg,tol = 0)
  }
  return(gg)
}




