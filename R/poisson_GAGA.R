
#' Fit a poisson model via the GAGA algorithm
#'
#' Fit a poisson model the Global Adaptive Generative Adjustment algorithm
#' @param X Input matrix, of dimension nobs*nvars; each row is an observation.
#' If the intercept term needs to be considered in the estimation process, then the first column of \code{X} must be all 1s.
#' In order to run the program stably, it is recommended that the value of X should not be too large. It is recommended to
#' preprocess all the items in X except the intercept item by means of preprocessing, so that the mean value of each column
#' is 0 and the standard deviation is \code{1/ colnum(X)}.
#' @param y Non-negative count response vector.
#' @param alpha Hyperparameter. The suggested value for alpha is 1 or 2.
#' When the collinearity of the load matrix is serious, the hyperparameters can be selected larger, such as 5.
#' @param itrNum The number of iteration steps. In general, 20 steps are enough.
#' If the condition number of \code{X} is large, it is recommended to greatly increase the
#' number of iteration steps.
#' @param flag It identifies whether to make model selection. The default is \code{TRUE}.
#' @param lamda_0 The initial value of the regularization parameter for ridge regression.
#' The running result of the algorithm is not sensitive to this value.
#' @param fdiag It identifies whether to use diag Approximation to speed up the algorithm.
#'
#'
#' @return Coefficient vector.
#' @export poisson_GAGA
#'
#' @examples
#' # Poisson
#' set.seed(2022)
#' library(mvtnorm)
#' p_size = 30
#' sample_size=300
#' R1 = 1/sqrt(p_size)
#' R2 = 5
#' rate = 0.5 #Proportion of value zero in beta
#' # Set true beta
#' zeroNum = round(rate*p_size)
#' ind = sample(1:p_size,zeroNum)#'
#' beta_true = runif(p_size,0,R2)
#' beta_true[ind] = 0
#' X = R1*matrix(rnorm(sample_size * p_size), ncol = p_size)
#' X[1:sample_size,1]=1
#' y = rpois(sample_size,lambda = as.vector(exp(X%*%beta_true)))
#' y = as.vector(y)
#' # Estimate
#' fit = GAGA(X,y,alpha = 2,family="poisson")
#' Eb = fit$beta
#' cat("\n err:", norm(Eb-beta_true,type="2")/norm(beta_true,type="2"))
#' cat("\n acc:", cal.w.acc(as.character(Eb!=0),as.character(beta_true!=0)))
#'
#'
poisson_GAGA = function(X,y,alpha=1,itrNum=30,flag=TRUE,lamda_0=0.5,fdiag=TRUE){

  vnames=colnames(X)
  if(is.null(vnames))vnames=paste("V",seq(ncol(X)),sep="")
  fit = list()
  class(fit) = c("GAGA","poisson")

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
  fit$beta = as.vector(beta)
  names(fit$beta) = vnames
  fit$alpha = alpha
  fit$itrNum = itrNum
  fit$fdiag = fdiag
  return(fit)
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
