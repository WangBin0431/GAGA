#' Fit a multinomial model via the GAGA algorithm
#'
#' Fit a multinomial model the Global Adaptive Generative Adjustment algorithm
#' @param X Input matrix, of dimension nobs*nvars; each row is an observation.
#' If the intercept term needs to be considered in the estimation process, then the first column of \code{X} must be all 1s.
#' @param y a One-hot response matrix or a \code{nc>=2} level factor
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
#' @return Coefficient matrix with K-1 columns beta_1,...,beta_{K-1} where K is the class number.
#' For k=1,..,K-1, the probability
#' \deqn{Pr(G=k|x)=exp(x^T beta_k) /(1+sum_{k=1}^{K-1}exp(x^T beta_k))}.
#' For k=K, the probability \deqn{Pr(G=K|x)=1/(1+sum_{k=1}^{K-1}exp(x^T beta_k))}.
#' @export multinomial_GAGA
#'
#' @examples
#'# multinomial
#' set.seed(2022)
#' library(mvtnorm)
#' cat("\n")
#' cat("Test multinomial GAGA\n")#'
#' p_size = 20
#' C = 3
#' classnames = c("C1","C2","C3","C4")
#' sample_size = 500
#' test_size = 1000#'
#' rate = 0.5 #Proportion of value zero in beta
#' Num = 10 # Total number of experiments
#' R1 = 1
#' R2 = 5#'
#' #Set true beta
#' beta_true = matrix(rep(0,p_size*C),c(p_size,C))
#' zeroNum = round(rate*p_size)
#' for(jj in 1:C){
#'   ind = sample(1:p_size,zeroNum)#'
#'   tmp = runif(p_size,0,R2)
#'   tmp[ind] = 0
#'   beta_true[,jj] = tmp
#' }#'
#' cov_mat=matrix(1:p_size*p_size,p_size,p_size)
#' for(i in 1:p_size){for(j in 1:p_size) {if(i!=j){cov_mat[i,j]=0.0}else{cov_mat[i,j]=1}}}#'
#' #Generate training samples
#' X = R1*rmvnorm(n=sample_size, mean=rep(0, nrow(cov_mat)), sigma=cov_mat)#'
#' X[1:sample_size,1]=1
#' z = X%*%beta_true
#' t = exp(z)/(1+rowSums(exp(z)))
#' t = cbind(t,1-rowSums(t))
#' tt = t(apply(t,1,cumsum))
#' tt = cbind(rep(0,sample_size),tt)#'
#' # y = matrix(rep(0,sample_size*(C+1)),c(sample_size,C+1))
#' y = rep(0,sample_size)
#' for(jj in 1:sample_size){
#'   tmp = runif(1,0,1)
#'   for(kk in 1:(C+1)){
#'     if((tmp>tt[jj,kk])&&(tmp<=tt[jj,kk+1])){
#'       # y[jj,kk] = 1
#'       y[jj] = kk
#'       break
#'     }
#'   }
#' }
#' y = classnames[y]
#' # Estimation
#' fit = GAGA(X, y,alpha=1,family = "multinomial")
#' Eb = fit$beta#'
#' #Prediction
#' #Generate test samples
#' X_t = R1*rmvnorm(n=test_size, mean=rep(0, nrow(cov_mat)), sigma=cov_mat)
#' X_t[1:test_size,1]=1
#' z = X_t%*%beta_true
#' t = exp(z)/(1+rowSums(exp(z)))
#' t = cbind(t,1-rowSums(t))
#' tt = t(apply(t,1,cumsum))
#' tt = cbind(rep(0,test_size),tt)
#'
#' y_t = rep(0,test_size)
#' #y_tt = rep(0,test_size)
#' for(jj in 1:test_size){
#'   tmp = runif(1,0,1)
#'   for(kk in 1:(C+1)){
#'     if((tmp>tt[jj,kk])&&(tmp<=tt[jj,kk+1])){
#'       y_t[jj] = kk
#'       break
#'     }
#'   }
#'   #y_tt[jj] = which.max(t[jj,])
#' }
#' y_t = classnames[y_t]#
#' Ey = predict.GAGA(fit,newx = X_t)
#
#' cat("\n--------------------")
#' cat("\n err:", norm(Eb-beta_true,type="2")/norm(beta_true,type="2"))
#' cat("\n acc:", cal.w.acc(as.character(Eb!=0),as.character(beta_true!=0)))
#' cat("\n pacc:", cal.w.acc(as.character(Ey),as.character(y_t)))
#' cat("\n")

multinomial_GAGA = function(X,y,alpha=2,itrNum=50,flag=TRUE,lamda_0=0.001,fdiag=TRUE){
  #browser()
  eps = 1.e-19
  N = nrow(X)
  P = ncol(X)

  # refers to "lognet.R" of glmnet
  vnames=colnames(X)
  if(is.null(vnames))vnames=paste("V",seq(P),sep="")

  if(!is.matrix(y)){
    y = as.factor(y)
    ntab = table(y)
    classnames=names(ntab)
    nc = as.integer(length(ntab))
    y = diag(nc)[as.numeric(y), ,drop=FALSE]
  }else{
    if(!(ncol(y)>=2&&min(y)==0&&max(y)==1))stop("y should be a one hot matrix or a nc>=2 level factor")
    classnames = colnames(y)
    if(is.null(classnames))classnames=paste("V",seq(ncol(y)),sep="")
  }

  K = ncol(y)
  C = K-1

  fit = list()
  fit$classnames = classnames
  class(fit) = c("GAGA","multinomial")

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
  fit$beta = beta
  rownames(fit$beta) = vnames
  colnames(fit$beta) = classnames[1:C]
  fit$alpha = alpha
  fit$itrNum = itrNum
  fit$fdiag = fdiag
  return(fit)
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
