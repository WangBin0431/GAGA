
#' Fit a linear model via the GAGA algorithm
#'
#' Fit a linear model with a Gaussian noise via the Global Adaptive Generative Adjustment algorithm
#' @param X Input matrix, of dimension nobs*nvars; each row is an observation.
#' If the intercept term needs to be considered in the estimation process, then the first column of \code{X} must be all 1s.
#' @param y Quantitative response vector.
#' @param alpha Hyperparameter. The suggested value for alpha is 2 or 3.
#' When the collinearity of the load matrix is serious, the hyperparameters can be selected larger, such as 5.
#' @param itrNum The number of iteration steps. In general, 20 steps are enough.
#' If the condition number of \code{X} is large, it is recommended to greatly increase the
#' number of iteration steps.
#' @param QR_flag It identifies whether to use QR decomposition to speed up the algorithm.
#' Currently only valid for linear models.
#' @param flag It identifies whether to make model selection. The default is \code{TRUE}.
#' @param lamda_0 The initial value of the regularization parameter for ridge regression.
#' The running result of the algorithm is not sensitive to this value.
#' @param fix_sigma It identifies whether to update the variance estimate of the Gaussian noise or not.
#' \code{fix_sigma=TRUE} uses the initial variance as the variance estimate in each loop.
#' \code{fix_sigma=FALSE} updates the variance estimate in each loop.
#' @param sigm2_0 The initial variance of the Gaussian noise.
#'
#' @return Coefficient vector.
#' @export LM_GAGA
#'
#' @examples
#' # Gaussian
#' set.seed(2022)
#' p_size = 30
#' sample_size=300
#' R1 = 3
#' R2 = 2
#' rate = 0.5 #Proportion of value zero in beta
#' # Set true beta
#' zeroNum = round(rate*p_size)
#' ind1 = sample(1:p_size,p_size)
#' ind2 = ind1[1:zeroNum]
#' beta_true = runif(p_size,0,R2)
#' X = R1*matrix(rnorm(sample_size * p_size), ncol = p_size)
#' y=X%*%beta_true + rnorm(sample_size,mean=0,sd=2)
#' # Estimate
#' Eb = GAGA(X,y,alpha = 3,family="gaussian")
#' cat("\n err:", norm(Eb-beta_true,type="2")/norm(beta_true,type="2"))
#' cat("\n acc:", cal.w.acc(as.character(Eb!=0),as.character(beta_true!=0)))
#'
#' # Binomial
#' set.seed(2022)
#' library(mvtnorm)
#' p_size = 30
#' sample_size=600
#' R1 = 1
#' R2 = 3
#' rate = 0.5 #Proportion of value zero in beta
#' # Set true beta
#' zeroNum = round(rate*p_size)
#' ind1 = sample(1:p_size,p_size)
#' ind2 = ind1[1:zeroNum]
#' beta_true = runif(p_size,R2*0.2,R2)
#' beta_true[ind2] = 0
#' X = R1*matrix(rnorm(sample_size * p_size), ncol = p_size)
#' X[1:sample_size,1]=1
#' t = 1/(1+exp(-X%*%beta_true))
#' tmp = runif(sample_size,0,1)
#' y = rep(0,sample_size)
#' y[t>tmp] = 1
#' # Estimate
#' Eb = GAGA(X,y,family = "binomial", alpha = 1)
#' cat("\n err:", norm(Eb-beta_true,type="2")/norm(beta_true,type="2"))
#' cat("\n acc:", cal.w.acc(as.character(Eb!=0),as.character(beta_true!=0)))
#'
LM_GAGA = function(X,y,alpha=3,itrNum=50,QR_flag=FALSE,flag=TRUE,lamda_0=0.001,fix_sigma=FALSE,sigm2_0 = 1){
  # print("This is LM_GAGA")
  # return(1)
  eps = 1.e-19

  if(!QR_flag){
    N = nrow(X);
    Q = ncol(X);

    if(fix_sigma)
      sigm2=sigm2_0
    else
      sigm2 = 1;

    B = diag(lamda_0*rep(1,Q))
    B_old = B

    Xty = t(X)%*%y
    XtX = t(X)%*%X
    yty = t(y)%*%y

    for (index in 1:itrNum) {

      if(sigm2 == 0)
        sigm2 = eps

      if(index==itrNum){
        dB = diag(B-B_old)
        B = B/alpha
      }
      beta = getEb.LM(XtX,Xty,sigm2,B)

      cov_beta = getEbb.LM(XtX,sigm2,B)

      E_pow_beta_M = cov_beta + beta%*%t(beta)

      E_pow_beta = diag(E_pow_beta_M)

      inv_b = E_pow_beta/alpha

      if(flag&&(index==itrNum)){
        tmpQ = sum(dB<=100)
        if(tmpQ == 0){
          beta = rep(0,Q)
          return(as.vector(beta))
        }
        if(tmpQ<=N){
          X_sub = X[,dB<=100]
          CRB_sub = diag(getEbb.LM(t(X_sub)%*%X_sub,sigm2,0))
          beta_sub = beta[dB<=100]
          E_pow_beta_sub = E_pow_beta[dB<=100];
          beta_sub[E_pow_beta_sub<CRB_sub] = 0;
          beta[dB<=100] = beta_sub;
          beta[dB>100] = 0;
          return(as.vector(beta))
        }
        else{
          pow_beta = beta^2
          aa = sort(pow_beta,T,index.return = T)
          idx = aa$ix

          sortbeta = beta[idx]
          sortbeta_sub = sortbeta[1:N]
          sortX = X[,idx]
          sortX_sub = sortX[,1:N]
          sortCRB_sub = diag(getEbb.LM(t(sortX_sub)%*%sortX_sub,sigm2,0))
          E_pow_beta_sort = E_pow_beta[idx]
          E_pow_beta_sort_sub = E_pow_beta_sort[1:N]
          sortbeta_sub[E_pow_beta_sort_sub<sortCRB_sub] = 0
          sortbeta[1:N] = sortbeta_sub
          sortbeta[N+1:length(sortbeta)] = 0
          beta[idx] = sortbeta
          return(as.vector(beta))
        }
      }#if(flag&&(index==itrNum))

      B_old = B
      B = diag(1/inv_b)

      if(fix_sigma){
        sigm2=sigm2_0
      }
      else{
        #sigm2 = (yty-2*t(beta)%*%Xty)/N + sum(diag(E_pow_beta_M%*%XtX))/N
        sigm2 = (yty-2*t(beta)%*%Xty)/N + sum(E_pow_beta_M*XtX)/N
        sigm2 = as.vector(sigm2)
      }

    }#for (index in 1:itrNum)
  }#if(!QR_flag)
  else{
    # use QR
    yty = t(y)%*%y

    N = nrow(X);
    Q = ncol(X);
    if(Q>N){
      beta = OMP(X,y,N)
      aa = sort(abs(beta)^2,T,index.return = T)
      idx1 = aa$ix
      X = X[,idx1[1:N]]
    #browser()
  }

  beta = solve(t(X)%*%X + lamda_0*diag(rep(1,ncol(X))),tol = 0)%*%t(X)%*%y

  pow_beta = abs(beta)^2
  aa = sort(pow_beta,T,index.return = T)
  idx = aa$ix
  sortX = X[,idx]

  tmp = qr(sortX)
  K = qr.Q(tmp)
  R = qr.R(tmp)

  tmp = diag(diag(R))

  R = solve(tmp,tol = 0)%*%R
  K = K%*%tmp

  N_K = nrow(K)
  Q_K = ncol(K)
  N_R = nrow(R)
  Q_R = ncol(R)

  if(fix_sigma)
    sigm2=sigm2_0
  else
    sigm2 = 1;


  b = lamda_0*rep(1,Q_K)
  b_old = b

  kk = diag(t(K)%*%K)
  Kty = t(K)%*%y

  for(index in 1:itrNum){
    if(sigm2 == 0){
      sigm2 = eps
    }
    if(index == itrNum){
      dB = b - b_old
      b = b/alpha
    }
    EW = getEb_orth.LM(kk,Kty,sigm2,b)

    covW = getEbb_orth.LM(kk,sigm2,b)

    EWW = covW + abs(EW)^2

    if(flag&&(index==itrNum)&&(sum(dB<=100)!=0)){

      K_sub = K[,dB<=100]
      kk_sub = diag(t(K_sub)%*%K_sub)
      CRB_sub = getEbb_orth.LM(kk_sub,sigm2,0)
      EW_sub = EW[dB<=100]
      EWW_sub = EWW[dB<=100]
      EW_sub[EWW_sub<CRB_sub] = 0
      EW[dB<=100] = EW_sub
      EW[dB>100] = 0
    }

    b_old = b
    b = alpha/EWW

    if(fix_sigma){
      sigm2=sigm2_0
    }
    else{
      sigm2 = (yty - 2*t(y)%*%K%*%EW)/N + sum(EWW*kk)/N
      sigm2 = as.vector(sigm2)
    }


  }#for(index in 1:itrNum)
  #browser()
  beta = solve(R,tol = 0)%*%EW
  beta[idx] = beta

  if(Q>N){
    tmp = rep(0,Q)
    tmp[1:N] = beta;
    tmp[idx1] = tmp
    beta = tmp
  }
  return(as.vector(beta))
  }#QR_flag == TRUE
}


getEb.LM = function(XtX,Xty,sigma2,W2){
  #browser()
  solve(XtX + sigma2*W2,tol = 0)%*%Xty
}
getEbb.LM = function(XtX,sigma2,W2){
  solve(XtX + sigma2*W2,tol = 0)*sigma2
}
getEb_orth.LM = function(kk,Kty,sigma2,b){
  Kty/(kk+sigma2*b)
}
getEbb_orth.LM = function(kk,sigma2,b){
  1/(kk/sigma2 + b)
}

#################################################################################
OMP = function(X,y,L,eps = 0){
  n = nrow(X)
  K = ncol(X)
  Xt = t(X)
  a = NULL
  residual = y
  indx = rep(0,L)

  for(j in 1:L){
    proj = Xt%*%residual
    pos = which.max(abs(proj))
    indx[j] = pos
    #browser()
    a = ginv(X[,indx[1:j]])%*%y
    residual = y - X[,indx[1:j]]%*%a

    if(sum(residual^2)<=eps){
      break
    }
  }#for(k in 1:L)
  temp = rep(0,K)
  temp[indx[1:j]] = a;
  return(temp)
}

