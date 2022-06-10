library(MASS)
GaGa = function(X,y,ratio=3,itrNum=50,QR_flag=FALSE,flag=TRUE,lamda_0=0.001,fix_sigma=FALSE,sigm2_0 = 1){
  
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
      #cat("index is ",index,"\n")
      #if(index == 35)
      #  browser()
      
      if(sigm2 == 0)
        sigm2 = eps
      
      if(index==itrNum){
        dB = diag(B-B_old)
        B = B/ratio
      }
      beta = getEW(XtX,Xty,sigm2,B)
      
      cov_beta = getEWW(XtX,sigm2,B)
      
      E_pow_beta_M = cov_beta + beta%*%t(beta)
      
      E_pow_beta = diag(E_pow_beta_M)
      
      inv_b = E_pow_beta/ratio
      
      if(flag&&(index==itrNum)){
        tmpQ = sum(dB<=100)
        if(tmpQ == 0){
          beta = rep(0,Q)
          return(as.vector(beta))
        }
        if(tmpQ<=N){
          X_sub = X[,dB<=100]
          CRB_sub = diag(getEWW(t(X_sub)%*%X_sub,sigm2,0))
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
          sortCRB_sub = diag(getEWW(t(sortX_sub)%*%sortX_sub,sigm2,0))
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
        b = b/ratio
      }
      EW = getEW_orth(kk,Kty,sigm2,b)
      
      covW = getEWW_orth(kk,sigm2,b)
      
      EWW = covW + abs(EW)^2
      
      if(flag&&(index==itrNum)&&(sum(dB<=100)!=0)){
        
        K_sub = K[,dB<=100]
        kk_sub = diag(t(K_sub)%*%K_sub)
        CRB_sub = getEWW_orth(kk_sub,sigm2,0)
        EW_sub = EW[dB<=100]
        EWW_sub = EWW[dB<=100]
        EW_sub[EWW_sub<CRB_sub] = 0
        EW[dB<=100] = EW_sub
        EW[dB>100] = 0
      }
      
      b_old = b
      b = ratio/EWW
      
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


getEW = function(XtX,Xty,sigma2,W2){
  #browser()
  solve(XtX + sigma2*W2,tol = 0)%*%Xty
}
getEWW = function(XtX,sigma2,W2){
  solve(XtX + sigma2*W2,tol = 0)*sigma2
}
getEW_orth = function(kk,Kty,sigma2,b){
  Kty/(kk+sigma2*b)
}
getEWW_orth = function(kk,sigma2,b){
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

