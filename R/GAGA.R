
#' Fit a GLM with GAGA algorithm
#'
#' Fits linear, logistic and multinomial, poisson, and Cox regression models via Global Adaptive Generative Adjustment algorithm.
#' @param X Input matrix, of dimension nobs*nvars; each row is an observation.
#' If the intercept term needs to be considered in the estimation process, then the first column of \code{X} must be all 1s.
#' @param y Response variable. Quantitative for \code{family="gaussian"}, or
#' \code{family="poisson"} (non-negative counts). For \code{family="binomial"} should be either a factor with two levels.
#' For \code{family="multinomial"} should be a one-hot matrix.
#' For \code{family="cox"} should be an n*2 matrix, one column should be named "time", indicating the survival time;
#' the other column must be named "status", and consists of 0 and 1, 0 indicates that the row of data is censored, 1 is opposite.
#' @param family Either a character string representing one of the built-in families,
#' \code{"gaussian"}, \code{"binomial"}, \code{"poisson"}, \code{"multinomial"} or \code{"cox"}.
#' @param alpha Hyperparameter. In general, alpha can be set to 1, 2 or 3.
#' for \code{family="gaussian"} and \code{family="cox"}, the suggested value for alpha is 2 or 3.
#' for \code{family="binomial"}, \code{family="poisson"} and \code{family="multinomial"}, the suggested value for alpha is 1 or 2.
#' but when the collinearity of the load matrix is serious, the hyperparameters can be selected larger, such as 5.
#' @param itrNum The number of iteration steps. In general, 20 steps are enough.
#' If the condition number of \code{X} is large, it is recommended to greatly increase the
#' number of iteration steps.
#' @param QR_flag It identifies whether to use QR decomposition to speed up the algorithm.
#' Currently only valid for linear models.
#' @param flag It identifies whether to make model selection. The default is \code{TRUE}.
#' @param lamda_0 The initial value of the regularization parameter for ridge regression.
#' The running result of the algorithm is not sensitive to this value.
#' @param fdiag It identifies whether to use diag Approximation to speed up the algorithm.
#'
#' @return Regression coefficients
#' @export GAGA
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
#' ind1 = sample(1:p_size,p_size)
#' ind2 = ind1[1:zeroNum]
#' beta_true = runif(p_size,0,R2)
#' beta_true[ind2] = 0
#' X = R1*matrix(rnorm(sample_size * p_size), ncol = p_size)
#' X[1:sample_size,1]=1
#' y = rpois(sample_size,lambda = as.vector(exp(X%*%beta_true)))
#' y = as.vector(y)
#' # Estimate
#' Eb = GAGA(X,y,alpha = 2,family="poisson")
#' cat("\n err:", norm(Eb-beta_true,type="2")/norm(beta_true,type="2"))
#' cat("\n acc:", cal.w.acc(as.character(Eb!=0),as.character(beta_true!=0)))
#'
#' # Multinomial
#' set.seed(2022)
#' library(mvtnorm)
#' p_size = 20
#' C = 3
#' sample_size = 500
#' test_size = 1000
#' rate = 0.5 #Proportion of value zero in beta
#' R1 = 1
#' R2 = 5
# Set true beta
#' beta_true = matrix(rep(0,p_size*C),c(p_size,C))
#' zeroNum = round(rate*p_size)
#' for(jj in 1:C){
#'   ind1 = sample(1:p_size,p_size)
#'   ind2 = ind1[1:zeroNum]
#'   tmp = runif(p_size,0,R2)
#'   tmp[ind2] = 0
#'   beta_true[,jj] = tmp
#' }
#' cov_mat=matrix(1:p_size*p_size,p_size,p_size)
#' for(i in 1:p_size){for(j in 1:p_size) {if(i!=j){cov_mat[i,j]=0.0}else{cov_mat[i,j]=1}}}
# Generate training samples
#' X = R1*rmvnorm(n=sample_size, mean=rep(0, nrow(cov_mat)), sigma=cov_mat)
#' X[1:sample_size,1]=1
#' z = X%*%beta_true
#' t = exp(z)/(1+rowSums(exp(z)))
#' t = cbind(t,1-rowSums(t))
#' tt = t(apply(t,1,cumsum))
#' tt = cbind(rep(0,sample_size),tt)
#' y = matrix(rep(0,sample_size*(C+1)),c(sample_size,C+1))
#' for(jj in 1:sample_size){
#'   tmp = runif(1,0,1)
#'   for(kk in 1:(C+1)){
#'     if((tmp>tt[jj,kk])&&(tmp<=tt[jj,kk+1])){
#'       y[jj,kk] = 1
#'       break
#'     }
#'   }
#' }
#' # Estimate
#' Eb = GAGA(X, y,alpha=1,family = "multinomial")
# Prediction
# Generate test samples
#' X_t = R1*rmvnorm(n=test_size, mean=rep(0, nrow(cov_mat)), sigma=cov_mat)
#' X_t[1:test_size,1]=1
#' z = X_t%*%beta_true
#' t = exp(z)/(1+rowSums(exp(z)))
#' t = cbind(t,1-rowSums(t))
#' tt = t(apply(t,1,cumsum))
#' tt = cbind(rep(0,test_size),tt)
#' y_t = rep(0,test_size)
#' for(jj in 1:test_size){
#'   tmp = runif(1,0,1)
#'   for(kk in 1:(C+1)){
#'     if((tmp>tt[jj,kk])&&(tmp<=tt[jj,kk+1])){
#'       y_t[jj] = kk
#'       break
#'     }
#'   }
#' }
#' Ey = rep(0,test_size)
#' z = X_t%*%Eb
#' t = exp(z)/(1+rowSums(exp(z)))
#' t = cbind(t,1-rowSums(t))
#' for(jj in 1:test_size){
#'   Ey[jj] = which.max(t[jj,])
#' }
#' cat("\n err:", norm(Eb-beta_true,type="2")/norm(beta_true,type="2"))
#' cat("\n acc:", cal.w.acc(as.character(Eb!=0),as.character(beta_true!=0)))
#' cat("\n pacc:", cal.w.acc(as.character(Ey!=0),as.character(y_t!=0)))
#'
#' # Cox
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
#' ind1 = sample(1:p_size,p_size)
#' ind2 = ind1[1:zeroNum]
#' beta_true = runif(p_size,-R2,R2)
#' beta_true[ind2] = 0
#' # Generate training samples
#' cov_mat=matrix(1:p_size*p_size,p_size,p_size) ##covariance matrix
#' for(i in 1:p_size){for(j in 1:p_size) {if(i!=j){cov_mat[i,j]=0.0}else{cov_mat[i,j]=1}}}
#' X = R1*rmvnorm(n=sample_size, mean=rep(0, nrow(cov_mat)), sigma=cov_mat)
#' z = X%*%beta_true
#' u = runif(sample_size,0,1)
#' t = ((-log(1-u)/(3*exp(z)))*100)^(0.1)
#' cs = rep(0,sample_size)
#' csNum = round(censoringRate*sample_size)
#' ind1 = sample(1:sample_size,sample_size)
#' ind2 = ind1[1:csNum]
#' cs[ind2] = 1
#' t[ind2] = runif(csNum,0,0.8)*t[ind2]
#' y = cbind(t,1 - cs)
#' colnames(y) = c("time", "status")
#' #Estimation
#' Eb = GAGA(X,y,alpha=2,family="cox")
#' cat("\n err:", norm(Eb-beta_true,type="2")/norm(beta_true,type="2"))
#' cat("\n acc:", cal.w.acc(as.character(Eb!=0),as.character(beta_true!=0)))
#'
#'

GAGA <- function(X,y,family=c("gaussian","binomial","poisson","multinomial","cox"),alpha=2,itrNum=20,QR_flag=FALSE,flag=TRUE,lamda_0=0.001,fdiag=TRUE) {

  if(!is.character(family)){
    print("Please check the input of family")
    family = "gaussian"
  }

  family=match.arg(family)
  beta=switch(family,
             "gaussian"=LM_GAGA(X,y,alpha=alpha,itrNum=itrNum,QR_flag=QR_flag,flag=flag,lamda_0=lamda_0,fix_sigma=FALSE,sigm2_0 = 1),
             "poisson"=poisson_GAGA(X,y,alpha=alpha,itrNum=itrNum,flag=flag,lamda_0=lamda_0,fdiag=fdiag),
             "binomial"=logistic_GAGA(X,y,alpha=alpha,itrNum=itrNum,flag=flag,lamda_0=lamda_0,fdiag=fdiag),
             "multinomial"=multinomial_GAGA(X,y,alpha=alpha,itrNum=itrNum,flag=flag,lamda_0=lamda_0,fdiag=fdiag),
             "cox"=cox_GAGA(X,y,alpha=alpha,itrNum=itrNum,flag=flag,lamda_0=lamda_0,fdiag=fdiag)
            )

  return(beta)
}
