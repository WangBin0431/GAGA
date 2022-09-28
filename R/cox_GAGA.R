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

  vnames=colnames(X)
  if(is.null(vnames))vnames=paste("V",seq(ncol(X)),sep="")
  fit = list()
  class(fit) = c("GAGA","cox")

  y = t[,"time"]
  cens = 1-t[,"status"]

  tmpfit = cpp_COX_gaga(X, as.matrix(y), as.matrix(cens), alpha, itrNum, thresh, flag, lamda_0, fdiag)

  fit$beta = as.vector(tmpfit$beta)
  names(fit$beta) = vnames
  fit$alpha = alpha
  fit$itrNum = tmpfit$itrNum
  fit$fdiag = fdiag
  return(fit)
}

