
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
#' @param thresh Convergence threshold for beta Change, if \code{max(abs(beta-beta_old))<threshold}, return.
#' @param QR_flag It identifies whether to use QR decomposition to speed up the algorithm.
#' Currently only valid for linear models.
#' @param flag It identifies whether to make model selection. The default is \code{TRUE}.
#' @param lamda_0 The initial value of the regularization parameter for ridge regression.
#' The running result of the algorithm is not sensitive to this value.
#' @param fix_sigma It identifies whether to update the variance estimate of the Gaussian noise or not.
#' \code{fix_sigma=TRUE} uses the initial variance as the variance estimate in each loop.
#' \code{fix_sigma=FALSE} updates the variance estimate in each loop.
#' @param sigm2_0 The initial variance of the Gaussian noise.
#' @param fdiag It identifies whether to use diag Approximation to speed up the algorithm.
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
#' ind = sample(1:p_size,zeroNum)
#' beta_true = runif(p_size,0,R2)
#' beta_true[ind] = 0
#' X = R1*matrix(rnorm(sample_size * p_size), ncol = p_size)
#' y=X%*%beta_true + rnorm(sample_size,mean=0,sd=2)
#' # Estimation
#' fit = GAGA(X,y,alpha = 3,family="gaussian")
#' Eb = fit$beta
#' #Create testing data
#' X_t = R1*matrix(rnorm(sample_size * p_size), ncol = p_size)
#' y_t=X_t%*%beta_true + rnorm(sample_size,mean=0,sd=2)
#' #Prediction
#' Ey = predict.GAGA(fit,newx=X_t)
#'
#' cat("\n err:", norm(Eb-beta_true,type="2")/norm(beta_true,type="2"))
#' cat("\n acc:", cal.w.acc(as.character(Eb!=0),as.character(beta_true!=0)))
#' cat("\n perr:", norm(Ey-y_t,type="2")/sqrt(sample_size))

LM_GAGA = function(X,y,alpha=3,itrNum=50,thresh=1.e-3,QR_flag=FALSE,flag=TRUE,lamda_0=0.001,fix_sigma=FALSE,sigm2_0 = 1,fdiag = FALSE){

  exitflag = FALSE
  vnames=colnames(X)
  if(is.null(vnames))vnames=paste("V",seq(ncol(X)),sep="")
  fit = list()
  class(fit) = c("GAGA","gaussian")

  tmpfit = rcpp_lm_gaga(X,as.matrix(y),alpha,itrNum,thresh,QR_flag,flag,lamda_0,fix_sigma,sigm2_0,fdiag)

  fit$beta = as.vector(tmpfit$beta)
  names(fit$beta) = vnames
  fit$alpha = alpha
  fit$itrNum = tmpfit$itrNum
  fit$fdiag = fdiag
  return(fit)

}

