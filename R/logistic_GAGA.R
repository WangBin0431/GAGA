#' Fit a logistic model via the Global Adaptive Generative Adjustment algorithm
#'
#' @param X Input matrix, of dimension nobs*nvars; each row is an observation.
#' If the intercept term needs to be considered in the estimation process, then the first column of \code{X} must be all 1s.
#' @param y should be either a factor with two levels.
#' @param alpha Hyperparameter. The suggested value for alpha is 1 or 2.
#' When the collinearity of the load matrix is serious, the hyperparameters can be selected larger, such as 5.
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
#' @export logistic_GAGA
#'
#' @examples
#' # Binomial
#' set.seed(2022)
#' library(mvtnorm)
#' cat("\n")
#' cat("Test binomial GAGA\n")
#' p_size = 30
#' sample_size=600
#' test_size=1000
#' R1 = 1
#' R2 = 3
#' rate = 0.5 #Proportion of value zero in beta
#' # Set true beta
#' zeroNum = round(rate*p_size)
#' ind = sample(1:p_size,zeroNum)#'
#' beta_true = runif(p_size,R2*0.2,R2)
#' beta_true[ind] = 0
#' # Generate training samples
#' X = R1*matrix(rnorm(sample_size * p_size), ncol = p_size)
#' X[1:sample_size,1]=1
#' t = 1/(1+exp(-X%*%beta_true))
#' tmp = runif(sample_size,0,1)
#' y = rep(0,sample_size)
#' y[t>tmp] = 1
#' # Estimation
#' fit = GAGA(X,y,family = "binomial", alpha = 1)
#' Eb = fit$beta
#' # Generate testing samples
#' X_t = R1*matrix(rnorm(test_size * p_size), ncol = p_size)
#' X_t[1:test_size,1]=1
#' t = 1/(1+exp(-X_t%*%beta_true))
#' tmp = runif(test_size,0,1)
#' y_t = rep(0,test_size)
#' y_t[t>tmp] = 1
#' #Prediction
#' Ey = predict.GAGA(fit,newx = X_t)
#' cat("\n--------------------")
#' cat("\n err:", norm(Eb-beta_true,type="2")/norm(beta_true,type="2"))
#' cat("\n acc:", cal.w.acc(as.character(Eb!=0),as.character(beta_true!=0)))
#' cat("\n pacc:", cal.w.acc(as.character(Ey),as.character(y_t)))
#'
logistic_GAGA = function(X,y,alpha=1,itrNum=30,thresh=1.e-3,flag=TRUE,lamda_0=0.001,fdiag=TRUE){


  vnames=colnames(X)
  if(is.null(vnames))vnames=paste("V",seq(ncol(X)),sep="")

  y=as.factor(y)
  ntab=table(y)
  minclass=min(ntab)
  if(minclass<=1)stop("one multinomial or binomial class has 1 or 0 observations; not allowed")
  if(minclass<8)warning("one multinomial or binomial class has fewer than 8  observations; dangerous ground")
  classnames=names(ntab)
  nc=as.integer(length(ntab))
  y=diag(nc)[as.numeric(y),]
  y = y[,2]


  fit = list()
  fit$classnames = classnames
  class(fit) = c("GAGA","binomial")

  tmpfit = cpp_logistic_gaga(X, as.matrix(y), alpha, itrNum, thresh, flag, lamda_0, fdiag)

  fit$beta = as.vector(tmpfit$beta)
  names(fit$beta) = vnames
  fit$alpha = alpha
  fit$itrNum = tmpfit$itrNum
  fit$fdiag = fdiag
  return(fit)
}
