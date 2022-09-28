
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
#' @param thresh Convergence threshold for beta Change, if \code{max(abs(beta-beta_old))<threshold}, return.
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
poisson_GAGA = function(X,y,alpha=1,itrNum=30,thresh=1.e-3,flag=TRUE,lamda_0=0.5,fdiag=TRUE){

  vnames=colnames(X)
  if(is.null(vnames))vnames=paste("V",seq(ncol(X)),sep="")
  fit = list()
  class(fit) = c("GAGA","poisson")

  tmpfit = cpp_poisson_gaga(X, as.matrix(y), alpha, itrNum, thresh, flag, lamda_0, fdiag)

  fit$beta = as.vector(tmpfit$beta)
  names(fit$beta) = vnames
  fit$alpha = alpha
  fit$itrNum = tmpfit$itrNum
  fit$fdiag = fdiag
  return(fit)
}


