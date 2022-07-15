

#' Title
#'
#' @param x
#' @param y
#' @param family
#' @param alpha
#' @param itrNum
#' @param QR_flag
#' @param flag
#' @param lamda_0
#' @param fdiag
#'
#' @return
#' @export
#'
#' @examples
GAGA <- function(X,y,family=c("gaussian","binomial","poisson","multinomial","cox"),alpha=2,itrNum=30,QR_flag=FALSE,flag=TRUE,lamda_0=0.001,fdiag=TRUE) {

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
