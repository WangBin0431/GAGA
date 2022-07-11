#For model 3
setwd("I:/GAGA/Github_project/GAGA/GAGA_linear_model")


library(ncvreg)
library(mvtnorm)
rm(list = ls())
set.seed(1234)

source("LM_GAGA.R")

Nlambda=100
p_size = 8
sample_size=50

sample_size_list = c(30,60,90,120,150)
SL=length(sample_size_list)

expnum = 100

Mean=0
Sd=1
rr = 0.6
mratio = 2

cov_mat = diag(rep(1,p_size))


ERR_ALASSO=matrix(1:expnum*SL,expnum,SL)
ERR_SCAD=matrix(1:expnum*SL,expnum,SL)
ERR_MCP=matrix(1:expnum*SL,expnum,SL)



ACC_ALASSO=matrix(1:expnum*SL,expnum,SL)
ACC_SCAD=matrix(1:expnum*SL,expnum,SL)
ACC_MCP=matrix(1:expnum*SL,expnum,SL)

ERR_GAGA = matrix(1:expnum*SL,expnum,SL)
ACC_GAGA = matrix(1:expnum*SL,expnum,SL)
ERR_GAGA_QR = matrix(1:expnum*SL,expnum,SL)
ACC_GAGA_QR = matrix(1:expnum*SL,expnum,SL)



for(ii in 1:length(sample_size_list)){
  
  cat("iteration = ", ii, "\n")
  
  for(iter in 1:expnum){
    sample_size = sample_size_list[ii]
    
    #set true beta
    zeroNum = round(rr*p_size)
    ind1 = sample(1:p_size,p_size)
    ind2 = ind1[1:zeroNum]
    beta_true = runif(p_size,0,5)
    beta_true[ind2] = 0
    ind3 = ind1[(zeroNum+1):p_size]
    dd = round(length(ind3)/3);
    beta_true[ind3] = runif(length(ind3),0,1)
    pos_true=ind3;
    pos_false=ind2;
  
    X=rmvnorm(n=sample_size, mean=rep(0, nrow(cov_mat)), sigma=cov_mat)
    raodong=rnorm(sample_size,mean=Mean,sd=Sd)
 
    y=X%*%beta_true+raodong
    
    
    ##Adaptive LASSO
    LSE=lm(y~X-1)## Linear Regression to create the Adaptive Weights Vector
    weight=abs(LSE$coefficients)^1# Using gamma = 1
    XW=X%*%diag(weight)
    cvfit_ALASSO<- cv.ncvreg(XW, y,family="gaussian", penalty="lasso",nlambda =Nlambda, nfolds=10)
    fit_ALASSO = cvfit_ALASSO$fit
    
    
    tmp = fit_ALASSO$beta[,cvfit_ALASSO$min][-1]
    tmp = weight*tmp
    ERR_ALASSO[iter,ii]=norm(as.matrix(beta_true-tmp),'f')
    pos1=which(tmp!=0);
    pos2=which(tmp==0);
    ACC_ALASSO[iter,ii]=(length(intersect(pos1,pos_true))+length(intersect(pos2,pos_false)))/length(beta_true)
   
    cvfit_SCAD <- cv.ncvreg(X, y,family="gaussian", penalty="SCAD",nlambda =Nlambda, nfolds=10)
    fit_SCAD = cvfit_SCAD$fit
    
    tmp = fit_SCAD$beta[,cvfit_SCAD$min][-1]
    ERR_SCAD[iter,ii]=norm(as.matrix(beta_true-tmp),'f')
    pos1=which(tmp!=0);
    pos2=which(tmp==0);
    ACC_SCAD[iter,ii]=(length(intersect(pos1,pos_true))+length(intersect(pos2,pos_false)))/length(beta_true)
    
    
    ## MCP
    cvfit_MCP <- cv.ncvreg(X, y,family="gaussian", penalty="MCP",nlambda =Nlambda, nfolds=10)
    fit_MCP = cvfit_MCP$fit
    
    tmp = fit_MCP$beta[,cvfit_MCP$min][-1]
    ERR_MCP[iter,ii]=norm(as.matrix(beta_true-tmp),'f')
    pos1=which(tmp!=0);
    pos2=which(tmp==0);
    ACC_MCP[iter,ii]=(length(intersect(pos1,pos_true))+length(intersect(pos2,pos_false)))/length(beta_true)
    
    ##GaGa
    EW = LM_GAGA(X,y,ratio = mratio)
    
    ERR_GAGA[iter,ii]=norm(as.matrix(beta_true-EW),'f')
    
    pos_GAGA=which(EW!=0);
    pos2_GAGA=which(EW==0);
    ACC_GAGA[iter,ii]=(length(intersect(pos_GAGA,pos_true))+length(intersect(pos2_GAGA,pos_false)))/p_size
    
    EW2 = LM_GAGA(X,y,ratio = mratio,QR_flag = T)
    
    ERR_GAGA_QR[iter,ii]=norm(as.matrix(beta_true-EW2),'f')
    
    pos_GAGA_QR=which(EW2!=0);
    pos2_GAGA_QR=which(EW2==0);
    ACC_GAGA_QR[iter,ii]=(length(intersect(pos_GAGA_QR,pos_true))+length(intersect(pos2_GAGA_QR,pos_false)))/p_size
    
  }#for(iter in 1:expnum)
  
}#for(ii in 1:length(sample_size_list))




Mean_ERR_ALASSO=colMeans(ERR_ALASSO)
Mean_ERR_SCAD=colMeans(ERR_SCAD)
Mean_ERR_MCP=colMeans(ERR_MCP)

Mean_ACC_ALASSO=colMeans(ACC_ALASSO)
Mean_ACC_SCAD=colMeans(ACC_SCAD)
Mean_ACC_MCP=colMeans(ACC_MCP)

Mean_ERR_GAGA = colMeans(ERR_GAGA)
Mean_ERR_GAGA_QR = colMeans(ERR_GAGA_QR)
Mean_ACC_GAGA = colMeans(ACC_GAGA)
Mean_ACC_GAGA_QR = colMeans(ACC_GAGA_QR)

Std_ACC_MCP=sqrt(apply(ACC_MCP, 2, var))
Std_ACC_SCAD=sqrt(apply(ACC_SCAD, 2, var))
Std_ACC_ALASSO=sqrt(apply(ACC_ALASSO, 2, var))
Std_ERR_MCP=sqrt(apply(ERR_MCP, 2, var))
Std_ERR_SCAD=sqrt(apply(ERR_SCAD, 2, var))
Std_ERR_ALASSO=sqrt(apply(ERR_ALASSO, 2, var))

Std_ERR_GAGA = sqrt(apply(ERR_GAGA, 2, var))
Std_ERR_GAGA_QR = sqrt(apply(ERR_GAGA_QR, 2, var))
Std_ACC_GAGA = sqrt(apply(ACC_GAGA, 2, var))
Std_ACC_GAGA_QR = sqrt(apply(ACC_GAGA_QR, 2, var))



ERR=c(Mean_ERR_ALASSO,Mean_ERR_SCAD,Mean_ERR_MCP,Mean_ERR_GAGA,Mean_ERR_GAGA_QR)
ACC=c(Mean_ACC_ALASSO,Mean_ACC_SCAD,Mean_ACC_MCP,Mean_ACC_GAGA,Mean_ACC_GAGA_QR)

xaes = c(sample_size_list,sample_size_list,sample_size_list,sample_size_list,sample_size_list)

Algorithms=factor(c(rep('ALASSO_CV',SL),rep('SCAD_CV',SL),rep('MCP_CV',SL),rep('GAGA',SL),rep('GAGA_QR',SL)),
            levels=c('GAGA','GAGA_QR','ALASSO_CV','SCAD_CV','MCP_CV'))


ERR_SL=data.frame(ERR,xaes,Algorithms)
ACC_SL=data.frame(ACC,xaes,Algorithms)


library(ggplot2)
g1=ggplot(ERR_SL,aes(xaes,ERR,shape=Algorithms,color=Algorithms))+ geom_point(size=3)+xlab("Sample Size") + geom_line()
g1

g2=ggplot(ACC_SL,aes(xaes,ACC,shape=Algorithms,color=Algorithms))+ geom_point(size=3)+xlab("Sample Size") + geom_line()
g2

# 
# 
# boxplot(ERR_GAGA,ERR_GAGA_QR,ERR_ALASSO_CV,ERR_SCAD_CV,ERR_MCP_CV,names=c("GAGA","GAGA_QR","ALASSO_CV","SCAD_CV","MCP_CV"),ylab = "ERR")
# boxplot(ACC_GAGA,ACC_GAGA_QR,ACC_ALASSO_CV,ACC_SCAD_CV,ACC_MCP_CV,names=c("GAGA","GAGA_QR","ALASSO_CV","SCAD_CV","MCP_CV"),ylab = "ACC")
