#设置工作路径
setwd("F:/投稿论文/JRSSB/GAGA_numericalexperiment/Experiment3")

#导入lasso计算???
library(ncvreg)
#导入多元正态分布库
library(mvtnorm)
rm(list = ls())



#设置了随机数种子后，每次产生的分布数都是相同???
set.seed(1234)


Nlambda=100##设置100个lambda
p_size = 8
sample_size=50##样本???

sample_size_list = c(30,60,90,120,150)
SL=length(sample_size_list)

expnum = 100##实验次数

Mean=0##??????
Sd=1##标准???
rr = 0.6
mratio = 2


#初始化设计阵行向量随机产生的协方差矩???
# cov_mat=matrix(1:p_size*p_size,p_size,p_size) ##协方差矩???
# for(i in 1:p_size){for(j in 1:p_size) {cov_mat[i,j]=0.5^{abs(i-j)}}}
cov_mat = diag(rep(1,p_size))

#初始化真参数与估计的误差
ERR_ALASSO=matrix(1:expnum*SL,expnum,SL)
ERR_SCAD=matrix(1:expnum*SL,expnum,SL)
ERR_MCP=matrix(1:expnum*SL,expnum,SL)


#初始化真参数与估计的ACC
ACC_ALASSO=matrix(1:expnum*SL,expnum,SL)
ACC_SCAD=matrix(1:expnum*SL,expnum,SL)
ACC_MCP=matrix(1:expnum*SL,expnum,SL)

ERR_GAGA = matrix(1:expnum*SL,expnum,SL)
ACC_GAGA = matrix(1:expnum*SL,expnum,SL)
ERR_GAGA_QR = matrix(1:expnum*SL,expnum,SL)
ACC_GAGA_QR = matrix(1:expnum*SL,expnum,SL)



for(ii in 1:length(sample_size_list)){
  
  cat("iteration = ", ii, "\n")
  
  for(iter in 1:expnum){#每个样本量下，实验重复expnum???
    sample_size = sample_size_list[ii]
    
    #随机产生非零信号
    #signal=c(runif(1,0,1),runif(1,0,1),runif(1,0,1))
    
    #真参数设???
    #beta_true=rep(0,p_size)
    #pos_true=c(1,2,5);#记录非零信号位置
    #pos_false=c(3,4,6,7,8);#记录零信号位???
    #beta_true[1] = signal[1]
    #beta_true[2] = signal[2]
    #beta_true[5] = signal[3]
    
    # if(iter == 1){
    zeroNum = round(rr*p_size)
    ind1 = sample(1:p_size,p_size)
    ind2 = ind1[1:zeroNum]
    beta_true = runif(p_size,0,5)
    
    beta_true[ind2] = 0
    
    ind3 = ind1[(zeroNum+1):p_size]
    
    dd = round(length(ind3)/3);
    
    # beta_true[ind3[1:dd]] = runif(dd,0,0.5)
    # beta_true[ind3[(dd+1):(dd+dd)]] = runif(dd,5,10)
    # beta_true[ind3[(dd+dd+1):length(ind3)]] = runif(length(ind3)-dd-dd,50,100)
    
    beta_true[ind3] = runif(length(ind3),0,1)
    pos_true=ind3;
    pos_false=ind2;
    # }  
    
    
    
    #保存真参???
    #write.table(beta_true,"beta_true.txt",row.names = FALSE,col.names = FALSE,quote=FALSE,append=TRUE)
    
    
    #产生随机设计???
    X=rmvnorm(n=sample_size, mean=rep(0, nrow(cov_mat)), sigma=cov_mat)
    
    #保存设计???
    #write.table(X,"designedmat_Mod.txt",row.names = FALSE,col.names = FALSE,quote=FALSE,append=TRUE)
    
    
    #产生正态随机扰???
    raodong=rnorm(sample_size,mean=Mean,sd=Sd)
    
    
    
    ##产生响应变量
    y=X%*%beta_true+raodong
    #保存响应变量
    #write.table(y,"response_Mod.txt",row.names = FALSE,col.names = FALSE,quote=FALSE,append=TRUE)
    
    ##Adaptive LASSO
    LSE=lm(y~X-1)## Linear Regression to create the Adaptive Weights Vector
    weight=abs(LSE$coefficients)^1# Using gamma = 1
    XW=X%*%diag(weight)##消除权重后的设计???
    # fit_ALASSO <- ncvreg(XW, y,family="gaussian", penalty="lasso",nlambda =Nlambd)
    cvfit_ALASSO<- cv.ncvreg(XW, y,family="gaussian", penalty="lasso",nlambda =Nlambda, nfolds=10)
    fit_ALASSO = cvfit_ALASSO$fit
    #记录各个罚值下，估计的参数和真是参数的二范???
    
    tmp = fit_ALASSO$beta[,cvfit_ALASSO$min][-1]
    tmp = weight*tmp
    ERR_ALASSO[iter,ii]=norm(as.matrix(beta_true-tmp),'f')
    pos1=which(tmp!=0);
    pos2=which(tmp==0);
    ACC_ALASSO[iter,ii]=(length(intersect(pos1,pos_true))+length(intersect(pos2,pos_false)))/length(beta_true)
    
    
    ## SCAD 
    # fit_SCAD <- ncvreg(X, y,family="gaussian", penalty="SCAD",nlambda =Nlambda)
    cvfit_SCAD <- cv.ncvreg(X, y,family="gaussian", penalty="SCAD",nlambda =Nlambda, nfolds=10)
    fit_SCAD = cvfit_SCAD$fit
    
    tmp = fit_SCAD$beta[,cvfit_SCAD$min][-1]
    ERR_SCAD[iter,ii]=norm(as.matrix(beta_true-tmp),'f')
    pos1=which(tmp!=0);
    pos2=which(tmp==0);
    ACC_SCAD[iter,ii]=(length(intersect(pos1,pos_true))+length(intersect(pos2,pos_false)))/length(beta_true)
    
    
    ## MCP
    # fit_MCP <- ncvreg(X, y,family="gaussian", penalty="MCP",nlambda =Nlambda)
    cvfit_MCP <- cv.ncvreg(X, y,family="gaussian", penalty="MCP",nlambda =Nlambda, nfolds=10)
    fit_MCP = cvfit_MCP$fit
    
    tmp = fit_MCP$beta[,cvfit_MCP$min][-1]
    ERR_MCP[iter,ii]=norm(as.matrix(beta_true-tmp),'f')
    pos1=which(tmp!=0);
    pos2=which(tmp==0);
    ACC_MCP[iter,ii]=(length(intersect(pos1,pos_true))+length(intersect(pos2,pos_false)))/length(beta_true)
    
    ##GaGa
    source("GaGa.R")
    EW = GaGa(X,y,ratio = mratio)
    
    ERR_GAGA[iter,ii]=norm(as.matrix(beta_true-EW),'f')
    
    pos_GAGA=which(EW!=0);
    pos2_GAGA=which(EW==0);
    ACC_GAGA[iter,ii]=(length(intersect(pos_GAGA,pos_true))+length(intersect(pos2_GAGA,pos_false)))/p_size
    
    EW2 = GaGa(X,y,ratio = mratio,QR_flag = T)
    
    ERR_GAGA_QR[iter,ii]=norm(as.matrix(beta_true-EW2),'f')
    
    pos_GAGA_QR=which(EW2!=0);
    pos2_GAGA_QR=which(EW2==0);
    ACC_GAGA_QR[iter,ii]=(length(intersect(pos_GAGA_QR,pos_true))+length(intersect(pos2_GAGA_QR,pos_false)))/p_size
    
    
  }#for(iter in 1:expnum)
  
}#for(ii in 1:length(sample_size_list))






#每列??????
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



#每列方差
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

#构造数据框


ERR=c(Mean_ERR_ALASSO,Mean_ERR_SCAD,Mean_ERR_MCP,Mean_ERR_GAGA,Mean_ERR_GAGA_QR)
ACC=c(Mean_ACC_ALASSO,Mean_ACC_SCAD,Mean_ACC_MCP,Mean_ACC_GAGA,Mean_ACC_GAGA_QR)

xaes = c(sample_size_list,sample_size_list,sample_size_list,sample_size_list,sample_size_list)

Algorithms=factor(c(rep('ALASSO_CV',SL),rep('SCAD_CV',SL),rep('MCP_CV',SL),rep('GAGA',SL),rep('GAGA_QR',SL)),
            levels=c('GAGA','GAGA_QR','ALASSO_CV','SCAD_CV','MCP_CV'))

# ERR=c(Mean_ERR_SCAD,Mean_ERR_MCP,Mean_ERR_GAGA,Mean_ERR_GAGA_QR)
# ACC=c(Mean_ACC_SCAD,Mean_ACC_MCP,Mean_ACC_GAGA,Mean_ACC_GAGA_QR)
# type=factor(c(rep('SCAD',Nlambda),rep('MCP',Nlambda),'GAGA','GAGA_QR'),
#             levels=c('GAGA','GAGA_QR','ALASSO','SCAD','MCP'))


ERR_SL=data.frame(ERR,xaes,Algorithms)
ACC_SL=data.frame(ACC,xaes,Algorithms)

#作图
library(ggplot2)
g1=ggplot(ERR_SL,aes(xaes,ERR,shape=Algorithms,color=Algorithms))+ geom_point(size=3)+xlab("Sample Size") + geom_line()
g1

g2=ggplot(ACC_SL,aes(xaes,ACC,shape=Algorithms,color=Algorithms))+ geom_point(size=3)+xlab("Sample Size") + geom_line()
g2

# 
# 
# boxplot(ERR_GAGA,ERR_GAGA_QR,ERR_ALASSO_CV,ERR_SCAD_CV,ERR_MCP_CV,names=c("GAGA","GAGA_QR","ALASSO_CV","SCAD_CV","MCP_CV"),ylab = "ERR")
# boxplot(ACC_GAGA,ACC_GAGA_QR,ACC_ALASSO_CV,ACC_SCAD_CV,ACC_MCP_CV,names=c("GAGA","GAGA_QR","ALASSO_CV","SCAD_CV","MCP_CV"),ylab = "ACC")
