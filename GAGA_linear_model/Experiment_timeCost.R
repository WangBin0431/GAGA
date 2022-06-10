
#设置工作路径
setwd("F:/投稿论文/JRSSB/GAGA_numericalexperiment/Experiment_timeCost")
#导入lasso计算库
library(ncvreg)
#导入多元正态分布库
library(mvtnorm)
rm(list = ls())
#设置了随机数种子后，每次产生的分布数都是相同的
set.seed(1234)


Nlambda=100##设置100个lambda
p_size_list = c(500,1000,2000)
sample_size=4000##样本量
expnum = 10##实验次数

Mean=0##均值
Sd=1##标准差

rr = 0.5






TC_GAGA = NULL
TC_GAGA_QR = NULL
TC_ALASSO_CV = NULL
TC_SCAD_CV = NULL
TC_MCP_CV = NULL

Mean_TC_GAGA=NULL
Mean_TC_GAGA_QR=NULL
Mean_TC_ALASSO_CV=NULL
Mean_TC_SCAD_CV=NULL
Mean_TC_MCP_CV=NULL

for(ii in 1:length(p_size_list)){
  
  cat("ip_size_list is ", ii, "\n")

  #初始化设计阵行向量随机产生的协方差矩阵
  p_size=p_size_list[ii];
  cov_mat=matrix(1:p_size*p_size,p_size,p_size) ##协方差矩阵
  for(i in 1:p_size){for(j in 1:p_size) {if(i!=j){cov_mat[i,j]=0.5}else{cov_mat[i,j]=1}}}

for(iter in 1:expnum){#每个样本量下，实验重复expnum次
  cat("expnum is:",iter,"\n");
  
  #真参数设置
  zeroNum = round(rr*p_size)
  ind1 = sample(1:p_size,p_size)
  ind2 = ind1[1:zeroNum]
  beta_true = runif(p_size,0,5)
  
  beta_true[ind2] = 0
  
  ind3 = ind1[(zeroNum+1):p_size]
  
  #dd = round(length(ind3)/3);
  # be = 1
  # beta_true[ind3[be:(be+dd-1)]] = runif(dd,0,1)
  # be = be+dd;
  # beta_true[ind3[be:(be+dd-1)]] = runif(dd,0,1)
  # be = be +dd;
  # beta_true[ind3[be:length(ind3)]] = runif(length(ind3)-be+1,20,50)
  #beta_true[ind3] = runif(length(ind3),0,1)
  
  pos_true=ind3;
  pos_false=ind2;
  
  
  
  
  #产生随机设计阵
  X=rmvnorm(n=sample_size, mean=rep(0, nrow(cov_mat)), sigma=cov_mat)
  
  
  
  
  #产生正态随机扰动
  raodong=rnorm(sample_size,mean=Mean,sd=Sd)
  
  
  
  ##产生响应变量
  y=X%*%beta_true+raodong
  #保存响应变量
  #write.table(y,"response_Mod.txt",row.names = FALSE,col.names = FALSE,quote=FALSE,append=TRUE)
  
  
  # #耗时测试
  # timeb = proc.time();
  # beta = solve(t(X)%*%X,tol = 0)%*%t(X)%*%y
  # timee = proc.time();
  # cat("Time of Bin:",timee[3]-timeb[3],"\n");
  # 
  # 
  # timeb = proc.time();
  # beta = solve(t(X)%*%X ,t(X)%*%y)
  # timee = proc.time();
  # cat("Time of Fei:",timee[3]-timeb[3],"\n");
  
  
  ##GaGa
  source("GaGa.R")
  mratio = 2
  
  timeb = proc.time();
  EW = GaGa(X,y,ratio = mratio)
  timee = proc.time();
  TC_GAGA[iter]=timee[3]-timeb[3];
  
  
  timeb = proc.time();
  EW2 = GaGa(X,y,ratio = mratio,QR_flag = T)
  timee = proc.time();
  TC_GAGA_QR[iter]=timee[3]-timeb[3];
  

  
  ##Adaptive LASSO
  timeb = proc.time();
  LSE=lm(y~X-1)## Linear Regression to create the Adaptive Weights Vector
  weight=abs(LSE$coefficients)^1# Using gamma = 1
  XW=X%*%diag(weight)##消除权重后的设计阵
  cvfit_ALASSO<- cv.ncvreg(XW, y,family="gaussian", penalty="lasso",nlambda =Nlambda, nfolds=10)
  fit_ALASSO = cvfit_ALASSO$fit
  tmp = fit_ALASSO$beta[,cvfit_ALASSO$min][-1]
  tmp = weight*tmp
  timee = proc.time();
  TC_ALASSO_CV[iter]=timee[3]-timeb[3];
  

  
  
  
  ## SCAD 
  timeb = proc.time();
  cvfit_SCAD <- cv.ncvreg(X, y,family="gaussian", penalty="SCAD",nlambda =Nlambda, nfolds=10)
  timee = proc.time();
  TC_SCAD_CV[iter]=timee[3]-timeb[3];
  
  
  ## MCP
  timeb = proc.time();
  cvfit_MCP <- cv.ncvreg(X, y,family="gaussian", penalty="MCP",nlambda =Nlambda, nfolds=10)
  timee = proc.time();
  TC_MCP_CV[iter]=timee[3]-timeb[3];
  

  
  
  
}#for(iter in 1:expnum)
  Mean_TC_GAGA[ii]=mean(TC_GAGA)
  Mean_TC_GAGA_QR[ii]=mean(TC_GAGA_QR)
  Mean_TC_ALASSO_CV[ii]=mean(TC_ALASSO_CV)
  Mean_TC_SCAD_CV[ii]=mean(TC_SCAD_CV)
  Mean_TC_MCP_CV[ii]=mean(TC_MCP_CV)
  
  
}

TimeCost=data.frame(Mean_TC_GAGA,Mean_TC_GAGA_QR,Mean_TC_ALASSO_CV,Mean_TC_SCAD_CV,Mean_TC_MCP_CV)
write.table(TimeCost,"timecost.txt",row.names = FALSE)

