#设置工作路径
setwd("F:/投稿论文/JRSSB/GAGA_numericalexperiment/Experiment_lambda")
library(MASS)
#导入lasso计算库
library(ncvreg)
#导入多元正态分布库
library(mvtnorm)
rm(list = ls())
#设置了随机数种子后，每次产生的分布数都是相同的
set.seed(12)





p_size =2;
sample_size=100;##样本量

ratio=0.5



expnum = 1##实验次数

Mean=0##均值
Sd=1##标准差


#初始化设计阵行向量随机产生的协方差矩阵
cov_mat=matrix(1:p_size*p_size,p_size,p_size) ##协方差矩阵
#for(i in 1:p_size){for(j in 1:p_size) {cov_mat[i,j]=0.5^{abs(i-j)}}}
for(i in 1:p_size){for(j in 1:p_size) {if(i!=j){cov_mat[i,j]=0}else{cov_mat[i,j]=1}}}



  
  #随机产生非零信号
  #signal=runif(p_size*ratio,0.5,1)
  signal=1
  #真参数设置
  beta_true=c(signal,rep(0,p_size*(1-ratio)))
  #保存真参数
  write.table(beta_true,"beta_true.txt",row.names = FALSE,col.names = FALSE,quote=FALSE,append=TRUE)



#保存真参数
#write.table(beta_true,"beta_true.txt",row.names = FALSE,col.names = FALSE,quote=FALSE,append=TRUE)


#产生随机设计阵
X=rmvnorm(n=sample_size, mean=rep(0, nrow(cov_mat)), sigma=cov_mat)
#qrresult=qr(X)
#X=sqrt(sample_size)*qr.Q(qrresult)


          
#产生正态随机扰动
raodong=rnorm(sample_size,mean=Mean,sd=Sd)



##产生响应变量
y=X%*%beta_true+raodong
#保存响应变量
#write.table(y,"response_Mod.txt",row.names = FALSE,col.names = FALSE,quote=FALSE,append=TRUE)

#browser()
##GaGa
source("GaGa.R")
mratio = 2
Num2=51
out = GaGa(X,y,ratio = mratio, itrNum=Num2)



EW=out[1];
BBs=out[[2]];


#browser()
BB1=BBs[1:p_size]

Num=Num2-41 #最后一次更新有调整，不看了
lambda_zero_com=matrix(0,nrow = Num,ncol = (p_size*ratio))
lambda_nonzero_com=matrix(0,nrow = Num,ncol = (p_size-p_size*ratio))
beta_mat=matrix(0,nrow = Num,ncol = p_size)

for(index in 1:Num){
  #cat('index=',index)
  nodes=((index-1)*p_size+1):((index-1)*p_size+p_size)
  lambda_nonzero_com[index,]=BBs[nodes[1:(p_size*ratio)]]
  lambda_zero_com[index,]=BBs[nodes[(p_size*ratio+1):p_size]]
}




Iteration=rep(1:Num,p_size*ratio)

lambda_nonzero_component=c(lambda_nonzero_com)
lambda_zero_component=c(lambda_zero_com)




Nonzero=data.frame(Iteration,lambda_nonzero_component)
Zero=data.frame(Iteration,lambda_zero_component)

library(ggplot2)

g=ggplot(Nonzero,aes(x=Iteration,y=lambda_nonzero_component))+ geom_line(col="blue",size=1.5)+geom_point(
  size=3,shape=21, fill="white")+scale_x_continuous(breaks=c(2,4,6,8,10))+ylab('Hyperparameter for the nonzero component')
g2=ggplot(Zero,aes(x=Iteration,y=lambda_zero_component))+ geom_line(col="red",size=1.5)+geom_point(
  size=3, shape=24, fill="white")+scale_x_continuous(breaks=c(2,4,6,8,10))+ylab('Hyperparameter for the zero component')


# g=ggplot(Nonzero,aes(x=Iteration,y=lambda_nonzero_component))+ geom_line(col="blue",size=1.5)+geom_point(
#   size=3,shape=21, fill="white")+scale_x_continuous(breaks=c(2,4,6,8,10))+ylab(expression(paste('Hyperparameter ',
#   lambda[1],' for the nonzero component')))
# g2=ggplot(Zero,aes(x=Iteration,y=lambda_zero_component))+ geom_line(col="red",size=1.5)+geom_point(
#   size=3, shape=24, fill="white")+scale_x_continuous(breaks=c(2,4,6,8,10))+ylab(expression(paste('Hyperparameter ',
#   lambda[2],' for the zero component')))




