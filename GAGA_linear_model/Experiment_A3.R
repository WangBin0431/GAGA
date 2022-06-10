#设置工作路径
setwd("F:/投稿论文/JRSSB/GAGA_numericalexperiment/Experiment_A3")
library(MASS)
#导入lasso计算库
library(ncvreg)
#导入多元正态分布库
library(mvtnorm)
rm(list = ls())
#设置了随机数种子后，每次产生的分布数都是相同的
set.seed(1234)


Nlambda=100##设置100个lambda

cishu=10
minratio_nonzero_components=NULL
maxratio_zero_components=NULL
for(index2 in 1:cishu){

p_size =20*index2;
sample_size=100*index2;##样本量

ratio=0.5



expnum = 1##实验次数

Mean=0##均值
Sd=1##标准差


#初始化设计阵行向量随机产生的协方差矩阵
cov_mat=matrix(1:p_size*p_size,p_size,p_size) ##协方差矩阵
#for(i in 1:p_size){for(j in 1:p_size) {cov_mat[i,j]=0.5^{abs(i-j)}}}
for(i in 1:p_size){for(j in 1:p_size) {if(i!=j){cov_mat[i,j]=0.5}else{cov_mat[i,j]=1}}}



  
  #随机产生非零信号
  signal=runif(p_size*ratio,0.1,1)
  #真参数设置
  beta_true=c(signal,rep(0,p_size*(1-ratio)))
  #保存真参数
  write.table(beta_true,"beta_true.txt",row.names = FALSE,col.names = FALSE,quote=FALSE,append=TRUE)



#保存真参数
#write.table(beta_true,"beta_true.txt",row.names = FALSE,col.names = FALSE,quote=FALSE,append=TRUE)


#产生随机设计阵
X=rmvnorm(n=sample_size, mean=rep(0, nrow(cov_mat)), sigma=cov_mat)



          
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
beta1=solve(t(X)%*%X+diag(BB1),t(X)%*%y)
a112=1/diag(solve(t(X)%*%X+diag(BB1)))
fenmu=a112*abs(beta1)#信号的分母

#browser()


Num=Num2-1 #最后一次更新有调整，不看了
fenzi_mat=matrix(0,nrow = Num,ncol = p_size)
beta_mat=matrix(0,nrow = Num,ncol = p_size)

for(index in 1:Num){
  #cat('index=',index)
  nodes=((index-1)*p_size+1):((index-1)*p_size+p_size)
  BB_index=BBs[nodes]
  beta_index=solve(t(X)%*%X+diag(BB_index),t(X)%*%y,tol = 0)
  
  #browser()
  
  #beta_index=solve(t(X)%*%X+diag(BB_index))%*%t(X)%*%y
  b112_index=1/diag(solve(t(X)%*%X+diag(BB_index),tol = 0))
  
  #browser()
  
  fenzi_mat[index,]=b112_index*abs(beta_index)#信号的分子
  
  
}

#browser()

signal_ratio=matrix(0,nrow=Num,ncol=p_size*ratio)
for(index in 1:Num){
  signal_ratio[index,]=fenzi_mat[index,1:(p_size*ratio)]/fenmu[1:(p_size*ratio)]#非零信号比值
}
zerosignal_ratio=matrix(0,nrow=Num,ncol=p_size*(1-ratio))
for(index in 1:Num){
  zerosignal_ratio[index,]=fenzi_mat[index,((p_size*ratio)+1):p_size]/max(fenmu[(p_size*ratio+1):p_size])#零信号比值
}

signal_minitr=rep(0,Num)
zerosignal_maxitr=rep(0,Num)
for(index in 1:Num){
  signal_minitr[index]=min(signal_ratio[index,])
  zerosignal_maxitr[index]=max(zerosignal_ratio[index,])
}

#browser()

minratio_nonzero_components=c(minratio_nonzero_components,signal_minitr)
maxratio_zero_components=c(maxratio_zero_components,zerosignal_maxitr)
}




Iter=rep(1:Num,cishu)
size=factor(c(rep('100*20',Num),rep('200*40',Num),rep('300*60',Num),rep('400*80',Num),rep('500*100',Num),
              rep('600*120',Num),rep('700*140',Num),rep('800*160',Num),rep('900*180',Num),rep('1000*200',Num)),
            levels=c('100*20','200*40','300*60','400*80','500*100',
                     '600*120','700*140','800*160','900*180','1000*200'))
minmaxratio=data.frame(Iter,minratio_nonzero_components,maxratio_zero_components,size)
library(ggplot2)
g=ggplot(minmaxratio,aes(x=Iter,y=minratio_nonzero_components,fill=size,
                         color=size))+ geom_line()+geom_point(size=3
                        )+xlab('Iteration')+ylab('Minimal ratios for nonzero components')+ylim(0.8,1.3)
g2=ggplot(minmaxratio,aes(x=Iter,y=maxratio_zero_components,fill=size,
                          color=size))+ geom_line()+geom_point(size=3
                       )+xlab('Iteration')+ylab('Maximal ratios for zero components')
