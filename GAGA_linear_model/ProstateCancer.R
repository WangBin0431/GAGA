rm(list = ls())
gc()
setwd("F:/投稿论文/JRSSB/GAGA_numericalexperiment/ProstateCancer")
library (lasso2)
data(Prostate)
dim(Prostate)
library(ncvreg)


#产生设计阵和响应向量
data_Pro=Prostate
X=model.matrix(lpsa~.,data_Pro)
#x=model.matrix(lpsa~.,Prostate)[,-1]#不使用截距项
y=Prostate$lpsa
samplesize=length(y)



#设置训练集，设置一个随机种子，保证实验结果可重复
set.seed(1234)
train_size=round(samplesize*0.9)
test_size=samplesize-train_size
expnum = 100##实验次数


test_error_GAGA=NULL
test_error_GAGA_QR=NULL
test_error_ALASSO=NULL
test_error_MCP=NULL
test_error_SCAD=NULL

for(iter in 1:expnum){
train = sample(1:samplesize,train_size)

##Adaptive LASSO
LSE=lm(y[train]~X[train,]-1)## Linear Regression to create the Adaptive Weights Vector
weight=abs(LSE$coefficients)^1# Using gamma = 1
XW=X[train,]%*%diag(weight)##消除权重后的设计???

# fit_ALASSO <- ncvreg(XW, y,family="gaussian", penalty="lasso",nlambda =Nlambd)
Nlambda=100##设置100个lambda
cvfit_ALASSO<- cv.ncvreg(XW, y[train],family="gaussian", penalty="lasso",nlambda =Nlambda, nfolds=10)
fit_ALASSO = cvfit_ALASSO$fit
tmp = fit_ALASSO$beta[,cvfit_ALASSO$min][-1]
tmp = weight*tmp
test_error_ALASSO[iter]=norm(y[-train]-X[-train,]%*%tmp,'f')/test_size


## SCAD 
# fit_SCAD <- ncvreg(X, y,family="gaussian", penalty="SCAD",nlambda =Nlambda)
cvfit_SCAD <- cv.ncvreg(X[train,], y[train],family="gaussian", penalty="SCAD",nlambda =Nlambda, nfolds=10)
fit_SCAD = cvfit_SCAD$fit

tmp = fit_SCAD$beta[,cvfit_SCAD$min][-1]
test_error_SCAD[iter]=norm(y[-train]-X[-train,]%*%tmp,'f')/test_size

## MCP
# fit_MCP <- ncvreg(X, y,family="gaussian", penalty="MCP",nlambda =Nlambda)
cvfit_MCP <- cv.ncvreg(X[train,], y[train],family="gaussian", penalty="MCP",nlambda =Nlambda, nfolds=10)
fit_MCP = cvfit_MCP$fit

tmp = fit_MCP$beta[,cvfit_MCP$min][-1]
test_error_MCP[iter]=norm(y[-train]-X[-train,]%*%tmp,'f')/test_size


##GaGa
source("GaGa.R")
mratio = 2
EW = GaGa(X[train,],y[train],ratio = mratio)
test_error_GAGA[iter]=norm(y[-train]-X[-train,]%*%EW,'f')/test_size

EW2 = GaGa(X[train,],y[train],ratio = mratio,QR_flag = T)
test_error_GAGA_QR[iter]=norm(y[-train]-X[-train,]%*%EW2,'f')/test_size

}

TEST_ERR=c(test_error_GAGA,test_error_GAGA_QR,test_error_ALASSO,test_error_SCAD,test_error_MCP)
Algorithms=factor(c(rep('GAGA',expnum),rep('GAGA_QR',expnum),rep('ALASSO_CV',expnum),rep('SCAD_CV',expnum),rep('MCP_CV',expnum)),
            levels=c('GAGA','GAGA_QR','ALASSO_CV','SCAD_CV','MCP_CV'))

#作图
library(ggplot2)
TEST_ERR_Graph=data.frame(TEST_ERR,Algorithms)
g=ggplot(TEST_ERR_Graph, aes(x=Algorithms, y=TEST_ERR,fill=Algorithms))+ylab("TEST ERR") + geom_boxplot()
g





