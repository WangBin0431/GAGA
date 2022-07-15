# multinomial
set.seed(2022)
library(mvtnorm)
cat("\n")
cat("Test multinomial GAGA\n")

p_size = 20
C = 3
sample_size = 500
test_size = 1000

rate = 0.5 #Proportion of value zero in beta
Num = 10 # Total number of experiments
R1 = 1
R2 = 5

#Set true beta
beta_true = matrix(rep(0,p_size*C),c(p_size,C))
zeroNum = round(rate*p_size)
for(jj in 1:C){
  ind1 = sample(1:p_size,p_size)
  ind2 = ind1[1:zeroNum]
  tmp = runif(p_size,0,R2)
  tmp[ind2] = 0
  beta_true[,jj] = tmp
}

cov_mat=matrix(1:p_size*p_size,p_size,p_size)
for(i in 1:p_size){for(j in 1:p_size) {if(i!=j){cov_mat[i,j]=0.0}else{cov_mat[i,j]=1}}}

#Generate training samples
X = R1*rmvnorm(n=sample_size, mean=rep(0, nrow(cov_mat)), sigma=cov_mat)

X[1:sample_size,1]=1
z = X%*%beta_true
t = exp(z)/(1+rowSums(exp(z)))
t = cbind(t,1-rowSums(t))
tt = t(apply(t,1,cumsum))
tt = cbind(rep(0,sample_size),tt)

y = matrix(rep(0,sample_size*(C+1)),c(sample_size,C+1))
for(jj in 1:sample_size){
  tmp = runif(1,0,1)
  for(kk in 1:(C+1)){
    if((tmp>tt[jj,kk])&&(tmp<=tt[jj,kk+1])){
      y[jj,kk] = 1
      break
    }
  }
}

Eb = GAGA(X, y,alpha=1,family = "multinomial")

#Prediction#######################################################################################################
#Generate test samples
X_t = R1*rmvnorm(n=test_size, mean=rep(0, nrow(cov_mat)), sigma=cov_mat)
X_t[1:test_size,1]=1
z = X_t%*%beta_true
t = exp(z)/(1+rowSums(exp(z)))
t = cbind(t,1-rowSums(t))
tt = t(apply(t,1,cumsum))
tt = cbind(rep(0,test_size),tt)

y_t = rep(0,test_size)
#y_tt = rep(0,test_size)
for(jj in 1:test_size){
  tmp = runif(1,0,1)
  for(kk in 1:(C+1)){
    if((tmp>tt[jj,kk])&&(tmp<=tt[jj,kk+1])){
      y_t[jj] = kk
      break
    }
  }
  #y_tt[jj] = which.max(t[jj,])
}

Ey = rep(0,test_size)
z = X_t%*%Eb
t = exp(z)/(1+rowSums(exp(z)))
t = cbind(t,1-rowSums(t))
for(jj in 1:test_size){
  Ey[jj] = which.max(t[jj,])
}


cat("\n--------------------")
cat("\n err:", norm(Eb-beta_true,type="2")/norm(beta_true,type="2"))
cat("\n acc:", cal.w.acc(as.character(Eb!=0),as.character(beta_true!=0)))
cat("\n pacc:", cal.w.acc(as.character(Ey!=0),as.character(y_t!=0)))
cat("\n")


