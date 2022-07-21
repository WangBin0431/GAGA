# binomial
set.seed(2022)
library(mvtnorm)
cat("\n")
cat("Test binomial GAGA\n")
p_size = 30
sample_size=600
test_size=1000
R1 = 1
R2 = 3
rate = 0.5 #Proportion of value zero in beta

#Set true beta
zeroNum = round(rate*p_size)
ind1 = sample(1:p_size,p_size)
ind2 = ind1[1:zeroNum]
beta_true = runif(p_size,R2*0.2,R2)
beta_true[ind2] = 0

X = R1*matrix(rnorm(sample_size * p_size), ncol = p_size)
X[1:sample_size,1]=1
t = 1/(1+exp(-X%*%beta_true))
tmp = runif(sample_size,0,1)
y = rep(0,sample_size)
y[t>tmp] = 1

fit = GAGA(X,y,family = "binomial", alpha = 1)
Eb = fit$beta

#Generate test samples
X_t = R1*matrix(rnorm(test_size * p_size), ncol = p_size)
X_t[1:test_size,1]=1
t = 1/(1+exp(-X_t%*%beta_true))
tmp = runif(test_size,0,1)
y_t = rep(0,test_size)
y_t[t>tmp] = 1
#Prediction
Ey = predict.GAGA(fit,newx = X_t)
cat("\n--------------------")
cat("\n err:", norm(Eb-beta_true,type="2")/norm(beta_true,type="2"))
cat("\n acc:", cal.w.acc(as.character(Eb!=0),as.character(beta_true!=0)))
cat("\n pacc:", cal.w.acc(as.character(Ey),as.character(y_t)))
cat("\n")
