# Poisson
set.seed(2022)
library(mvtnorm)
cat("\n")
cat("Test poisson GAGA\n")
p_size = 30
sample_size=300
R1 = 1/sqrt(p_size)
R2 = 5
rate = 0.5 #Proportion of value zero in beta

#Set true beta
zeroNum = round(rate*p_size)
ind1 = sample(1:p_size,p_size)
ind2 = ind1[1:zeroNum]
beta_true = runif(p_size,0,R2)
beta_true[ind2] = 0

X = R1*matrix(rnorm(sample_size * p_size), ncol = p_size)
X[1:sample_size,1]=1
y = rpois(sample_size,lambda = as.vector(exp(X%*%beta_true)))
y = as.vector(y)

Eb = GAGA(X,y,alpha = 2,family="poisson")

cat("\n--------------------")
cat("\n err:", norm(Eb-beta_true,type="2")/norm(beta_true,type="2"))
cat("\n acc:", calculate.w.accuracy(as.character(Eb!=0),as.character(beta_true!=0)))
cat("\n")
