# Gaussian
set.seed(2022)
cat("\n")
cat("Test LM GAGA\n")
p_size = 30
sample_size=300
R1 = 3
R2 = 2
rate = 0.5 #Proportion of value zero in beta

#Set true beta
zeroNum = round(rate*p_size)
ind1 = sample(1:p_size,p_size)
ind2 = ind1[1:zeroNum]
beta_true = runif(p_size,0,R2)

X = R1*matrix(rnorm(sample_size * p_size), ncol = p_size)
y=X%*%beta_true + rnorm(sample_size,mean=0,sd=2)

Eb = GAGA(X,y,alpha = 3,family="gaussian")

cat("\n--------------------")
cat("\n err:", norm(Eb-beta_true,type="2")/norm(beta_true,type="2"))
cat("\n acc:", calculate.w.accuracy(as.character(Eb!=0),as.character(beta_true!=0)))
cat("\n")






