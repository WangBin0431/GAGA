# Poisson
set.seed(2022)

cat("\n")
cat("Test poisson GAGA\n")
p_size = 30
sample_size=300
test_size = 1000
R1 = 1/sqrt(p_size)
R2 = 5
rate = 0.5 #Proportion of value zero in beta

#Set true beta
zeroNum = round(rate*p_size)
ind = sample(1:p_size,zeroNum)
beta_true = runif(p_size,0,R2)
beta_true[ind] = 0

X = R1*matrix(rnorm(sample_size * p_size), ncol = p_size)
X[1:sample_size,1]=1
y = rpois(sample_size,lambda = as.vector(exp(X%*%beta_true)))
y = as.vector(y)

fit = GAGA(X,y,alpha = 2,family="poisson")
Eb = fit$beta

#Generate test samples
X_t = R1*matrix(rnorm(test_size * p_size), ncol = p_size)
X_t[1:test_size,1]=1
y_t = rpois(test_size,lambda = as.vector(exp(X_t%*%beta_true)))
y_t = as.vector(y_t)
#Prediction
Ey = predict.GAGA(fit,newx = X_t)

cat("\n--------------------")
cat("\n err:", norm(Eb-beta_true,type="2")/norm(beta_true,type="2"))
cat("\n acc:", cal.w.acc(as.character(Eb!=0),as.character(beta_true!=0)))
cat("\n perr:", norm(Ey-y_t,type="2")/sqrt(sample_size))
cat("\n")
