# cox
set.seed(2022)

cat("\n")
cat("Test cox GAGA\n")

p_size = 50
sample_size = 500
test_size = 1000
test_size = 1000

R1 = 3
R2 = 1
rate = 0.5 #Proportion of value zero in beta
censoringRate = 0.25 #Proportion of censoring data in observation data



#Set true beta
zeroNum = round(rate*p_size)
ind = sample(1:p_size,zeroNum)
beta_true = runif(p_size,-R2,R2)
beta_true[ind] = 0

#Generate training samples
cov_mat=matrix(1:p_size*p_size,p_size,p_size) ##covariance matrix
for(i in 1:p_size){for(j in 1:p_size) {if(i!=j){cov_mat[i,j]=0.0}else{cov_mat[i,j]=1}}}
X = R1*rmvnorm(n=sample_size, mean=rep(0, nrow(cov_mat)), sigma=cov_mat)
z = X%*%beta_true
u = runif(sample_size,0,1)
t = ((-log(1-u)/(3*exp(z)))*100)^(0.1)
cs = rep(0,sample_size)
csNum = round(censoringRate*sample_size)
ind = sample(1:sample_size,csNum)
cs[ind] = 1
t[ind] = runif(csNum,0,0.8)*t[ind]
y = cbind(t,1 - cs)
colnames(y) = c("time", "status")

#Estimation
fit = GAGA(X,y,alpha=2,family="cox")
Eb = fit$beta
# cat("\n exit itrNum: ", fit$itrNum)
#Generate testing samples
X_t = R1*rmvnorm(n=test_size, mean=rep(0, nrow(cov_mat)), sigma=cov_mat)
z = X_t%*%beta_true
u = runif(sample_size,0,1)
t = ((-log(1-u)/(3*exp(z)))*100)^(0.1)
cs = rep(0,test_size)
csNum = round(censoringRate*test_size)
ind = sample(1:test_size,csNum)
cs[ind] = 1
t[ind] = runif(csNum,0,0.8)*t[ind]
y_t = cbind(t,1 - cs)
colnames(y_t) = c("time", "status")

#Prediction
pred = predict.GAGA(fit,newx=X_t)

cat("\n--------------------")
cat("\n err:", norm(Eb-beta_true,type="2")/norm(beta_true,type="2"))
cat("\n acc:", cal.w.acc(as.character(Eb!=0),as.character(beta_true!=0)))
cat("\n Cindex:", cal.cindex(pred,y_t))
cat("\n")



