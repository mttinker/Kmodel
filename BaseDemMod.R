# Examine base demographic matrix for sex-structured pop model for southern sea otters --------------
lam = numeric()
fs = numeric()
ms = numeric()
R = numeric() 
PropK = seq(.01,1,by=.01)
nd = length(PropK)
for (i in 1:nd){
  pK = PropK[i]
  R[i] = (2^(-3*pK))/10 + 0.167
  fs[i] = exp(-.01*exp(1.3+1.5*pK))
  ms[i] = exp(-.01*exp(1.3+.5+1.5*pK))
  M = matrix(c(fs[i]*(1+R[i]),   0,
               fs[i]*R[i],      ms[i]),byrow = TRUE,ncol = 2)
  lam[i] = eigen(M)$values[1]; lam # Growth rate
  W=eigen(M)$vectors          # W=matrix of right eigenvectors 
  w=abs(W[,1])					      # w=stable distribution, unscaled
  ssd = w/sum(w)              # w=stable distribution, scaled
}
par(mfrow=c(2,2))
plot(PropK,log(lam),col="black",type = "l",xlab = "Proportion of K",ylab = "Log-Lambda", ylim = c(0,.2))
plot(PropK,R,col="purple",type = "l",xlab = "Proportion of K",ylab = "Reproduction", ylim = c(.15,.3))
plot(PropK,fs,col="blue",type = "l",xlab = "Proportion of K",ylab = "Female Survival", ylim = c(.75,1))
plot(PropK,ms,col="red",type = "l",xlab = "Proportion of K",ylab = "Male Survival", ylim = c(.75,1))
# Test function ------------------------------------------------------------------------------------
source("somatrixHz.r")
lam = numeric()
R = numeric()
for (i in 1:5000){
  rslt = somatrixHz(1,0,0)
  lam[i] = rslt$lam
  R[i] = rslt$vrts[1]
}
hist(lam)
hist(R)
mean(lam)
mean(R)  
