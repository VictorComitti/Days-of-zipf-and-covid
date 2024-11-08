rm(list=ls())
library(gravity)
library(broom)
library(lmtest)
library(janitor)
library(sads)

set.seed(1234)

alpha<-1.0

sigma2<-1
n=50
L=1000
gamma<-0.0001
coefOLS<-NULL
coefPPML<-NULL
for(l in 1:L){
  eta=NULL
  zipf<-rzipf(60000, n, alpha)
  C<-tabyl(zipf, sort = TRUE)$n
  for(i in 1:n){
    eta[i]<-rlnorm(1, -gamma*C[i]*sigma2/2, ((gamma*C[i])^2)*sigma2)
  }
  lC<-log(C)
  R<-(10000/(C^alpha))*eta
  lR<-log(R-0.5)
  data<-as.data.frame(cbind(R, C))
  coefOLS[l]<-coef(lm(lR~lC))[2]
  coefPPML[l]<-coef(ppml("R", "C",  additional_regressors = c(), data=data))[2]
  print(l)
}


round(mean(coefOLS), 2)
round(mean(coefPPML),2)

#length(C)

rbiasOLS<-mean(((abs(coefOLS)-alpha)/alpha)*100)
rbiasPPML<-mean(((abs(coefPPML)-alpha)/alpha)*100)
round(rbiasOLS, 2)
round(rbiasPPML,2)

