rm(list = ls())

source("rSMSN.r")

library(sn)
library(ggplot2)
library(rstan)
library(dplyr) 
library(ggthemes)

N <- 200
G <-2
p <-2

seed <- 14
set.seed(seed)
x <- matrix(c(rep(1,N),abs(runif(N,1,6))),N,p)

beta <- matrix(c(2,3,2,5),p,G)
# mu <- c(3, 6)
# gam2 <- c(0.1, 0.1)
# Delta <- c(-2,2)
Theta <- c(.7, .3)
sigma2 <- c(0.8, 0.8)
alpha <- c(10,-10)

k1 <- 1
b <- -sqrt(2/pi)*k1

# sigma2=Delta^2+gam2;
delta <- alpha / sqrt(1 + alpha^2)
Delta = delta*sqrt(sigma2)
# gam2 = sigma^2*(1-delta^2); 

# var = (sigma2)*(1-2*(delta^2)/pi)
# mean = mu + sqrt(sigma2)*delta*sqrt(2/pi)
# skewness = ((4-pi)/2)*((delta*sqrt(2/pi))^3)/((1-(2*delta^2/pi))^(3/2))

#for(i in 11:100){
#print(i)  
#set.seed(i+10)
set.seed(34)
z <- sample(1:G, size = N, prob = Theta, replace = T)
#y <- gen.SN(N,mu[z],gam2[z],Delta[z])
#y <- rSMSN(N, mu[z], sigma[z]^2, alpha[z],dist="SN")

y=NULL
for(j in 1:N){
  #y[j] <- rSMSN(1,x[j,]%*%beta[,z[j]] + b*Delta[z[j]], sigma[z[j]]^2, alpha[z[j]], dist="SN")
  #y[j] <- rSMSN(1,x[j,]%*%beta[,z[j]] , sigma2[z[j]], alpha[z[j]], dist="SN")
  y[j] <- x[j,]%*%beta[,z[j]] + rSMSN(1, b*Delta[z[j]], sigma2[z[j]], alpha[z[j]], dist="SN")
}

yc <- y

trun=quantile(yc,0.1)
yc[y<trun]=trun
ind_cen = which(yc==trun)
ind_obs = which(yc!=trun)
n_obs = length(ind_obs)
n_cen = length(ind_cen)

par(mfrow=c(1,2))
boxplot(y~z,ylim=c(min(y),max(y)),main="Data Real")
boxplot(yc~z,ylim=c(min(y),max(y)),main="Data Censured")
abline(h=trun)

stanDso = rstan::stan_model(file='StanMixtureRegressionCensured-SN.stan')
lSData = list(n_obs=n_obs,n_cen=n_cen,trun=trun,y_obs=yc[ind_obs],x_obs=x[ind_obs,],x_cen=x[ind_cen,],G=G,y=y,mualpha=alpha,p=p)

#ini = function(){list(mu=mu,alpha=alpha,sigma2=sigma2,Theta=Theta)}
#ini = function(){list(mu=mu,Delta=Delta,gam2=gam2,Theta=Theta)}

## give initial values function to stan
fit.stan = sampling(stanDso, data=lSData, iter=10000,chains=1,thin=1, cores=4,control = list(adapt_delta = 0.99,max_treedepth = 15))

print(fit.stan)
pairs(fit.stan,pars=c("sigma2","alpha"))

#save(list = ls(),file = "rdata.RData")
