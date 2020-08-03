library(dplyr) 
library(rstan)
library(ggplot2)
library(ggthemes)

N <- 200
p <- 2
seed <- 420

#121 betas bons thetas nao muito
#122 um beta ruim o resto bom
set.seed(seed)
x <- matrix(c(rep(1,N),rnorm(N,0,1)),N,p)

#tree groups-----
# G <- 3
# beta <- matrix(c(2,3,2,6,2,10),p,G)
# sigma <- c(0.8, 0.5, 0.7)
# Theta <- c(.2, .3, .5)
# nu <- 2

#two groups-----
G <- 2
beta <- matrix(c(2,3,2,10),p,G)
sigma <- c(0.5, 0.5)
Theta <- c(.3, .7)
nu <- c(2,2)

source("RanCensoredRegressionMixture - T.r")
set.seed(seed)
data = RcensRegMix.T(N,x,beta,sigma,nu,Theta,truncated="left",trunc.quantile=0.1)

plot(x[,2],data$Y,col=data$z)
abline(beta[,1])
abline(beta[,2],col=2)
#abline(beta[,3],col=3)

#--------------------------------------------------------------------------------------
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7')

stanDso = rstan::stan_model(file='StanCensoredRegressionMixture - T.stan')

ind_cen = which(data$Y==data$trunc)
ind_obs = which(data$Y!=data$trunc)
n_obs = length(ind_obs)
n_cen = length(ind_cen)
y_obs = data$Y[ind_obs]

## take a subset of the data
lSData = list(p=p,n_obs=n_obs,n_cen=n_cen, y_obs=y_obs, G=data$n_groups, x_obs=data$X[ind_obs,],x_cen=data$X[ind_cen,],trun=as.numeric(data$trunc))

## give initial values function to stan
fit.stan = sampling(stanDso, data=lSData,control = list(adapt_delta = 0.9999, max_treedepth = 15), iter=5000, thin=1 ,chains=1, 
                    warmup=2000,refresh=1)
print(fit.stan) 
pairs(fit.stan,pars=c("beta","sigma","Theta","nu"))
#pairs(fit.stan,pars=c("beta","nu"))
save(list = ls(),file = "rdata.RData")

#pairs(fit.stan,pars=c("beta","sigma","Theta","nu"))
