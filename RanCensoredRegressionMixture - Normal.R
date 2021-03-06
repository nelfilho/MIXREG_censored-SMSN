#rm(list = ls())
#r.N - function for generation of regression models in the SMSM family - N distribution

#input
### N : number of observations
### x : planning matrix
### beta: parameters matrix pxG
### G: number of groups in model
### sigma: parameter vector Gx1
### Theta: probability vector for the G groups, Gx1
### truncated: truncation side "left" or "right"
### trunc.quantile: truncating above or below the specified quantil, iterval 0 - 1
### trunc : exact truncation value

#output
### list containing: response variable, planning matrix, and truncation value
library(mnormt)

RcensRegMix.N <- function(N,x,beta,sigma,Theta,truncated="left",trunc.quantile=0.1,trunc=NULL){
  
  n_groups <- ncol(beta)
  p <- nrow(beta)
  
  # Draw which model each belongs to
  z <- sample(1:n_groups, size = N, prob = Theta, replace = T)
  
  # Some white noise
  epsilon <- rnorm(N)
  y <- apply(x*t(beta[,z]),1,"sum") + sigma[z]*epsilon
  
  if(is.numeric(trunc) & is.numeric(trunc.quantile)){print("ERRO")}else{
    if(is.numeric(trunc)){
      if(truncated == "left"){
        y[y<trunc] = trunc
      }
      if(truncated == "right"){
        print("ok")
        y[y>trunc] = trunc
      }  
      return(list(Y=y,X=x,trunc = trunc ))
    }
    if(is.numeric(trunc.quantile)){
      if(truncated == "left"){
        qt = quantile(y,trunc.quantile)
        y[y< quantile(y,trunc.quantile)] = qt
      }
      if(truncated == "right"){
        qt = quantile(y,trunc.quantile)
        y[y> quantile(y,trunc.quantile)] = qt
      }  
      return(list(Y=y,X=x,trunc = qt,z=z,n_groups=n_groups))
    }
  }
  
}

# N <- 200
# p <- 2
# n_groups <- 3
# x <- matrix(c(rep(1,N),runif(N,1,12)),N,p)
# beta <- matrix(c(2,3,2,6,2,10),p,n_groups)
# sigma <- c(0.8, 0.5, 0.7)
# Theta <- c(.2, .3, .5)
# 
# data = RcensRegMix.N(N,x,beta,sigma,Theta,truncated="left",trunc.quantile=0.1)
# 
# plot(x[,2],data$Y,col=data$z)
# abline(beta[,1])
# abline(beta[,2])
# abline(beta[,3])
