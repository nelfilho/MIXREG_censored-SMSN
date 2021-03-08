# Fun��o para gerar observa��es do modelo de regress�o linear com com erros nas
# vari�veis
# onde a vari�vel latente � uma mistura de normais inversas Gaussianas  
#========================================================================
# Input 
#========================================================================
# n: tamanho amostral
# alpha, beta: vetores de regress�o de  dimens�o r; 
# Omega: matriz p x p associada � escala  dos erros de observa��o, onde p=r+1; 
# mu: vetor de  par�metros associados � loca��o da vari�vel latente para cada
# sub-popula��o. dim(mu)=G, onde
# G � o n�mero de sub-popula��es e mu[i] � o par�metro de loca��o associado �
# sub-popula��o i;
# lambda: analogamente, � um vetor de dimens�o G contendo os par�metros reguladores
# de assimetria para cada sub-popula��o;
# gama, delta: constantes positivas, par�metros do fator de escala - que tem
# distribui��o Gaussiana inversa
# pii: vetor de dimens�o G com os pesos na mistura
#---------------------------------------------------------------------------
#  In�cio da fun��o
#---------------------------------------------------------------------------
reiv.NIGmix <- function(n,alpha,beta,Omega,mu,lambda,gama,delta,pii){
library(mnormt)
library(SuppDists)
p <- length(alpha) + 1
G <- length(pii)        # n�mero de componentes na mistura 
#-------------------------------------------------------------------------------
x <- matrix(0,n,1)            # vetor que conter� os valores gerados da vari�vel
                              # latente
Z <- matrix(0,n,p)            # vetor que conter� as observa��es geradas
#===============================================================================
a <- matrix(c(0,alpha))
B1 <- matrix(c(1,beta))
u <- vector()
w <- matrix(sample(G,size=n,replace=TRUE,prob=pii)) # vetor de aloca��es. Se z[i]=j
                                                 #o indiv�duo i � alocado � classe j
nu <- delta/gama
psi <- delta^2 
for (i in 1:n){
u[i] <- rinvGauss(1,nu,psi)
x[i] <- rnorm(1,mu[w[i]]+u[i]*lambda[w[i]],sqrt(u[i]))
Z[i,] <- matrix(rmnorm(1,a+B1*x[i],u[i]*Omega))
}
return(list(Z=Z,w=w,x=x))
    }
    

# Exemplo
mu <- c(-3,1)
#lambda <- c(-0.2,0.1)
#lambda<-c(-0.5,0.4)
lambda<-c(-2,1)
delta <- 0.7
gama <- 1 
pii <- c(0.4,0.6)
Omega <- diag(c(1,1,1)) 
alpha <- c(0.4,0.1)
beta <- c(0.8,0.9)
n <- 200

library(robustbase)    
y.teste <- reiv.NIGmix(n,alpha,beta,Omega,mu,lambda,gama,delta,pii)
#par(mfrow=c(1,5))
#w<-y.teste$w
#plot(density(y.teste$x[w==1]),main='Group 1')
#plot(density(y.teste$x[w==2]),main='Group 2')
#plot(density(y.teste$x),main='latent variable')
#adjbox(y.teste$x[w==1],main='Group 1')
#adjbox(y.teste$x[w==2],main='Group 2')

#Z=y.teste$Z