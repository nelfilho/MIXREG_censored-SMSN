data {
  int<lower=0> n_obs;
  int<lower=0> n_cen;
  int<lower=1> p;
  real y_obs[n_obs]; 
  int<lower=2> G;
  matrix[n_obs,p] x_obs;
  matrix[n_cen,p] x_cen;
  real b;
  real trun;
}
parameters {
  vector[p] beta[G]; // G vetores de tamanho p
  //real Delta[G]; 
  //real<lower=0> gam2[G];
  real<lower=0> sigma2[G];
  real alpha[G];
  simplex[G] Theta;
  real<lower=0> f;
  real mualpha;
}
transformed parameters{
  real Delta[G]; 
  for(k in 1:G){
   Delta[k] = (sqrt(sigma2[k])*alpha[k])/sqrt(1 + alpha[k]^2);
  }
}
model {
  real ps[G];
  real pp[G];
  Theta ~ dirichlet(rep_vector(2.0, G));
  f ~ gamma(0.01,0.01);
  mualpha ~ normal(0,5);
  for(k in 1:G){
    beta[k] ~ multi_normal(rep_vector(0,p), diag_matrix(rep_vector(10,p)));
    //mu[k] ~ normal(0,10);
    sigma2[k] ~ inv_gamma(2,f);
    //gam2[k] ~ inv_gamma(2,0.01);
    //alpha[k] ~ normal(mualpha[k],10);
    alpha[k] ~ normal(mualpha,10);
    //Delta[k] ~ normal(0,10);
  }
  for(n in 1:n_obs){
    for (k in 1:G){
      ps[k] = log(Theta[k]) + skew_normal_lpdf(y_obs[n] |x_obs[n,]*beta[k] + b*Delta[k], sqrt(sigma2[k]), alpha[k]);
      //ps[k] = log(Theta[k]) + skew_normal_lpdf(y_obs[n] |x_obs[n,]*beta[k], sqrt(sigma2[k]), alpha[k]);
    }
    target += log_sum_exp(ps);
  }
    for(n in 1:n_cen){
    for (k in 1:G){
      pp[k] = log(Theta[k]) + skew_normal_lcdf(trun |x_cen[n,]*beta[k] + b*Delta[k], sqrt(sigma2[k]), alpha[k]);
      //pp[k] = log(Theta[k]) + skew_normal_lcdf(trun |x_cen[n,]*beta[k],  sqrt(sigma2[k]), alpha[k]);
    }
    target += log_sum_exp(pp);
  }
}