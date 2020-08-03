data {
  int p;
  int G;
  int<lower=0> n_obs;
  int<lower=0> n_cen;
  real y_obs[n_obs]; 
  real<upper=min(y_obs)> trun;
  matrix[n_obs,p] x_obs;
  matrix[n_cen,p] x_cen;
}
parameters {
  //matrix[p,G] beta;
  vector[p] beta[G]; // G vetores de tamanho p
  real<lower=0> sigma[G];
  simplex[G] Theta;
  //simplex[G] phi;
  //real<lower=0> kappa;
}
//transformed parameters {
  //vector[G] alpha = kappa * phi;
  //vector[n_obs] mu_obs[G];
  //vector[n_cen] mu_cen[G];
  //for(k in 1:G){
  //mu_obs[k] = x_obs*beta[k];
  //mu_cen[k] = x_cen*beta[k];
  //}
//}
model {
  real ps[G];
  real pp[G];
  sigma ~ cauchy(0, 2.5);
  Theta ~ dirichlet(rep_vector(2, G));
  //Theta ~ dirichlet(alpha);
  for(k in 1:G){
    beta[k]~multi_normal(rep_vector(0,p),diag_matrix(rep_vector(100,p)));
      }
  for(n in 1:n_obs){
    for(k in 1:G){
        ps[k] = log(Theta[k]) + normal_lpdf(y_obs[n] | x_obs[n,]*beta[k], sigma[k]); //mistura
      }
    target += log_sum_exp(ps);
    //target += log_mix(Theta[1], normal_lpdf(y_obs[n] | x_obs[n,]*beta[1], sigma[1]),normal_lpdf(y_obs[n] | x_obs[n,]*beta[2], sigma[2]));
  }
  //target += log_mix(Theta[1], normal_lpdf(y_obs | mu_obs[1], sigma[1]),normal_lpdf(y_obs | mu_obs[2], sigma[2]));
  for(n in 1:n_cen){
    for(k in 1:G){
        pp[k] = log(Theta[k]) + normal_lcdf(trun | x_cen[n,]*beta[k], sigma[k]); // truncamento
      }
    target += log_sum_exp(pp);
    //target += log_mix(Theta[1],normal_lcdf(trun | x_cen[n,]*beta[1], sigma[1]),normal_lcdf(trun | x_cen[n,]*beta[2], sigma[2]));
  }
  //target += log_mix(Theta[1],normal_lcdf(trun | mu_cen[1], sigma[1]),normal_lcdf(trun | mu_cen[2], sigma[2]));
}