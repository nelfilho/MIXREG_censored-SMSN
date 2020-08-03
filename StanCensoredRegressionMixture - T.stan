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
  vector[p] beta[G]; // G vetores de tamanho p
  real<lower=0> sigma[G];
  simplex[G] Theta;
  real<lower=0> nu;
}
model {
  real ps[G];
  real pp[G];
  sigma ~ cauchy(0, 5);
  Theta ~ dirichlet(rep_vector(2, G));
  nu ~ exponential(2);
  for(k in 1:G){
      beta[k]~multi_normal(rep_vector(0,p),diag_matrix(rep_vector(100,p)));
      }
  for(n in 1:n_obs){
    for(k in 1:G){
        ps[k] = log(Theta[k]) + student_t_lpdf(y_obs[n] |nu,x_obs[n,]*beta[k],sigma[k]); //mistura
      }
    target += log_sum_exp(ps);
  }
  for(n in 1:n_cen){
    for(k in 1:G){
        pp[k] = log(Theta[k]) + student_t_lcdf(trun |nu,x_cen[n,]*beta[k],sigma[k]); // truncamento
      }
    target += log_sum_exp(pp);
  }
}
