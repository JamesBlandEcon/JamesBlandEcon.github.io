
data {
  int<lower=0> N; // number of observations
  int<lower=0> n; //number of subjects
  int y[N]; // =1 iff q-lottery is chosen
  matrix[N,4] Q; // probabilities for Q-lottery
  matrix[N,4] P; // probabilities for P-lottery
  int<lower=0> id[N];
  real prizes[4];
  real prior_r[4];
  real prior_l[4];
  real prior_corr[2];
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector[2] mu; 
  vector<lower=0>[2] sigma;
  
  real<lower=0,upper=1> trans_corr;
  vector[2] theta[n];
}

transformed parameters {
  real<lower=-1,upper=1> corr = 2*trans_corr-1;
}


model {
  matrix[2,2] SIGMA;
  
  real DU[N] = rep_array(0.0,N);
  
  
  // priors
  mu[1] ~ normal(prior_r[1],prior_r[2]);
  mu[2] ~ normal(prior_l[1],prior_l[2]);
  sigma[1] ~ lognormal(prior_r[3],prior_r[4]);
  sigma[2] ~ lognormal(prior_l[3],prior_l[4]);
  trans_corr ~ beta(prior_corr[1],prior_corr[2]);
  
  
  // likelihood
  for (ii in 1:N) {
    real r = exp(theta[id[ii]][1]);
    real l = exp(theta[id[ii]][2]);
    
    for (kk in 1:4) {
    
      DU[ii] = DU[ii] + (Q[ii,kk]-P[ii,kk])*prizes[kk]^r; 
    
    
  }
  
  DU[ii] = l*DU[ii];
}
  y ~ bernoulli(inv_logit(DU));
  
  SIGMA[1,1] = sigma[1]^2;
  SIGMA[2,2] = sigma[2]^2;
  SIGMA[1,2] = sigma[1]*sigma[2]*corr;
  SIGMA[2,1] = sigma[1]*sigma[2]*corr;
  
  for (ii in 1:n) {
    theta[ii] ~ multi_normal(mu,SIGMA);
  }
    

}