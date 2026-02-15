
data {
  
  // sample sizes
  int<lower = 0> N1, N2;
  // first sample
  vector[N1] Y1;
  
  // prior
  vector[2] prior_mu;
  real<lower=0> prior_sigma;
  
  // cost of second sample
  real<lower=0> m;
  
  // CARA coefficient
  real r;

}


parameters {
  
  // posterior draws after first sample
  real mu1;
  real<lower=0> sigma1;
  
  // posterior draws after second sample
  real mu2;
  real<lower=0> sigma2;
  
  // data in second sample
  vector[N2] Y2;
  
  
}


model {
  
  mu1 ~ normal(prior_mu[1],prior_mu[2]);
  sigma1 ~ cauchy(0.0, prior_sigma);
  
  // likelihood of first sample
  Y1 ~ normal(mu1, sigma1);
  
  // data in second sample
  Y2 ~ normal(mu1,sigma1);
  
  // updating for 2nd sample
  // prior
  mu2 ~ normal(prior_mu[1],prior_mu[2]);
  sigma2 ~ cauchy(0.0, prior_sigma);
  // likelihood
  Y1 ~ normal(mu2, sigma2);
  Y2 ~ normal(mu2, sigma2);
  
}

generated quantities {
  
  vector[3] EVchoices;
    // take the bet now
    EVchoices[1] = (1-exp(-r*normal_rng(mu1,sigma1)))/r;
    // stop now
    EVchoices[2] = 0.0;
    // pay for more data
    EVchoices[3] = mu2>0 ? 
        (1-exp(-r*(normal_rng(mu2,sigma2)-m)))/r : (1-exp(-r*(-m)))/r;
  
  
}

