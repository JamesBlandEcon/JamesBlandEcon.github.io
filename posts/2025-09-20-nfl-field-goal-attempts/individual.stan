
data {
  int<lower=0> N;
  vector[N] dist;
  int good[N];
  
  vector[2] prior_sigma;
  vector[2] prior_nu;
  
  int<lower=0,upper=1> UseData;
  
}

transformed data {
  
  
  real goalpost_width = 18.5/3.0; // in yards
  
  vector[N] angle = atan(goalpost_width/(2.0*dist));
  
  
}

parameters {
  real<lower=0> sigma;
  real<lower=0> nu;
  
}

model {
  
  if (UseData==1) {
  
    vector[N] pr_angle;
  
    for (ii in 1:N) {
      pr_angle[ii] = 2.0*student_t_cdf(angle[ii], nu, 0.0, sigma)-1.0;
    }
  
    target += bernoulli_lpmf(good | pr_angle);
  }
  
  target += lognormal_lpdf(sigma | prior_sigma[1],prior_sigma[2]);
  target += lognormal_lpdf(nu | prior_nu[1],prior_nu[2]);
  
  
}

generated quantities {
  
  // predicted probability
  
  vector[100] pr;
  for (ii in 1:100) {
    
    pr[ii] = 2.0*student_t_cdf(atan(goalpost_width/(2.0*ii)), nu, 0.0, sigma)-1.0;
    
  }
  
}

