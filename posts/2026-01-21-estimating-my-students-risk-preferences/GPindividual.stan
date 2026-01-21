
data {
  int<lower=0> N; // number of observations
  array[N] int<lower=0, upper=25> y; // integer amount invested
  
  vector[N] Return; // high return of risky asset (chi in blog post)
  
  vector[2] prior_r; // prior for r (normal, mean and sd)
  vector[2] prior_lambda; // prior for lambda (lognormal, mean and sd)
  
}

transformed data {
  
  real low_return = 0.2; // low return of risky asset
  
  // the 26 numbers participants could choose between
  vector[26] y_grid;
  
  for (yy in 1:26) {
    y_grid[yy] = yy-1;
  }
  
}

parameters {
  
  real r; // CRRA risk aversion
  real<lower=0> lambda; // logit choice precision
  
}

model {
  
  // prior
  r ~ normal(prior_r[1],prior_r[2]);
  lambda ~ lognormal(prior_lambda[1],prior_lambda[2]);
  
  
  // likelihood
  for (ii in 1:N) { // loop over each round of the task
    
    // utility for each of the 26 possible choices
    vector[26] U = 0.5*pow(25.0-y_grid+low_return*y_grid,1.0-r)/(1.0-r)+0.5*pow(25.0-y_grid+Return[ii]*y_grid,1.0-r)/(1.0-r);
    
    target += categorical_logit_lpmf(y[ii]+1 | lambda*U);
    
  }
  
}

generated quantities {
  
  vector[N] prediction;
  
  for (ii in 1:N) {
    // utility for each of the 26 possible choices
    vector[26] U = 0.5*pow(25.0-y_grid+low_return*y_grid,1.0-r)/(1.0-r)+0.5*pow(25.0-y_grid+Return[ii]*y_grid,1.0-r)/(1.0-r);
    // choice probabilities for each of the 26 possible choices
    vector[26] p = softmax(lambda*U);
    
    // mean prediction. 
    prediction[ii] = p'*y_grid;
    
  }
  
  
}

