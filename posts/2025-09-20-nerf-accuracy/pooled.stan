
data {
  int<lower=0> N; // number of rows of data
  vector[N] count; // count of the number of times that outcome occurred
  vector[N] angle_x_lb; // lower bound for horizontal angle
  vector[N] angle_x_ub; // upper bound of horizontal angle
  vector[N] angle_y_lb; // lower bound of vertical angle
  vector[N] angle_y_ub; // upper bound of vertical angle
  vector[N] miss; // =1 if the netire target was missed, zero otherwise
  
  vector[2] whLetter; // dimensions of a letter-sized piece of paper, in meters. Long side is element 1
  
  // priors for parameters, both log-normal
  vector[2] prior_sigma;
  vector[2] prior_nu;
  
  // indicator to use the data or not (useful for prior calibration)
  int<lower=0,upper=1> UseData;
}


parameters {
  
  real<lower=0> sigma_x; // scale in horizontal direction
  real<lower=0> nu_x; // df in horizontal direction
  
  real<lower=0> sigma_y; // scale in vertical direction
  real<lower=0> nu_y; // df in vertical direction
  
}

transformed parameters {
  
  /* Compute log probability of each outcome. I did this in the transformed parameters 
  block because I was debugging the code and wanted to see what the predictions were.
  */
   vector[N] lpr;
    for (ii in 1:N) {
     lpr[ii] = log(student_t_cdf(angle_x_ub[ii], nu_x, 0.0, sigma_x)-student_t_cdf(angle_x_lb[ii], nu_x, 0.0, sigma_x))
                    +
                    log(student_t_cdf(angle_y_ub[ii], nu_y, 0.0, sigma_y)-student_t_cdf(angle_y_lb[ii], nu_y, 0.0, sigma_y))
                    ;
    }                
    // update for misses
      lpr = lpr.*(1.0-miss)+log(1.0-exp(lpr)).*miss;
      
      vector[N] prob = exp(lpr);
  
  
}


model {
  
  // if we are using the data, increment the likelihood
  if (UseData==1) {
      target += count.*lpr;
  }
  // priors
  target += lognormal_lpdf(sigma_x | prior_sigma[1],prior_sigma[2]);
  target += lognormal_lpdf(sigma_y | prior_sigma[1],prior_sigma[2]);
  target += lognormal_lpdf(nu_x | prior_nu[1],prior_nu[2]);
  target += lognormal_lpdf(nu_y | prior_nu[1],prior_nu[2]);
  
}

generated quantities {
  
  /* Make some predictions. Here we will compute the probability of hitting a 
  letter-sized piece of paper and a half-letter-size piece of paper at increments 
  of 10cm from the target. I.e. pr_letter[20] is the probability of hitting a 
  letter-sized piece of paper from 2m away.
  */
  vector[200] pr_letter;
  vector[200] pr_halfletter;
  
  
  
  for (ii in 1:200) {
    
    real dist = (ii+0.0)/10.0;
    
    real angle_x = atan(whLetter[1]/2/dist);
    real angle_y = atan(whLetter[2]/2/dist);
    
    pr_letter[ii] = exp(log(student_t_cdf(angle_x, nu_x, 0.0, sigma_x)-student_t_cdf(-angle_x, nu_x, 0.0, sigma_x))
                    +
                    log(student_t_cdf(angle_y, nu_y, 0.0, sigma_y)-student_t_cdf(-angle_y, nu_y, 0.0, sigma_y))
                    )
                    ;
                    
    angle_x = atan(whLetter[1]/4/dist);
    angle_y = atan(whLetter[2]/4/dist);
    
    pr_halfletter[ii] = exp(log(student_t_cdf(angle_x, nu_x, 0.0, sigma_x)-student_t_cdf(-angle_x, nu_x, 0.0, sigma_x))
                    +
                    log(student_t_cdf(angle_y, nu_y, 0.0, sigma_y)-student_t_cdf(-angle_y, nu_y, 0.0, sigma_y))
                    )
                    ;
    
    
  }
  
  
}

