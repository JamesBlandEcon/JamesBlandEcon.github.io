
data {
  int<lower=0> N;
  vector[N] count;
  vector[N] angle_x_lb;
  vector[N] angle_x_ub;
  vector[N] angle_y_lb;
  vector[N] angle_y_ub;
  vector[N] miss;
  
  vector[2] whLetter;
  
  vector[2] prior_sigma;
  
  int ndarts;
  int dartid[N];
  
  int nblasters;
  int blasterid[N];
  
  int npeople;
  int personid[N];

}


parameters {
  
  vector[ndarts] Bdart_sigma_x;
  vector[nblasters] Bblaster_sigma_x;
  vector[npeople] Bperson_sigma_x;
  
  
  vector[ndarts] Bdart_sigma_y;
  vector[nblasters] Bblaster_sigma_y;
  vector[npeople] Bperson_sigma_y;
  
  
  

}

transformed parameters {
  
   
  
  
}


model {
  
  vector[N] sigma_x = exp(
    Bdart_sigma_x[dartid]+Bblaster_sigma_x[blasterid]+Bperson_sigma_x[personid]
  );
 
  vector[N] sigma_y = exp(
    Bdart_sigma_y[dartid]+Bblaster_sigma_y[blasterid]+Bperson_sigma_y[personid]
  );
 
  
    vector[N] lpr;
    for (ii in 1:N) {
     lpr[ii] = log(normal_cdf(angle_x_ub[ii], 0.0, sigma_x[ii])-normal_cdf(angle_x_lb[ii], 0.0, sigma_x[ii]))
                    +
                    log(normal_cdf(angle_y_ub[ii],0.0, sigma_y[ii])-normal_cdf(angle_y_lb[ii],  0.0, sigma_y[ii]))
                    ;
    }                
      lpr = lpr.*(1.0-miss)+log(1.0-exp(lpr)).*miss;
    
      target += count.*lpr;
  
  
  target += normal_lpdf(Bdart_sigma_x | prior_sigma[1],prior_sigma[2]);
  target += normal_lpdf(Bperson_sigma_x | prior_sigma[1],prior_sigma[2]);
  target += normal_lpdf(Bblaster_sigma_x | prior_sigma[1],prior_sigma[2]);
  target += normal_lpdf(Bdart_sigma_y | prior_sigma[1],prior_sigma[2]);
  target += normal_lpdf(Bperson_sigma_y | prior_sigma[1],prior_sigma[2]);
  target += normal_lpdf(Bblaster_sigma_y | prior_sigma[1],prior_sigma[2]);
  
}

generated quantities {
  
  /* Predictions for a distance of 5m on the large target
  */
  
  vector[nblasters] pr_blaster;
  vector[ndarts] pr_dart;
  vector[npeople] pr_person;
  
  {
    real dist = 5.0;
    
    real angle_x = atan(whLetter[1]/2/dist);
    real angle_y = atan(whLetter[2]/2/dist);
    
    real sigma_x;
    real nu_x;
    real sigma_y;
    real nu_y;
    
    for (bb in 1:nblasters) {
      
      sigma_x = exp(Bblaster_sigma_x[bb]+mean(Bdart_sigma_x)+mean(Bperson_sigma_x));
      sigma_y = exp(Bblaster_sigma_y[bb]+mean(Bdart_sigma_y)+mean(Bperson_sigma_y));
      
      pr_blaster[bb] = normal_cdf(angle_x,  0.0, sigma_x)-normal_cdf(-angle_x_lb, 0.0, sigma_x);
      
    }
    
    for (dd in 1:ndarts) {
      
      sigma_x = exp(Bdart_sigma_x[dd]+mean(Bblaster_sigma_x)+mean(Bperson_sigma_x));
      sigma_y = exp(Bdart_sigma_y[dd]+mean(Bblaster_sigma_y)+mean(Bperson_sigma_y));
      
      pr_dart[dd] = normal_cdf(angle_x,  0.0, sigma_x)-normal_cdf(-angle_x_lb, 0.0, sigma_x);
    }
    
    for (pp in 1:npeople) {
      
      sigma_x = exp(Bperson_sigma_x[pp]+mean(Bdart_sigma_x)+mean(Bblaster_sigma_x));
      sigma_y = exp(Bperson_sigma_y[pp]+mean(Bdart_sigma_y)+mean(Bblaster_sigma_y));
      
      pr_person[pp] = normal_cdf(angle_x,  0.0, sigma_x)-normal_cdf(-angle_x_lb, 0.0, sigma_x);
    }
    
  }
  
  
}

