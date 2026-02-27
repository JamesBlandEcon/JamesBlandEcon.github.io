
functions {
  
  int which_max(vector x) {
    
    int d = dims(x)[1];
    
    real m = -99999999.999999999;
    
    int i = -1;
    
    for (tt in 1:d) {
      
      i = x[tt]>m ? tt : i;
      m = x[tt]>m ? x[tt] : m;
      
    }
    
    return i;
    
  }
  
  vector ufun(vector x, real r) {
    
    return (1.0-exp(-r*x))/r;
    
    
  }
  
}

data {
  int<lower=0> N;
  
  matrix[N,3] pL, pR;
  
  vector[3] prizes;
  
  array[N] int choiceR;
  
}

transformed data {
  
  matrix[N,3] dPR = pR-pL;
  
  int nprobs = 99;
  vector[nprobs] probs;
  for (pp in 1:nprobs) {
    probs[pp] = (pp+0.0)/100.0;
  }
  int nreports = 9999;
  vector[nreports] reports;
  for (rr in 1:nreports) {
    reports[rr] = (rr+0.0)/(10000.0);
  }
  
}

parameters {
  
  real r; // CARA coefficient
  real<lower=0> lambda; // choice precision
  
}

transformed parameters {
  
   // utility of each prize
  vector[3] u = (1.0-exp(-r*prizes))/r;
  // contextual utility normalization
  u = u/(u[3]-u[1]); 
  
}


model {
  
  // prior ---------------------------------------------------------------------
  
  r ~ normal(0,1);
  lambda ~ lognormal(0,5);
  
  // likelihood ----------------------------------------------------------------
  
  choiceR ~ bernoulli_logit(lambda*dPR*u);
  
}

generated quantities {
  
  real choicePR = inv_logit(lambda*[-0.5, 1, -0.5]*u);
  
  
  // determine the optimal report given r
  vector[nprobs] REPORT;
  
  for (pp in 1:nprobs) {
    
    real p = probs[pp];
    
    vector[nreports] U = p*ufun(1+log2(reports),r)+(1.0-p)*ufun(1+log2(1.0-reports),r);
    int ii = which_max(U);
    
    REPORT[pp] = reports[ii];
    
  }
  
}
