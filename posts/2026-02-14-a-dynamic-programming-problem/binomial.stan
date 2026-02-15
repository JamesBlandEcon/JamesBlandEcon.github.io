
functions {
  
  // CARA utility function
  real ufun(real x, real r) {
    
    return r==0 ? x : (1.0-exp(-r*x))/r;
  }
  
  // inverse CARA utility function 
  real iufun(real u, real r) {
    
    /* some working
    
    u = (1-exp(-rx))/r
    ru = 1-exp(-rx)
    -rx = log(1-ru)
    x = -log(1-ru)/r
    */
    
    return r==0 ? u : -log(1.0-r*u)/r;
  }
  
}

data {
  // first and second sample sizes
  int<lower=0> N1, N2;
  
  // number of heads observed in first sample
  int<lower=0,upper=N1> nHeads1;
  
  // prior for theta (Beta)
  vector[2] prior_theta;
  
  // benefit/cost values
  real<lower=0> W, L, m;
  
  // CARA coefficient
  real r; 
  
}


parameters {
  
  // probability of flipping heads
  real<lower=0, upper=1> theta;
  
}

model {
  
  // prior
  theta ~ beta(prior_theta[1],prior_theta[2]);
  
  // likelihood of first sample
  nHeads1 ~ binomial(N1, theta);
  
}

generated quantities {
  
  // posterior after observing first sample
  vector[2] posterior1 =[prior_theta[1]+nHeads1, 
          prior_theta[2]+N1-nHeads1]';
  
  // prediction after observing first sample
  real prHeads1 = posterior1[1]/sum(posterior1);
  
  // second sample
  int nHeads2 = binomial_rng(N2,theta);
  
  // posterior for beta after observing the second sample
  vector[2] posterior2 = [prior_theta[1]+nHeads1+nHeads2, 
          prior_theta[2]+N1+N2-nHeads1-nHeads2]';
  
  
  // prediction for next coin flip conditional on seeing both samples
  real prHeads2 = posterior2[1]/sum(posterior2);
  
  vector[3] EUchoices;
  
  // take the bet now
  EUchoices[1] = ufun(W,r)*prHeads1+ufun(-L,r)*(1.0-prHeads1);
  // stop now
  EUchoices[2] = ufun(0.0,r);
  // pay for more data
  EUchoices[3] = fmax(ufun(-m,r),ufun(W-m, r)*prHeads2+ufun(-L-m,r)*(1.0-prHeads2));
  
  vector[3] CEchoices;
    for (cc in 1:3) {
      
      CEchoices[cc] = iufun(EUchoices[cc],r);
      
    }

}

