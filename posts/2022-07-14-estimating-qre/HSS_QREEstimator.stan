functions {
  
  
  /*
  QRE fixed point condition
  This function is equal to zero when evaluated at the logit of the QRE probabilities
  Stan requires that functions passed to the algebra solver have a certain signature, so
  I had to add in the x_r and x_i inputs to get it to work, even though I don't
  need these
  */
  vector QREFixedPoint(vector lp, vector pars,  data real[] x_r, data int[] x_i) {
    vector[2] z;
    
    real lambda = pars[1]; // lambda is the first element of pars
    real cost = pars[2];   // attack cost is the second element of pars
    
    // convert logit probabilities into actual probabilities
    vector[2] p = 1.0 ./(1.0+exp(-lp)); 
    
    // fixed point condition for row player
    z[1] = lp[1]-lambda*(6*p[2]-3);
    // fixed point condition for column player
    z[2] = lp[2]-lambda*(6*(1-p[1])-cost-1);
    
    return z;
    
  }
}

data {
  int<lower=0> n; // number of rows of data
  vector[n] cost; // attack cost
  vector[n] H;    // number of times High Alert is played in a treatment
  vector[n] A;    // Number of times Attack is played
  vector[n] N;    // total number of actions played
  
  vector[2] prior_lambda; // prior for lambda (lognormal)
}

transformed data {
  // The algebra solver needs to have some data to pass to it, so let's greate some
  real x_r[1] = {3.14};
  int x_i[1]  = {42};
}

parameters {
  // logit choice precision parameter
  real<lower=0> lambda;
  
}

transformed parameters {
  
  /* sometimes we need to store some intermediate values from our calculations.
  In particular here, I want to store the model's predictions to plot later.
  I also store the log-likelihood, as it is useful for model evaluation ( although
  I', not doing one here).
  */
  vector[n] log_like;
  vector[2] predictions[n];
  for (tt in 1:n) {
    
    // set up some inputs for the fixed point calculation
    vector[2] lp0 = to_vector({0.0,0.0}); // initial guess
    vector[2] pars = to_vector({lambda,cost[tt]}); // parameters to pass to solver
    
    // solve QRE 
    vector[2] lp = algebra_solver(QREFixedPoint,lp0,pars,x_r,x_i);
    
    // convert logit probabilities to actual probabilities
    vector[2] p = 1.0 ./(1+exp(-lp));
    predictions[tt] = p;
    
    
    // likelihood contribution
    log_like[tt]=  H[tt]*log(p[1])+(N[tt]-H[tt])*log(1.0-p[1])
            + A[tt]*log(p[2])+(N[tt]-A[tt])*log(1.0-p[2]);
    
  }
}


model {
  
  // prior
  lambda ~ lognormal(prior_lambda[1],prior_lambda[2]);
  
  // likelihood contribution
  target+=log_like;
  
  
}

