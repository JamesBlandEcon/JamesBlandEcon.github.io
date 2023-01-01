// saved as Brisket.stan
data {
  
  int n; // number of observations
  vector[n] meat;
  vector[n] grate;
  vector[n] dt;
  vector[n] dTempMeat;
  
}
transformed data {

}
parameters {
  real<lower=0> h1;     
  real<lower=0> h2;
  real          Tstall;
  real<lower=0> T_error;
}
transformed parameters {
  
}
model {
  vector[n] T_hat;
  vector[n] L;
  vector[n] dT;
  T_hat[1] = meat[1];
  L[1] = 1;
  dT[1] = 0;
  
  for (ii in 2:n) { // temperature below stall
	  if (meat[ii-1]<Tstall) {
	    T_hat[ii] = T_hat[ii-1] + h1*(grate[ii-1]-meat[ii-1])*dt[ii-1];
	    L[ii] = L[ii-1];
	  }
	  else { // temperature at or above stall
	    if (L[ii-1]>0 ){ // still some latent heat to overcome
	      L[ii] = L[ii-1]-h2*(grate[ii-1]-meat[ii-1])*dt[ii-1];
	      T_hat[ii] = T_hat[ii-1];
	    }
	    else { // sufficient latent heat to overcome stall
	      T_hat[ii] = T_hat[ii-1] + h1*(grate[ii-1]-meat[ii-1])*dt[ii-1];
	      L[ii] = L[ii-1];
	    }
	  }
	  dT[ii] = T_hat[ii]-T_hat[ii-1];
	  meat[ii] ~ normal(T_hat[ii],T_error);
	}
	
	// priors 
	h1 ~ exponential(100);
	h2 ~ exponential(100);
	Tstall ~ normal(50,10);

	T_error ~ exponential(5);
	
	
	
	//dTempMeat ~ normal(h1*(grate-meat).*dt,T_error);
  
}

generated quantities {
  
 
	
	
}
