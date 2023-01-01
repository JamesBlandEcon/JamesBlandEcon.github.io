// saved as AFLAngleT.stan
data {
  vector[2] priork;
  vector[2] priornu;
  vector[3] priorDistance;
  int n;
  vector[n] goal;
  vector[n] behind;
  vector[n] MissShort;
  vector[n] bounds;
  vector[n] aGoal;
  vector[n] aLBehind;
  vector[n] aRBehind;
  vector[n] dist;
}
transformed data {
	vector[n] ones;
	for (ii in 1:n) {
		ones[ii] = 1;
	}
}
parameters {
  real<lower=0> k;     
  real<lower=0> nu;
  real ldm;
  real<lower=0> ldsd;
}
transformed parameters {
  
}
model {
	vector[n] pGoal;
	vector[n] pLBehind;
	vector[n] pRBehind;
	vector[n] pDistance;
	
	
	// prior
	k ~ lognormal(priork[1],priork[2]);
	nu ~ lognormal(priornu[1],priornu[2]);
	ldm ~ normal(priorDistance[1],priorDistance[2]);
	ldsd ~ exponential(priorDistance[3]);
	

	
	// student_t_cdf(reals y, reals nu, reals mu, reals sigma)
	for (ii in 1:n) {
		pGoal[ii]       = student_t_cdf(aGoal[ii] / (2.0),nu,0.0,k)-student_t_cdf(-aGoal[ii] / (2.0),nu,0.0,k);
        pLBehind[ii]    = student_t_cdf(-aGoal[ii]/(2.0),nu,0.0,k)-student_t_cdf(-aGoal[ii]/(2.0)-aLBehind[ii],nu,0.0,k);
        pRBehind[ii]    = student_t_cdf(aRBehind[ii]+aGoal[ii]/(2.0),nu,0.0,k)-student_t_cdf(aGoal[ii]/(2.0),nu,0.0,k);
		pDistance[ii] = normal_cdf(log(dist[ii]),ldm,ldsd);
	}
	
	
	
	// likelihood contribution for goals
	target += goal .* log(pGoal .* (1.0-pDistance));
	// likelihood contribution for behinds
	target += behind .* log((pLBehind+pRBehind)  .* (1.0-pDistance));
	// likelihood contribution for miss, didn't make the distance
	target += MissShort .* log(pDistance);
	target += bounds .* log((1.0-pDistance) .* (1.0-pGoal - pRBehind-pLBehind));
	
	// old code for G/B/other
	// likelihood contribution for goals
	//target += goal .* log(pGoal);
	// likelihood contribution for behinds
	//target += behind .* log((pLBehind+pRBehind));
	// likelihood contribution for the rest
	//target += (MissShort+bounds) .* log(1.0-pLBehind-pRBehind-pGoal);
}

generated quantities {
	
}
