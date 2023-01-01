// saved as AFLAngleT.stan
data {
  vector[2] priork;
  vector[2] priornu;
  int n;
  vector[n] goal;
  vector[n] behind;
  vector[n] aGoal;
  vector[n] aLBehind;
  vector[n] aRBehind;
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
}
transformed parameters {
  
}
model {
	vector[n] pGoal;
	vector[n] pLBehind;
	vector[n] pRBehind;
	
	
	// prior
	k ~ lognormal(priork[1],priork[2]);
	nu ~ lognormal(priornu[1],priornu[2]);
	

	
	// student_t_cdf(reals y, reals nu, reals mu, reals sigma)
	for (ii in 1:n) {
		pGoal[ii]   	= student_t_cdf(aGoal[ii] / (2.0),nu,0.0,k)-student_t_cdf(-aGoal[ii] / (2.0),nu,0.0,k);
		pLBehind[ii] 	= student_t_cdf(-aGoal[ii]/(2.0),nu,0.0,k)-student_t_cdf(-aGoal[ii]/(2.0)-aLBehind[ii],nu,0.0,k);
		pRBehind[ii] 	= student_t_cdf(aRBehind[ii]/+aGoal[ii]/(2.0),nu,0.0,k)-student_t_cdf(aGoal[ii]/(2.0),nu,0.0,k);
	}
	// goals only model
	target += goal .* log(pGoal)+(1.0-goal) .* log(1.0-pGoal);
	
	// likelihood contribution for goals
	//target += goal .* log(pGoal);
	// likelihood contribution for behinds
	//target += behind .* log(pLBehind+pRBehind);
	// likelihood contribution for everything else
	//target += (1.0-goal-behind) .* log(1.0-pGoal-pLBehind-pRBehind);
}

generated quantities {
	
}
