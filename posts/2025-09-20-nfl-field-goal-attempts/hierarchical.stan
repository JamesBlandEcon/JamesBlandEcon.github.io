
data {
  int<lower=0> N;
  int nkickers;
  int id[N];
  
  
  vector[N] dist;
  int good[N];
  
  vector[2] prior_MU[2];
  vector[2] prior_TAU;
  real prior_Omega;
  
  
  
}

transformed data {
  
  
  real goalpost_width = 18.5/3.0; // in yards
  
  vector[N] angle = atan(goalpost_width/(2.0*dist));
  
  
}

parameters {
  vector[2] MU;
  vector<lower=0.0>[2] TAU;
  cholesky_factor_corr[2] L_Omega;
  
  matrix[2,nkickers] z;
  
}

transformed parameters {
  
  vector[nkickers] sigma;
  vector[nkickers] nu;
  
  {
    matrix[2,nkickers]  theta = rep_matrix(MU,nkickers)
                    +diag_pre_multiply(TAU,L_Omega)*z;
    sigma = exp(theta[1,]');
    nu = exp(theta[2,]');
  }
  
}

model {
  
  
  vector[N] pr_angle;
  
  for (ii in 1:N) {
    pr_angle[ii] = 2.0*student_t_cdf(angle[ii], nu[id[ii]], 0.0, sigma[id[ii]])-1.0;
  }
  
  target += bernoulli_lpmf(good | pr_angle);
  
  for (pp in 1:2) {
    target += normal_lpdf(MU[pp] | prior_MU[pp][1],prior_MU[pp][2]);
    target += cauchy_lpdf(TAU[pp] | 0.0, prior_TAU[pp]);
  }
  target += lkj_corr_cholesky_lpdf(L_Omega| prior_Omega);
  target += std_normal_lpdf(to_vector(z));
  
  
  
}

generated quantities {
  
  // predicted probability
  
  vector[nkickers] pr40;
  for (ii in 1:nkickers) {
    pr40[ii] = 2.0*student_t_cdf(atan(goalpost_width/(2.0*40.0)), nu[ii], 0.0, sigma[ii])-1.0;
  }
  
}

