
data {
  int<lower=0> N; // number of observations
  int<lower=0> n; // number of participants
  vector[N] chi;
  vector[N] Part;
  vector[N] UseData;
  int<lower=1> ID[N];
  real priorR[4];
  real priorL[2];
  
  int nInvestGrid;
  int InvestGrid[nInvestGrid];
  int investL[N];
  int investR[N];
  real investFrac[N];
}

transformed data {
  
}

//Define the model's parameters
parameters {
  real mu_r;
  real<lower=0> sigma_r;
  vector<lower=0>[n] r;
  real<lower=0> lambda;
  
  
}

transformed parameters {
  
  
}


model {
  
  vector[N] log_lik;
  {
    for (ii in 1:N) {
      vector[nInvestGrid] U;
      vector[nInvestGrid] lpY;
      for (uu in 1:nInvestGrid) {
        U[uu] = 0.5*(100-InvestGrid[uu]+InvestGrid[uu] * chi[ii])^r[ID[ii]]+0.5*(100.0-InvestGrid[uu])^r[ID[ii]];
      }
      lpY=log_softmax(lambda*U/(max(U)-min(U)));
      log_lik[ii] = investFrac[ii]*lpY[investL[ii]]+(1.0-investFrac[ii])*lpY[investR[ii]];
    }
  }
  
  target+= log_lik;
  
  mu_r ~ normal(priorR[1],priorR[2]);
  sigma_r~ lognormal(priorR[3],priorR[4]);
  
  r ~ lognormal(mu_r,sigma_r);
  lambda ~ lognormal(priorL[1],priorL[2]);
  
 
}

