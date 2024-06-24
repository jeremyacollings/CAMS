data {
  int N; // number of plants
  int S; // number of species
  int sp[N]; // species of focal plants
  int comp[N]; // species of competitor plants
  int count[N]; // number of competitor plants
  real fec[N]; // fecundity of focal plants
}

parameters { 
  vector<lower=0>[S] lambda; // population lambda per sp
  vector<lower=0>[S] sigma; // sd of fec per sp
}

transformed parameters {
  real mu[N]; //
  matrix[S,S] ND; 
  real CR[S,S,2];
  real DR[S,S,2];
  real FR[S,S,2];
  matrix[S,S] FI;
for(n in 1:N){
  mu[n] = lambda[sp[n]];
  }
}

model {
  for(s in 1:S){
  lambda[s] ~ normal(1000, 250);
  sigma[s] ~ uniform(.1, 100000);
  }
  for(n in 1:N){
  fec[n] ~ normal(mu[n], sigma[sp[n]]);
  }
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = normal_lpdf(fec[n] | mu[n], sigma[sp[n]]);
  }
}
