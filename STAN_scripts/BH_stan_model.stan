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
  matrix<lower=0>[S,S] alpha; // population alphas
}

transformed parameters {
  real mu[N]; //
  matrix[S,S] ND; 
  real CR[S,S,2];
  real DR[S,S,2];
  real FR[S,S,2];
  matrix[S,S] FI;
for(n in 1:N){
  mu[n] = lambda[sp[n]]/(1 + alpha[sp[n], comp[n]]*count[n]);
  }
  for(i in 1:S){
  for(j in 1:S){
    ND[i,j] = 1 - sqrt((alpha[i,j]*alpha[j,i])/(alpha[i,i]*alpha[j,j]));
    CR[i, j, 1] = sqrt((alpha[j,j]*alpha[j, i])/(alpha[i,i]*alpha[i,j]));
    CR[i, j, 2] = sqrt((alpha[i,i]*alpha[i, j])/(alpha[j,j]*alpha[j,i]));
    DR[i, j, 1] = lambda[i]/lambda[j];
    DR[i, j, 2] = lambda[j]/lambda[i];
    FR[i,j,1] = CR[i,j,1]*DR[i,j,1];
    FR[i,j,2] = CR[i,j,2]*DR[i,j,2];
    FI[i, j] = max(FR[i,j,]);
  }
}
}

model {
  for(s in 1:S){
  lambda[s] ~ normal(1000, 250);
  sigma[s] ~ uniform(.1, 100000);
  }
  for(fs in 1:S){
    for(cs in 1:S){
      alpha[fs, cs] ~ normal(1,1);
    }
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
