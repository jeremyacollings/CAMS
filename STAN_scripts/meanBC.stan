//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int D; // number of soil samples
  matrix[D,D] dist; // distance matrix for soil samples
  int soil_sp[D]; // species associated with soil samples
  int S; // number of species
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  matrix<lower=0, upper=1>[S,S] mean_dist; // average distance between species
  matrix<lower=0>[S,S] sd_dist; // variance in distance between species
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  for(i in 1:D){
    for(j in 1:D){
      dist[i, j] ~ normal(mean_dist[soil_sp[i], soil_sp[j]], sd_dist[soil_sp[i], soil_sp[j]]);
    }
  }
}

