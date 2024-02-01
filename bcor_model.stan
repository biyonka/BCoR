//
data{
	  int<lower=0> K;
	  int<lower=0> N;
    int<lower=0> d;
    int<lower=0> time_horizon;
    int<lower=0> t_c; //current t
	  matrix[N,K] X;
	  matrix[time_horizon,d] B;
	  int<lower=0, upper=1> y[N, t_c];
	  int<lower=0> S[N, t_c];
    int<lower=0> A[N, t_c]; 
}

// The parameters accepted by the model.
parameters {
  vector[d] eta[2,2]; 
  real<lower=0> tau[2, 2]; //tau is a 2x2 matrix where the ijth entry is the variance for alpha_ij
	real b0;
	real b1;
	vector[K] mu_beta;
	vector[K] beta_noise[2, 2];
	vector[N] alpha_prime[2, 2]; //alpha is a 2x2 array of vectors of length N representing the random effects of each arm for each s,a
}
transformed parameters {
  vector[N] alpha[2, 2];
  vector[K] beta[2,2]; 
  alpha[1,1] = alpha_prime[1,1]*sqrt(tau[1, 1]);
  alpha[2,1] = alpha_prime[2,1]*sqrt(tau[2, 1]);
  alpha[1,2] = alpha_prime[1,2]*sqrt(tau[1,2]) + (b0*(alpha_prime[1,1]*sqrt(tau[1,1])) + b1*(alpha_prime[2,1]*sqrt(tau[2,1])));
  alpha[2,2] = alpha_prime[2,2]*sqrt(tau[2,2]) + (b0*(alpha_prime[1,1]*sqrt(tau[1,1])) + b1*(alpha_prime[2,1]*sqrt(tau[2,1])));
  for (s in 1:2){
    beta[s, 1]=mu_beta+beta_noise[s, 1]; 
    beta[s, 2]=mu_beta+beta_noise[s, 2];
  }
}
// The model to be estimated.
model {
  b0~normal(0, 0.1);
  b1~normal(0, 0.1);
  mu_beta~normal(0, 0.3);
  for (s in 1:2){
    tau[s, 1]~inv_gamma(100, 1); // shape (tau), rate(sigma) parameterization of inv-gamma
    tau[s, 2]~inv_gamma(100, 1);
    eta[s, 1]~normal(0, 0.3);
    eta[s, 2]~normal(0, 0.3);
    alpha_prime[s, 1] ~ normal(0, 1);
    alpha_prime[s, 2] ~ normal(0, 1);
    beta_noise[s, 1] ~ normal(0, 0.1);
    beta_noise[s, 2] ~ normal(0, 0.1);
  }
	for(n in 1:N)
	  for(i in 1:t_c){
	    y[n, i] ~ bernoulli(Phi(alpha[S[n, i], A[n, i], n] + X[n]*beta[S[n, i], A[n, i]] + B[i]*eta[S[n, i], A[n, i]])); 
	}
}
