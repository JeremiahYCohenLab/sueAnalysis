data {
	int<lower=1> N; // number of neurons
  int<lower=1> T; // longest session
  int<lower=1, upper=T> Tsesh[N]; // vector of trial count for each session
  int<lower=0, upper=1> outcome[N, T]; // 0 for no reward, 1 for reward
  real spikes[N, T]; // spike data
}
parameters {
  // neuron-level raw parameters
  vector[N] R_pr;
  vector[N] aPE_pr;
  vector[N] slope; 
  vector[N] intercept;
  vector<lower=0>[N] sigma;
}
transformed parameters{
  // transform from hierarchical level to neuron-level
 vector<lower=0, upper=1>[N] R;
 vector<lower=0, upper=1>[N] aPE;

  //this will do a normal cdf transformation from [-inf,inf] to [0,1]
  for (n in 1:N) {
    R[n]   = Phi_approx(R_pr[n]);
    aPE[n] = Phi_approx(aPE_pr[n]);
  }
}
model {
  R_pr      ~ normal(0,1);
  aPE_pr    ~ normal(0,1);
  slope     ~ normal(0,10);
  intercept ~ normal(0,10);
  sigma     ~ cauchy(0,5);

  // neuron loop and trial loop
  for (n in 1:N) {
    real pe;                // prediction error
    vector[Tsesh[n]] peBar; // expected average value

    peBar[1] = 0;

    for (t in 1:(Tsesh[n] - 1)) {
      pe = outcome[n, t] - R[n];
      peBar[t+1] = peBar[t] + aPE[n] * (peBar[t] - pe);
    }
    spikes[n, 1:Tsesh[n]] ~ normal(slope[n]*peBar + intercept[n], sigma[n]); 
  }
}

generated quantities{
  matrix[N,T] spikes_pred;
  real pe;
  matrix[N,T] peBar;

  for (n in 1:N){
    peBar[n,1] = 0;

    for (t in 1:(Tsesh[n] - 1)) {
      pe = outcome[n, t] - R[n];
      peBar[n, t+1] = peBar[n, t] + aPE[n] * (peBar[n, t] - pe);
    }

    spikes_pred[n, 1:Tsesh[n]] = slope[n]*peBar[n,1:Tsesh[n]] + intercept[n];       //might want to convert this to median and std before matlab to save memory/time
  }
}
