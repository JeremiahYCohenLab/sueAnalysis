data {
	int<lower=1> N; // number of neurons
  int<lower=1> T; // longest session
  int<lower=1, upper=T> Tsesh[N]; // vector of trial count for each session
  int<lower=0, upper=1> outcome[N, T]; // 0 for no reward, 1 for reward
  real spikes[N, T]; // spike data
}
parameters {
  // neuron-level raw parameters
  vector[N] v_pr;
  vector[N] Rmin; 
  vector[N] slope; 
  vector[N] intercept;
  vector<lower=0>[N] sigma;
}
transformed parameters{
  // transform from hierarchical level to neuron-level
 vector<lower=0, upper=1>[N] v;
 vector<lower=0, upper=1>[N] Rmin;
 vector<lower=0, upper=1>[N] aPE;

  //this will do a normal cdf transformation from [-inf,inf] to [0,1]
  for (n in 1:N) {
    v[n] = Phi_approx(v_pr[n]);
    Rmin[n] = Phi_approx(Rmin_pr[n]);
    aPE[n] = Phi_approx(aPE_pr[n]);
  }
}
model {
  v_pr      ~ normal(0,1);
  Rmin_pr   ~ normal(0,1);
  aPE_pr    ~ normal(0,1);
  slope     ~ normal(0,10);
  intercept ~ normal(0,10);
  sigma     ~ cauchy(0,5);

  // neuron loop and trial loop
  for (n in 1:N) {
    vector[Tsesh[n]] R;     // reward history term
    real pe;                // prediction error
    vector[Tsesh[n]] peBar; // expected average value

    peBar[1] = 0;
    R[1] = Rmin[n];                 // hard-coded start value

    for (t in 1:(Tsesh[n] - 1)) {
      pe = outcome[n, t] - R[t];
      peBar[t+1] = peBar[t] + aPE[n] * (peBar[t] - pe);
      R[t+1] = v[n] * outcome[n, t] + (1-v[n]) * R[t];
    }
    spikes[n, 1:Tsesh[n]] ~ normal(slope[n]*peBar + intercept[n], sigma[n]); 
  }
}

generated quantities{
  matrix[N,T] R;            // reward history term
  matrix[N,T] spikes_pred;
  real pe;
  matrix[N,T] peBar;

  for (n in 1:N){
    peBar[n,1] = 0;
    R[n,1] = 0;               // hard-coded start value

    for (t in 1:(Tsesh[n] - 1)) {
      pe = outcome[n, t] - R[t];
      peBar[n, t+1] = peBar[n, t] + aPE[n] * (peBar[n, t] - pe);
      R[n, t+1] = v[n] * outcome[n, t] + (1-v[n]) * R[n,t];
    }

    spikes_pred[n, 1:Tsesh[n]] = slope[n]*peBar[n,1:Tsesh[n]] + intercept[n];       //might want to convert this to median and std before matlab to save memory/time
  }
}
