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
  vector[N] v2_pr;
  vector[N] slope; 
  vector[N] intercept;
  vector<lower=0>[N] sigma;
}
transformed parameters{
  // transform from hierarchical level to neuron-level
 vector<lower=0, upper=1>[N] v;
 vector<lower=0, upper=1>[N] v2;

  //this will do a normal cdf transformation from [-inf,inf] to [0,1]
  for (n in 1:N) {
    v[n] = Phi_approx(v_pr[n])/10;
    v2[n] = Phi_approx(v2_pr[n]);
  }
}
model {
  v_pr      ~ normal(0,1);
  v2_pr      ~ normal(0,1);
  slope     ~ normal(0,10);
  intercept ~ normal(0,10);
  sigma     ~ cauchy(0,5);

  // neuron loop and trial loop
  for (n in 1:N) {
    vector[Tsesh[n]] R;       // reward history term
    vector[Tsesh[n]] R2;       // reward history term
    R[1] = 0;                 // hard-coded start value
    R2[1] = 0;                 // hard-coded start value

    for (t in 1:(Tsesh[n] - 1)) {
      R[t+1] = v[n] * outcome[n, t] + (1-v[n]) * R[t];
      R2[t+1] = v2[n] * outcome[n, t] + (1-v2[n]) * R2[t];
    }
    spikes[n, 1:Tsesh[n]] ~ normal(slope[n]*(R+R2) + intercept[n], sigma[n]); 
  }
}

generated quantities{
  matrix[N,T] R;            // reward history term
  matrix[N,T] R2;            // reward history term
  matrix[N,T] spikes_pred;

  for (n in 1:N){
    R[n,1] = 0;               // hard-coded start value
    R2[n,1] = 0;               // hard-coded start value

    for (t in 1:(Tsesh[n] - 1)) {
      R[n, t+1] = v[n] * outcome[n, t] + (1-v[n]) * R[n,t];
      R2[n, t+1] = v2[n] * outcome[n, t] + (1-v2[n]) * R2[n,t];
    }

    spikes_pred[n, 1:Tsesh[n]] = slope[n]*(R[n,1:Tsesh[n]] + R2[n,1:Tsesh[n]]) + intercept[n];       //might want to convert this to median and std before matlab to save memory/time
  }
}
