data {
	int<lower=1> N; // number of neurons
  int<lower=1> T; // longest session
  int<lower=1, upper=T> Tsesh[N]; // vector of trial count for each session
  int<lower=0, upper=1> outcome[N, T]; // 0 for no reward, 1 for reward
  real spikes[N, T]; // spike data
}
parameters {
  // neuron-level raw parameters
  vector[N] slope; 
  vector[N] intercept;
  vector<lower=0>[N] sigma;
}
model {
  slope     ~ normal(0,10);
  intercept ~ normal(0,10);
  sigma     ~ cauchy(0,5);

  // neuron loop and trial loop
  for (n in 1:N) {
    vector[Tsesh[n]] V;       // reward history term
    V[1] = 0.5;                 // hard-coded start value
    V[2] = 0.5;

    for (t in 2:(Tsesh[n] - 1)) {
      if(t > 10){
        V[t+1] = sd(outcome[n,t-10:t]);
      }else{
        V[t+1] = sd(outcome[n,1:t]);
      }
    }
    spikes[n, 1:Tsesh[n]] ~ normal(slope[n]*V + intercept[n], sigma[n]); 
  }
}

generated quantities{
  matrix[N,T] outcome;            // reward history term
  matrix[N,T] spikes_pred;

  for (n in 1:N){
    V[n,1] = 0.5;                 // hard-coded start value
    V[n,2] = 0.5;

    for (t in 2:(Tsesh[n] - 1)) {
      if(t > nT){
        V[t+1] = sd(outcome[n,t-10:t]);
      }else{
        V[t+1] = sd(outcome[n,1:t]);
      }
    }
    spikes_pred[n, 1:Tsesh[n]] = slope[n]*outcome[n,1:Tsesh[n]] + intercept[n];       //might want to convert this to median and std before matlab to save memory/time
  }
}
