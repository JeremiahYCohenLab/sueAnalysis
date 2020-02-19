data {
	int<lower=1> N; // number of neurons
  int<lower=1> T; // longest session
  int<lower=1, upper=T> Tsesh[N]; // vector of trial count for each session
  int<lower=0, upper=1> outcome[N, T]; // 0 for no reward, 1 for reward
  int<lower=0, upper=2> choice[N, T];
  real spikes[N, T]; // spike data
}
transformed data {
  vector[2] initQ;  // initial values for Q
  initQ = rep_vector(0.0, 2);
}
parameters {
  // neuron-level raw parameters
  vector[N] slope; 
  vector[N] intercept;
  vector<lower=0>[N] sigma;

  vector[N] aNscale_pr;        // learning rate for NPE
  vector[N] aNmin_pr;        // learning rate for NPE
  vector[N] aP_pr;        // learning rate for PPE
  vector[N] aF_pr;        // forgetting rate
  vector[N] v_pr;  // inverse temperature updating rate
}
transformed parameters {
  // session-level parameters
  vector<lower=0, upper=1>[N] aNscale;
  vector<lower=0, upper=1>[N] aNmin;
  vector<lower=0, upper=1>[N] aP;
  vector<lower=0, upper=1>[N] aF;
  vector<lower=0, upper=1>[N] v;

  for (n in 1:N) {
    aNscale[n]   = Phi_approx(aNscale_pr[n]);
    aNmin[n]   = Phi_approx(aNmin_pr[n]);
    aP[n]   = Phi_approx(aP_pr[n]);
    aF[n]   = Phi_approx(aF_pr[n]);
    v[n] = Phi_approx(v_pr[n]) ;
  }
}
model {
  slope     ~ normal(0,10);
  intercept ~ normal(0,10);
  sigma     ~ cauchy(0,5);

  aNscale_pr  ~ normal(0, 1);
  aNmin_pr    ~ normal(0, 1);
  aP_pr       ~ normal(0, 1);
  aF_pr       ~ normal(0, 1);
  v_pr        ~ normal(0, 1);

  // neuron loop and trial loop
  for (n in 1:N) {
    vector[2] Q;               // expected value
    vector[Tsesh[n]] Qdiff;
    real PE;                   // prediction error
    real aN;
    vector[Tsesh[n]] V;        // volatility history term
    
    Q = initQ;
    V[1] = 0.5;                // hard-coded start value
    aN = aNmin[n];


    for (t in 1:(Tsesh[n] - 1)) {
      Qdiff[t] = Q[2] - Q[1];

      if (choice[n,t] == 1) {
        PE = outcome[n, t] - Q[2];
        if (PE < 0){
          Q[2] += aN * PE;
        }else{
          Q[2] += aP[n] * PE;
        }
        Q[1] = Q[1] * aF[n];
      }else{
        PE = outcome[n, t] - Q[1];
        if (PE < 0){
          Q[1] += aN * PE;
        }else{
          Q[1] += aP[n] * PE;
        }
        Q[2] = Q[2] * aF[n];
      }
      V[t+1] = v[n] * abs(PE) + (1-v[n]) * V[t];
      aN = V[t+1] * aNscale[n] + aNmin[n];
    }
    spikes[n, 1:Tsesh[n]] ~ normal(slope[n]*V + intercept[n], sigma[n]); 
  }
}

generated quantities{
  matrix[N,T] V;            // reward history term
  matrix[N,T] spikes_pred;

  for (n in 1:N) {
    vector[2] Q;               // expected value
    vector[Tsesh[n]] Qdiff;
    real PE;                   // prediction error
    real aN;
    
    Q = initQ;
    V[n, 1] = 0.5;                // hard-coded start value
    aN = aNmin[n];


    for (t in 1:(Tsesh[n] - 1)) {
      Qdiff[t] = Q[2] - Q[1];

      if (choice[n,t] == 1) {
        PE = outcome[n, t] - Q[2];
        if (PE < 0){
          Q[2] += aN * PE;
        }else{
          Q[2] += aP[n] * PE;
        }
        Q[1] = Q[1] * aF[n];
      }else{
        PE = outcome[n, t] - Q[1];
        if (PE < 0){
          Q[1] += aN * PE;
        }else{
          Q[1] += aP[n] * PE;
        }
        Q[2] = Q[2] * aF[n];
      }
      V[n, t+1] = v[n] * abs(PE) + (1-v[n]) * V[n,t];
      aN = V[n, t+1] * aNscale[n] + aNmin[n];
    }
    spikes_pred[n, 1:Tsesh[n]] = slope[n]*V[n,1:Tsesh[n]] + intercept[n]; 
  }
}
