data {
  int<lower=1> N;
  int<lower=1> T;
  int<lower=1> S;
  int s[S+1];
  int<lower=1, upper=T> Tsesh[N];
  int<lower=0, upper=2> choice[N, T];
  int<lower=0, upper=1> laser[N,T];
  real outcome[N, T];  // no lower and upper bounds
}
transformed data {
  vector[2] initQ;  // initial values for Q
  initQ = rep_vector(0.0, 2);
}
parameters {
// Declare all parameters as vectors for vectorizing
  // Hyper(animal)-parameters
  real diff_mu;
  // vector<lower=0>[S, 5] sigma;
  real<lower=0> sigma_diff;
  // Session-level raw parameters
  vector[S] a_pr;    // learning rate for RPE
  vector[S] aF_pr;    // forgetting rate
  vector[S] beta_pr;  // inverse temperature
  vector[S] k_pr;        // choice autocorrelation
  vector[S] diff_pr;   // rpe change
  vector[N] bias;    //choice bias for each session
}
transformed parameters {
// Transform session-level raw parameters
  vector<lower=0, upper=1>[S] a;
  vector<lower=0, upper=1>[S] aF;
  vector<lower=0, upper=20>[S] beta;
  vector<lower=0, upper=3>[S] k;
  vector<lower=-1, upper=1>[S] diff;

  for (ani in 1:S) {
    a[ani]   = Phi_approx(a_pr[ani]);
    aF[ani]   = Phi_approx(aF_pr[ani]);
    beta[ani] = Phi_approx(beta_pr[ani]) * 20;
    k[ani]    = Phi_approx(k_pr[ani]) * 3;
    diff[ani] = Phi_approx(diff_mu + sigma_diff * diff_pr[ani]) * 2 - 1;
  }
}
model {
  // Hyperparameter
  diff_mu ~ normal(0, 1);
  sigma_diff ~ cauchy(0, 1);

  // individual parameters
  a_pr   ~ normal(0, 1);
  aF_pr   ~ normal(0, 1);
  beta_pr ~ normal(0, 1);
  k_pr    ~ normal(0, 1);
  diff_pr ~ normal(0, 1);
  bias   ~ normal(0, 20);

  // session loop and trial loop
  for (ani in 1:S) {
    for (n in s[ani]:(s[ani+1]-1)) {
      vector[2] Q; // expected value
      real prevChoice;
      real PE;      // prediction error
      vector[Tsesh[n]] Qdiff;

      Q = initQ;
      prevChoice = 0;

      for (t in 1:(Tsesh[n])) {
        Qdiff[t] = Q[2] - Q[1];
        choice[n, t] ~ bernoulli_logit(beta[ani] * Qdiff[t] + k[ani] * prevChoice + bias[n]);

        if (choice[n,t] == 1) {
          PE = outcome[n, t] - Q[2];
          if(laser[n, t] == 1){
            PE = PE + diff[ani];
          }
          Q[2] = Q[2] + a[ani] * PE;
          Q[1] = Q[1] * aF[ani];
          prevChoice = 1;
        }else{
          PE = outcome[n, t] - Q[1];
          if(laser[n, t] == 1){
            PE = PE + diff[ani];
          }
          Q[1] = Q[1] + a[ani] * PE;
          Q[2] = Q[2] * aF[ani];
          prevChoice = -1;
        }
      }
    }
  }
}

generated quantities {
  // For group level parameters
  real<lower=-1, upper=1> mu_diff;

  // For log likelihood calculation
  real log_lik[N];
  real log_likMean[N];


  mu_diff = 2*Phi_approx(diff_mu)-1;

  {
    for (ani in 1:S) {
      for (n in s[ani]:(s[ani+1]-1)) {
        vector[2] Q; // expected value
        real prevChoice;
        real PE;      // prediction error
        vector[Tsesh[n]] Qdiff;

        Q = initQ;
        prevChoice = 0;
        log_lik[n] = 0;

        for (t in 1:(Tsesh[n])) {
          Qdiff[t] = Q[2] - Q[1];

          // compute log likelihood of current trial
          log_lik[n] = log_lik[n] + bernoulli_logit_lpmf(choice[n, t] | beta[ani] * Qdiff[t] + k[ani] * prevChoice + bias[n]);

          if (choice[n,t] == 1) {
            PE = outcome[n, t] - Q[2];
            if(laser[n, t] == 1){
              PE = PE + diff[ani];
            }
            Q[2] = Q[2] + a[ani] * PE;
            Q[1] = Q[1] * aF[ani];
            prevChoice = 1;
          }else{
            PE = outcome[n, t] - Q[1];
            if(laser[n, t] == 1){
              PE = PE + diff[ani];
            }
            Q[1] = Q[1] + a[ani] * PE;
            Q[2] = Q[2] * aF[ani];
            prevChoice = -1;
          }
        }
        log_likMean[n] = log_lik[n]/Tsesh[n];
      }
    }
  }
}
