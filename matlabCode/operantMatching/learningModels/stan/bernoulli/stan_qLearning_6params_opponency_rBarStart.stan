data {
  int<lower=1> N;
  int<lower=1> T;
  int<lower=1, upper=T> Tsesh[N];
  int<lower=0, upper=2> choice[N, T];
  int<lower=0, upper=1> outcome[N, T];
}
transformed data {
  vector[2] initQ;  // initial values for Q
  initQ = rep_vector(0.0, 2);
}
parameters {
// Declare all parameters as vectors for vectorizing
  // Hyper(animal)-parameters
  vector[6] mu_p;
  vector<lower=0>[6] sigma;

  // Session-level raw parameters
  vector[N] aN_pr;    // learning rate for NPE
  vector[N] aP_pr;    // learning rate for PPE
  vector[N] aF_pr;    // forgetting rate
  vector[N] beta_pr;  // inverse temperature
  vector[N] v_pr;     // expected average value learning rate
  vector[N] rBarStart_pr;     // expected average value start value 
}
transformed parameters {
  // session-level parameters
  vector<lower=0, upper=1>[N] aN;
  vector<lower=0, upper=1>[N] aP;
  vector<lower=0, upper=1>[N] aF;
  vector<lower=0, upper=10>[N] beta;
  vector<lower=0, upper=1>[N] v;
  vector<lower=0, upper=1>[N] rBarStart;

  for (n in 1:N) {
    aN[n]        = Phi_approx(mu_p[1]  + sigma[1]  * aN_pr[n]);
    aP[n]        = Phi_approx(mu_p[2]  + sigma[2]  * aP_pr[n]);
    aF[n]        = Phi_approx(mu_p[3]  + sigma[3]  * aF_pr[n]);
    beta[n]      = Phi_approx(mu_p[4] + sigma[4] * beta_pr[n]) * 10;
    v[n]         = Phi_approx(mu_p[5]  + sigma[5]  * v_pr[n]);
    rBarStart[n] = Phi_approx(mu_p[6]  + sigma[6]  * v_pr[n]);
  }
}
model {
  // Hyperparameters
  mu_p  ~ normal(0, 1);
  sigma ~ cauchy(0, 1);

  // individual parameters
  aN_pr        ~ normal(0, 1);
  aP_pr        ~ normal(0, 1);
  aF_pr        ~ normal(0, 1);
  beta_pr      ~ normal(0, 1);
  v_pr         ~ normal(0, 1);
  rBarStart_pr ~ normal(0, 1);

  // session loop and trial loop
  for (n in 1:N) {
    vector[2] Q; // expected value
    real PE;      // prediction error
    real rBar; // expected average value
    vector[Tsesh[n]] Qdiff;

    Q = initQ;
    rBar = rBarStart[n];

    for (t in 1:(Tsesh[n])) {
      Qdiff[t] = Q[2] - Q[1];
      choice[n, t] ~ bernoulli_logit(beta[n] * Qdiff[t]);

      if (choice[n,t] == 1) {
        PE = outcome[n, t] - Q[2] - rBar;
        if (PE < 0){
          Q[2] += aN[n] * PE;
        }else{
          Q[2] += aP[n] * PE;
        }
        Q[1] = Q[1] * aF[n];
      }else{
        PE = outcome[n, t] - Q[1] - rBar;
        if (PE < 0){
          Q[1] += aN[n] * PE;
        }else{
          Q[1] += aP[n] * PE;
        }
        Q[2] = Q[2] * aF[n];
      }
      rBar = v[n] * outcome[n, t] + (1-v[n]) * rBar;
    }
  }
}
generated quantities {
  // For group level parameters
  real<lower=0, upper=1> mu_aN;
  real<lower=0, upper=1> mu_aP;
  real<lower=0, upper=1> mu_aF;
  real<lower=0, upper=10> mu_beta;
  real<lower=0, upper=1> mu_v;
  real<lower=0, upper=1> mu_rBarStart;

  // For log likelihood calculation
  real log_lik[N];

  mu_aN        = Phi_approx(mu_p[1]);
  mu_aP        = Phi_approx(mu_p[2]);
  mu_aF        = Phi_approx(mu_p[3]);
  mu_beta      = Phi_approx(mu_p[4]) * 20;
  mu_v         = Phi_approx(mu_p[5]);
  mu_rBarStart = Phi_approx(mu_p[6]);

{ // local section, this saves time and space
    for (n in 1:N) {
      vector[2] Q; // expected value
      real PE;      // prediction error
      real rBar; // expected average value
      vector[Tsesh[n]] Qdiff;

      // Initialize values
      Q = initQ;
      rBar = rBarStart[n];

      log_lik[n] = 0;

      for (t in 1:(Tsesh[n])) {
        Qdiff[t] = Q[2] - Q[1];

        // compute log likelihood of current trial
        log_lik[n] = log_lik[n] + bernoulli_logit_lpmf(choice[n, t] | beta[n] * Qdiff[t]);

        if (choice[n,t] == 1) {
          PE = outcome[n, t] - Q[2] - rBar;
          if (PE < 0){
            Q[2] += aN[n] * PE;
          }else{
            Q[2] += aP[n] * PE;
          }
          Q[1] = Q[1] * aF[n];
        }else{
          PE = outcome[n, t] - Q[1] - rBar;
          if (PE < 0){
            Q[1] += aN[n] * PE;
          }else{
            Q[1] += aP[n] * PE;
          }
          Q[1] = Q[1] * aF[n];
        }
        rBar = v[n] * outcome[n, t] + (1-v[n]) * rBar;
      }
    }
  }
}
