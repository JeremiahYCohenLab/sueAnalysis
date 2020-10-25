data {
  int<lower=1> N;
  int<lower=1> T;
  int<lower=1, upper=T> Tsesh[N];
  int<lower=0, upper=2> choice[N, T];
  int<lower=0, upper=1> outcome[N, T];
  real ITI[N, T]; //no lower or upper bounds
}
transformed data {
  vector[2] initQ;  // initial values for Q
  initQ = rep_vector(0.0, 2);
}
parameters {
// Declare all parameters as vectors for vectorizing
  // Hyper(animal)-parameters
  vector[4] mu_p;
  vector<lower=0>[4] sigma;

  // Session-level raw parameters
  vector[N] aN_pr;    // learning rate for NPE
  vector[N] aP_pr;    // learning rate for PPE
  vector[N] tF_pr;    // forgetting rate
  vector[N] beta_pr;  // inverse temperature
}
transformed parameters {
// Transform session-level raw parameters
  vector<lower=0, upper=1>[N] aN;
  vector<lower=0, upper=1>[N] aP;
  vector<lower=0, upper=1>[N] tF;
  vector<lower=0, upper=40>[N] beta;

  for (n in 1:N) {
    aN[n]   = Phi_approx(mu_p[1] + sigma[1] * aN_pr[n]);
    aP[n]   = Phi_approx(mu_p[2] + sigma[2] * aP_pr[n]);
    tF[n]   = Phi_approx(mu_p[3] + sigma[3] * tF_pr[n]);
    beta[n] = Phi_approx(mu_p[4] + sigma[4] * beta_pr[n]) * 40;
  }
}
model {
  // Hyperparameters
  mu_p  ~ normal(0, 1);
  sigma ~ cauchy(0, 1);

  // individual parameters
  aN_pr   ~ normal(0, 1);
  aP_pr   ~ normal(0, 1);
  tF_pr   ~ normal(0, 1);
  beta_pr ~ normal(0, 1);

  // session loop and trial loop
  for (n in 1:N) {
    vector[2] Q; // expected value
    real PE;      // prediction error
    vector[Tsesh[n]] Qdiff;

    Q = initQ;

    for (t in 1:(Tsesh[n])) {
      Q = Q*(tF[n]^ITI[n,t]);
      //Q = Q*exp((-1)*tF[n]*ITI[n,t]);
      Qdiff[t] = Q[2] - Q[1];
      choice[n, t] ~ bernoulli_logit(beta[n] * Qdiff[t]);

      if (choice[n,t] == 1) {
        PE = outcome[n, t] - Q[2];
        if (PE < 0){
          Q[2] = Q[2] + aN[n] * PE;
        }else{
          Q[2] = Q[2] + aP[n] * PE;
        }
      }else{
        PE = outcome[n, t] - Q[1];
        if (PE < 0){
          Q[1] = Q[1] + aN[n] * PE;
        }else{
          Q[1] = Q[1] + aP[n] * PE;
        }
      }
    }
  }
}
generated quantities {
  // For group level parameters
  real<lower=0, upper=1> mu_aN;
  real<lower=0, upper=1> mu_aP;
  real<lower=0, upper=1> mu_tF;
  real<lower=0, upper=40> mu_beta;

  // For log likelihood calculation
  real log_lik[N];

  // For posterior predictive check
  real y_pred[N, T];

  // Set all posterior predictions to 0 (avoids NULL values)
  for (n in 1:N) {
    for (t in 1:T) {
      y_pred[n, t] = -1;
    }
  }

  mu_aN   = Phi_approx(mu_p[1]);
  mu_aP   = Phi_approx(mu_p[2]);
  mu_tF   = Phi_approx(mu_p[3]);
  mu_beta = Phi_approx(mu_p[4])*40;

  { // local section, this saves time and space
    for (n in 1:N) {
      vector[2] Q; // expected value
      real PE;      // prediction error
      vector[Tsesh[n]] Qdiff;

      // Initialize values
      Q = initQ;

      log_lik[n] = 0;

      for (t in 1:(Tsesh[n])) {
        Q = Q*(tF[n]^ITI[n,t]);
        //Q = Q*exp((-1)*tF[n]*ITI[n,t]);
        Qdiff[t] = Q[2] - Q[1];
        // compute log likelihood of current trial
        log_lik[n] = log_lik[n] + bernoulli_logit_lpmf(choice[n, t] | beta[n] * Qdiff[t]);

        // generate posterior prediction for current trial
        y_pred[n, t] = categorical_rng(softmax(beta[n] * Q));

        if (choice[n,t] == 1) {
          PE = outcome[n, t] - Q[2];
          if (PE < 0){
            Q[2] = Q[2] + aN[n] * PE;
          }else{
            Q[2] = Q[2] + aP[n] * PE;
          }
        }else{
          PE = outcome[n, t] - Q[1];
          if (PE < 0){
            Q[1] = Q[1] + aN[n] * PE;
          }else{
            Q[1] = Q[1] + aP[n] * PE;
          }
        }
      }
    }
  }
}
