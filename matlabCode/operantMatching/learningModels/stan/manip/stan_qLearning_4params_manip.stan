data {
  int<lower=1> N;   // number of pre-manipulation sessions
  int<lower=1> M;   // number of post-manipulation sessions
  int<lower=1> T;   // max number of trials across all sessions
  int<lower=1, upper=T> Tsesh[N];         // number of trials in each pre-manip session
  int<lower=0, upper=2> choice[N, T];     // choice array for pre-manip 
  int<lower=0, upper=1> outcome[N, T];    // outcome array for pre-manip
  int<lower=1, upper=T> TseshM[M];        // number of trials in each post-manip session
  int<lower=0, upper=2> choiceM[M, T];    // choice array for post-manip 
  int<lower=0, upper=1> outcomeM[M, T];   // outcome array for post-manip
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
  vector[4] d_mu_p;
  vector<lower=0>[4] d_sigma;

  // Session-level raw parameters
  vector[N] aN_pr;    // learning rate for negative prediction error (NPE)
  vector[N] aP_pr;    // learning rate for positive prediction error (PPE)
  vector[N] aF_pr;    // forgetting rate
  vector[N] beta_pr;  // inverse temperature for softmx (decision) function

  vector[M] d_aN_pr;    // change in learning rate for NPE
  vector[M] d_aP_pr;    // change in learning rate for PPE
  vector[M] d_aF_pr;    // change in forgetting rate
  vector[M] d_beta_pr;  // change in inverse temperature
}
transformed parameters {
  // session-level parameters
  vector<lower=0, upper=1>[N] aN;
  vector<lower=0, upper=1>[N] aP;
  vector<lower=0, upper=1>[N] aF;
  vector<lower=0, upper=10>[N] beta;

  vector<lower=-1, upper=1>[M] d_aN;
  vector<lower=-1, upper=1>[M] d_aP;
  vector<lower=-1, upper=1>[M] d_aF;
  vector<lower=-10, upper=10>[M] d_beta;

  for (i in 1:N) {
    aN[i]   = Phi_approx(mu_p[1]  + sigma[1]  * aN_pr[i]);
    aP[i]   = Phi_approx(mu_p[2]  + sigma[2]  * aP_pr[i]);
    aF[i]   = Phi_approx(mu_p[3]  + sigma[3]  * aF_pr[i]);
    beta[i] = Phi_approx(mu_p[4] + sigma[4] * beta_pr[i]) * 10;
    }
  for (i in 1:M){
    d_aN[i]   = -1 + 2*Phi_approx(d_mu_p[1]  + d_sigma[1]  * d_aN_pr[i]);
    d_aP[i]   = -1 + 2*Phi_approx(d_mu_p[2]  + d_sigma[2]  * d_aP_pr[i]);
    d_aF[i]   = -1 + 2*Phi_approx(d_mu_p[3]  + d_sigma[3]  * d_aF_pr[i]);
    d_beta[i] = -10 + 20*Phi_approx(d_mu_p[4]  + d_sigma[4]  * d_beta_pr[i]);
  }
}
model {
  // Hyperparameters
  mu_p  ~ normal(0, 1);
  sigma ~ cauchy(0, 1);
  d_mu_p ~ normal(0,1);
  d_sigma ~ cauchy(0,1);

  // individual parameters
  aN_pr   ~ normal(0, 1);
  aP_pr   ~ normal(0, 1);
  aF_pr   ~ normal(0, 1);
  beta_pr ~ normal(0, 1);
  d_aN_pr   ~ normal(0, 1);
  d_aP_pr   ~ normal(0, 1);
  d_aF_pr   ~ normal(0, 1);
  d_beta_pr ~ normal(0, 1);



  // session loop and trial loop
  for (i in 1:N) {
    vector[2] Q; // expected value
    real PE;      // prediction error

    Q = initQ;

    for (t in 1:(Tsesh[i])) {
      // compute action probabilities
      choice[i, t] ~ categorical_logit(beta[i] * Q);

      // prediction error
      PE = outcome[i, t] - Q[choice[i, t]];

      // value updating (learning)
      if (PE < 0){
        Q[choice[i, t]] = Q[choice[i, t]] + aN[i] * PE;
      }
      else{
        Q[choice[i, t]] = Q[choice[i, t]] + aP[i] * PE;
      }
      if (choice[i, t] == 1){
        Q[2] = Q[2] * aF[i];
      }else{
        Q[1] = Q[1] * aF[i];
      }
    }
  }
  for (i in 1:M) {
    vector[2] Q; // expected value
    real PE;      // prediction error

    Q = initQ;

    for (t in 1:(TseshM[i])) {
      // compute action probabilities
      choiceM[i, t] ~ categorical_logit((beta[i] + d_beta[i]) * Q);

      // prediction error
      PE = outcomeM[i, t] - Q[choiceM[i, t]];

      // value updating (learning)
      if (PE < 0){
        Q[choiceM[i, t]] = Q[choiceM[i, t]] + (aN[i] + d_aN[i]) * PE;
      }
      else{
        Q[choiceM[i, t]] = Q[choiceM[i, t]] + (aP[i] + d_aP[i]) * PE;
      }
      if (choiceM[i, t] == 1){
        Q[2] = Q[2] * (aF[i] + d_aF[i]);
      }else{
        Q[1] = Q[1] * (aF[i] + d_aF[i]);
      }
    }
  }
}
generated quantities {
  // For group level parameters
  real<lower=0, upper=1> mu_aN;
  real<lower=0, upper=1> mu_aP;
  real<lower=0, upper=1> mu_aF;
  real<lower=0, upper=10> mu_beta;
  real<lower=-1, upper=1> d_mu_aN;
  real<lower=-1, upper=1> d_mu_aP;
  real<lower=-1, upper=1> d_mu_aF;
  real<lower=-10, upper=10> d_mu_beta;

  // For log likelihood calculation
  real log_lik[N];
  real log_likM[M];

  // For posterior predictive check
  real y_pred[N, T];
  real y_predM[M,T];

  // Set all posterior predictions to 0 (avoids NULL values)
  for (i in 1:N) {
    for (t in 1:T) {
      y_pred[i, t] = -1;
    }
  }
  for (i in 1:M) {
    for (t in 1:T) {
      y_predM[i, t] = -1;
    }
  }

  mu_aN   = Phi_approx(mu_p[1]);
  mu_aP   = Phi_approx(mu_p[2]);
  mu_aF   = Phi_approx(mu_p[3]);
  mu_beta = Phi_approx(mu_p[4]) * 10;

  d_mu_aN   = -1 + 2*(Phi_approx(d_mu_p[1]));
  d_mu_aP   = -1 + 2*(Phi_approx(d_mu_p[2]));
  d_mu_aF   = -1 + 2*(Phi_approx(d_mu_p[3]));
  d_mu_beta = -10 + 20*(Phi_approx(d_mu_p[4]));

  { // local section, this saves time and space
    for (i in 1:N) {
      vector[2] Q; // expected value
      real PE;      // prediction error

      // Initialize values
      Q = initQ;
      log_lik[i] = 0;

      for (t in 1:(Tsesh[i])) {
        // compute log likelihood of current trial
        log_lik[i] = log_lik[i] + categorical_logit_lpmf(choice[i, t] | beta[i] * Q);

        // generate posterior prediction for current trial
        y_pred[i, t] = categorical_rng(softmax(beta[i] * Q));

        // prediction error
        PE = outcome[i, t] - Q[choice[i, t]];

        // value updating (learning)
        if (PE < 0){
          Q[choice[i, t]] = Q[choice[i, t]] + aN[i] * PE;
        }
        else{
          Q[choice[i, t]] = Q[choice[i, t]] + aP[i] * PE;
        }
        if (choice[i, t] == 1){
          Q[2] = Q[2] * aF[i];
        }else{
          Q[1] = Q[1] * aF[i];
        }
      }
    }
    for (i in 1:M) {
      vector[2] Q; // expected value
      real PE;      // prediction error

      // Initialize values
      Q = initQ;
      log_likM[i] = 0;

      for (t in 1:(TseshM[i])) {
        // compute log likelihood of current trial
        log_likM[i] = log_likM[i] + categorical_logit_lpmf(choiceM[i, t] | (beta[i] + d_beta[i]) * Q);

        // generate posterior prediction for current trial
        y_predM[i, t] = categorical_rng(softmax((beta[i] + d_beta[i]) * Q));

        // prediction error
        PE = outcomeM[i, t] - Q[choiceM[i, t]];

        // value updating (learning)
        if (PE < 0){
          Q[choiceM[i, t]] = Q[choiceM[i, t]] + (aN[i] + d_aN[i]) * PE;
        }
        else{
          Q[choiceM[i, t]] = Q[choiceM[i, t]] + (aP[i] + d_aP[i]) * PE;
        }
        if (choiceM[i, t] == 1){
          Q[2] = Q[2] * (aF[i] + d_aF[i]);
        }else{
          Q[1] = Q[1] * (aF[i] + d_aF[i]);
        }
      }
    }
  }
}
