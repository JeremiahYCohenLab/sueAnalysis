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
  // Hyper-hyperparameters
  vector[2] h_mu_p;
  vector<lower=0>[2] h_sigma;

  // Hyperparameters
  vector[5] mu_p_pr;
  vector<lower=0>[5] sigma_pr;
  vector[5] d_mu_p_pr;
  vector<lower=0>[5] d_sigma_pr;

  // Session-level raw parameters
  vector[N] aN_pr;    // learning rate for negative prediction error (NPE)
  vector[N] aP_pr;    // learning rate for positive prediction error (PPE)
  vector[N] aF_pr;    // forgetting rate
  vector[N] v_pr;    // average reward update rate
  vector[N] betaScale_pr;  // inverse temperature for softmax (decision) function

  vector[M] d_aN_pr;    // change in learning rate for NPE
  vector[M] d_aP_pr;    // change in learning rate for PPE
  vector[M] d_aF_pr;    // change in forgetting rate
  vector[M] d_v_pr;    // change in average reward update rate
  vector[M] d_betaScale_pr;  // change in inverse temperature
}
transformed parameters {
  //hyperparameters
  vector[5] mu_p;
  vector[5] sigma;
  vector[5] d_mu_p;
  vector[5] d_sigma;

  // session-level parameters
  vector<lower=0, upper=1>[N] aN;
  vector<lower=0, upper=1>[N] aP;
  vector<lower=0, upper=1>[N] aF;
  vector<lower=0, upper=1>[N] v;
  vector<lower=0, upper=20>[N] betaScale;

  vector<lower=0, upper=1>[M] d_aN;
  vector<lower=0, upper=1>[M] d_aP;
  vector<lower=0, upper=1>[M] d_aF;
  vector<lower=0, upper=1>[M] d_v;
  vector<lower=0, upper=20>[M] d_betaScale;

  for (i in 1:5){
    mu_p[i] = h_mu_p[1]  + h_sigma[1]  * mu_p_pr[i];
    sigma[i] = h_mu_p[2]  + h_sigma[2]  * sigma_pr[i];
    d_mu_p[i] = h_mu_p[1]  + h_sigma[1]  * d_mu_p_pr[i];
    d_sigma[i] = h_mu_p[2]  + h_sigma[2]  * d_sigma_pr[i];
  }

  for (i in 1:N) {
    aN[i]   = Phi_approx(mu_p[1]  + sigma[1]  * aN_pr[i]);
    aP[i]   = Phi_approx(mu_p[2]  + sigma[2]  * aP_pr[i]);
    aF[i]   = Phi_approx(mu_p[3]  + sigma[3]  * aF_pr[i]);
    v[i]   = Phi_approx(mu_p[4]  + sigma[4]  * v_pr[i]);
    betaScale[i] = 20*Phi_approx(mu_p[5] + sigma[5] * betaScale_pr[i]);
    }
  for (i in 1:M){
    d_aN[i]   = Phi_approx(d_mu_p[1]  + d_sigma[1]  * d_aN_pr[i]);
    d_aP[i]   = Phi_approx(d_mu_p[2]  + d_sigma[2]  * d_aP_pr[i]);
    d_aF[i]   = Phi_approx(d_mu_p[3]  + d_sigma[3]  * d_aF_pr[i]);
    d_v[i]   = Phi_approx(d_mu_p[4]  + d_sigma[4]  * d_v_pr[i]);
    d_betaScale[i] = 20*Phi_approx(d_mu_p[5] + d_sigma[5] * d_betaScale_pr[i]);
  }
}
model {
  //hyper-hyperparameters
  h_mu_p  ~ normal(0, 0.5);
  h_sigma ~ cauchy(0, 0.5);

  // Hyperparameters
  mu_p_pr  ~ normal(0, 1);
  sigma_pr ~ cauchy(0, 1);
  d_mu_p_pr  ~ normal(0, 1);
  d_sigma_pr ~ cauchy(0, 1);

  // individual parameters
  aN_pr   ~ normal(0, 1);
  aP_pr   ~ normal(0, 1);
  aF_pr   ~ normal(0, 1);
  v_pr   ~ normal(0, 1);
  betaScale_pr ~ normal(0, 1);
  d_aN_pr   ~ normal(0, 1);
  d_aP_pr   ~ normal(0, 1);
  d_aF_pr   ~ normal(0, 1);
  d_v_pr   ~ normal(0, 1);
  d_betaScale_pr ~ normal(0, 1);



  // session loop and trial loop
  for (i in 1:N) {
    vector[2] Q; // expected value
    real PE;      // prediction error
    real R;
    real beta;

    Q = initQ;
    R = 0;
    beta = 0;

    for (t in 1:(Tsesh[i])) {
      // compute action probabilities
      choice[i, t] ~ categorical_logit(beta * Q);

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
      R = v[i] * outcome[i, t] + (1-v[i]) * R;
      beta = R * betaScale[i];
    }
  }
  for (i in 1:M) {
    vector[2] Q; // expected value
    real PE;      // prediction error
    real R;
    real beta;

    Q = initQ;
    R = 0;
    beta = 0;

    for (t in 1:(TseshM[i])) {
      // compute action probabilities
      choiceM[i, t] ~ categorical_logit(beta * Q);

      // prediction error
      PE = outcomeM[i, t] - Q[choiceM[i, t]];

      // value updating (learning)
      if (PE < 0){
        Q[choiceM[i, t]] = Q[choiceM[i, t]] + d_aN[i] * PE;
      }
      else{
        Q[choiceM[i, t]] = Q[choiceM[i, t]] + d_aP[i] * PE;
      }
      if (choiceM[i, t] == 1){
        Q[2] = Q[2] * d_aF[i];
      }else{
        Q[1] = Q[1] * d_aF[i];
      }
      R = d_v[i] * outcomeM[i, t] + (1-d_v[i]) * R;
      beta = R * d_betaScale[i];
    }
  }
}
generated quantities {
  // For group level parameters
  real<lower=0, upper=1> mu_aN;
  real<lower=0, upper=1> mu_aP;
  real<lower=0, upper=1> mu_aF;
  real<lower=0, upper=1> mu_v;
  real<lower=0, upper=20> mu_betaScale;
  real<lower=0, upper=1> d_mu_aN;
  real<lower=0, upper=1> d_mu_aP;
  real<lower=0, upper=1> d_mu_aF;
  real<lower=0, upper=1> d_mu_v;
  real<lower=0, upper=20> d_mu_betaScale;

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
  mu_v   = Phi_approx(mu_p[4]);
  mu_betaScale = 20*Phi_approx(mu_p[5]);

  d_mu_aN   = Phi_approx(d_mu_p[1]);
  d_mu_aP   = Phi_approx(d_mu_p[2]);
  d_mu_aF   = Phi_approx(d_mu_p[3]);
  d_mu_v   = Phi_approx(d_mu_p[4]);
  d_mu_betaScale = 20*(Phi_approx(d_mu_p[5]));

  { // local section, this saves time and space
    for (i in 1:N) {
      vector[2] Q; // expected value
      real PE;      // prediction error
      real R;
      real beta;

      // Initialize values
      Q = initQ;
      R = 0;
      beta = 0;
      log_lik[i] = 0;

      for (t in 1:(Tsesh[i])) {
        // compute log likelihood of current trial
        log_lik[i] = log_lik[i] + categorical_logit_lpmf(choice[i, t] | beta * Q);

        // generate posterior prediction for current trial
        y_pred[i, t] = categorical_rng(softmax(beta * Q));

        // prediction error
        PE = outcome[i, t] - Q[choice[i, t]];

        // value updating (learning)
        if (PE < 0){
          Q[choice[i, t]] +=aN[i] * PE;
        }
        else{
          Q[choice[i, t]] += aP[i] * PE;
        }
        if (choice[i, t] == 1){
          Q[2] = Q[2] * aF[i];
        }else{
          Q[1] = Q[1] * aF[i];
        }
        R = v[i] * outcome[i, t] + (1-v[i]) * R;
        beta = R * betaScale[i];
      }
    }
    for (i in 1:M) {
      vector[2] Q; // expected value
      real PE;      // prediction error
      real R;
      real beta;

      // Initialize values
      Q = initQ;
      R = 0;
      beta = 0;
      log_likM[i] = 0;

      for (t in 1:(TseshM[i])) {
        // compute log likelihood of current trial
        log_likM[i] = log_likM[i] + categorical_logit_lpmf(choiceM[i, t] | beta * Q);

        // generate posterior prediction for current trial
        y_predM[i, t] = categorical_rng(softmax(beta * Q));

        // prediction error
        PE = outcomeM[i, t] - Q[choiceM[i, t]];

        // value updating (learning)
        if (PE < 0){
          Q[choiceM[i, t]] += d_aN[i] * PE;
        }
        else{
          Q[choiceM[i, t]] += d_aP[i] * PE;
        }
        if (choiceM[i, t] == 1){
          Q[2] = Q[2] * d_aF[i];
        }else{
          Q[1] = Q[1] * d_aF[i];
        }
        R = d_v[i] * outcomeM[i, t] + (1-d_v[i]) * R;
        beta = R * d_betaScale[i];
      }
    }
  }
}
