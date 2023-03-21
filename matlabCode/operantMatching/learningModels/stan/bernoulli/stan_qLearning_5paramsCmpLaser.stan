data {
  int<lower=1> N; //number of sessions
  int<lower=1> T; //max number of trials in all session
  int<lower=1, upper=T> Tsesh[N];
  int<lower=0, upper=1> choice[N, T];
  int<lower=0, upper=1> laser[N, T];
  real outcome[N, T];  // no lower and upper bounds
  int<lower=1> cmpParam; 
}
transformed data {
  vector[2] initQ;  // initial values for Q
  initQ = rep_vector(0.0, 2);
}
parameters {
// Declare all parameters as vectors for vectorizing
  // Hyper(animal)-parameters
  vector[3] mu_p;
  real mu_diff;
  real mu_common;

  vector<lower=0>[3] sigma;
  real<lower=0> sigmaDiff;
  real<lower=0> sigmaCommon;

  // Session-level raw parameters
  vector[N] aN_pr;    // learning rate for NPE
  vector[N] aP_pr;    // learning rate for PPE
  vector[N] aF_pr;    // forgetting rate
  vector[N] beta_pr;  // inverse temperature
  vector[N] diff_pr;  //difference between sessions
  vector[N] bias;  // side bias
}
transformed parameters {
// Transform session-level raw parameters
  vector<lower=0, upper=1>[N] aN_L;
  vector<lower=0, upper=1>[N] aP_L;
  vector<lower=0, upper=1>[N] aF_L;
  vector<lower=0, upper=10>[N] beta_L;
  vector<lower=0, upper=1>[N] aN_H;
  vector<lower=0, upper=1>[N] aP_H;
  vector<lower=0, upper=1>[N] aF_H;
  vector<lower=0, upper=10>[N] beta_H;
  vector [N] diff;
  for (n in 1:N){
    diff[n] = mu_diff + sigmaDiff * diff_pr[n];
  }
  for (n in 1:N) {

    if (cmpParam == 1){
      aP_L[n]   = Phi_approx(mu_p[1] + sigma[1] * aP_pr[n]); aP_H[n] = aP_L[n];
      aF_L[n]   = Phi_approx(mu_p[2] + sigma[2] * aF_pr[n]); aF_H[n] = aF_L[n];
      beta_L[n] = Phi_approx(mu_p[3] + sigma[3] * beta_pr[n]) * 10; beta_H[n] = beta_L[n];

      aN_L[n] = Phi_approx(mu_common + sigmaCommon * aN_pr[n] - diff[n]);
      aN_H[n] = Phi_approx(mu_common + sigmaCommon * aN_pr[n] + diff[n]);
    }else if (cmpParam == 2){
     // print(mu_common);
     // print(aP_pr[n]);
     // print(sigmaCommon);
     // print(diff[n]);
      aN_L[n]   = Phi_approx(mu_p[1] + sigma[1] * aN_pr[n]); aN_H[n] = aN_L[n];
      aF_L[n]   = Phi_approx(mu_p[2] + sigma[2] * aF_pr[n]); aF_H[n] = aF_L[n];
      beta_L[n] = Phi_approx(mu_p[3] + sigma[3] * beta_pr[n]) * 10; beta_H[n] = beta_L[n];

      aP_L[n] = Phi_approx(mu_common + sigmaCommon * aP_pr[n] - diff[n]);
      aP_H[n] = Phi_approx(mu_common + sigmaCommon * aP_pr[n] + diff[n]);
    }else if (cmpParam == 3){
      aN_L[n]   = Phi_approx(mu_p[1] + sigma[1] * aN_pr[n]); aN_H[n] = aN_L[n];
      aP_L[n]   = Phi_approx(mu_p[2] + sigma[2] * aP_pr[n]); aP_H[n] = aP_L[n];
      beta_L[n] = Phi_approx(mu_p[3] + sigma[3] * beta_pr[n]) * 10; beta_H[n] = beta_L[n];

      aF_L[n] = Phi_approx(mu_common + sigmaCommon * aF_pr[n] - diff[n]);
      aF_H[n] = Phi_approx(mu_common + sigmaCommon * aF_pr[n] + diff[n]);
    }else{
      aN_L[n]   = Phi_approx(mu_p[1] + sigma[1] * aN_pr[n]); aN_H[n] = aN_L[n];
      aP_L[n]   = Phi_approx(mu_p[2] + sigma[2] * aP_pr[n]); aP_H[n] = aP_L[n];
      aF_L[n]   = Phi_approx(mu_p[3] + sigma[3] * aF_pr[n]); aF_H[n] = aF_L[n];

      beta_L[n] = Phi_approx(mu_common + sigmaCommon * beta_pr[n] - diff[n])*10;
      beta_H[n] = Phi_approx(mu_common + sigmaCommon * beta_pr[n] + diff[n])*10;
    }

  }
}

model {
  // Hyperparameters
  mu_p        ~ normal(0, 1);
  mu_diff     ~ normal(0, pow(0.5, 0.5));
  mu_common   ~ normal(0, pow(0.5, 0.5));
  sigma       ~ cauchy(0, 0.2);
  sigmaDiff   ~ cauchy(0, 0.2);
  sigmaCommon ~ cauchy(0, 0.1);

  // individual parameters
  aN_pr   ~ normal(0, 1);
  aP_pr   ~ normal(0, 1);
  aF_pr   ~ normal(0, 1);
  beta_pr ~ normal(0, 1);
  diff_pr ~ normal(0, 1);
  bias    ~ normal(0, 10);

  // session loop and trial loop
  for (n in 1:N) {
    vector[2] Q; // expected value
    real PE;      // prediction error
    vector[Tsesh[n]] Qdiff;

    Q = initQ;

    for (t in 1:(Tsesh[n])) {
      Qdiff[t] = Q[2] - Q[1];
      if (laser[n,t] == 0) {
        choice[n, t] ~ bernoulli_logit(beta_L[n] * Qdiff[t] + bias[n]);
        if (choice[n,t] == 1) {
          PE = outcome[n, t] - Q[2];
          if (PE < 0){
            Q[2] = Q[2] + aN_L[n] * PE;
          }else{
            Q[2] = Q[2] + aP_L[n] * PE;
          }
          Q[1] = Q[1] * aF_L[n];
        }else{
          PE = outcome[n, t] - Q[1];
          if (PE < 0){
            Q[1] = Q[1] + aN_L[n] * PE;
          }else{
            Q[1] = Q[1] + aP_L[n] * PE;
          }
          Q[2] = Q[2] * aF_L[n];
        }
      }else{
        choice[n, t] ~ bernoulli_logit(beta_H[n] * Qdiff[t] + bias[n]);
        if (choice[n,t] == 1) {
          PE = outcome[n, t] - Q[2];
          if (PE < 0){
            Q[2] = Q[2] + aN_H[n] * PE;
          }else{
            Q[2] = Q[2] + aP_H[n] * PE;
          }
          Q[1] = Q[1] * aF_H[n];
        }else{
          PE = outcome[n, t] - Q[1];
          if (PE < 0){
            Q[1] = Q[1] + aN_H[n] * PE;
          }else{
            Q[1] = Q[1] + aP_H[n] * PE;
          }
          Q[2] = Q[2] * aF_H[n];
        }
      }
    }
  }
}

generated quantities {
  // For group level parameters
  real<lower=0, upper=1> mu_aN_L;
  real<lower=0, upper=1> mu_aP_L;
  real<lower=0, upper=1> mu_aF_L;
  real<lower=0, upper=10> mu_beta_L;

  real<lower=0, upper=1> mu_aN_H;
  real<lower=0, upper=1> mu_aP_H;
  real<lower=0, upper=1> mu_aF_H;
  real<lower=0, upper=10> mu_beta_H;

  if (cmpParam == 1) {
    mu_aP_L   = Phi_approx(mu_p[1]); mu_aP_H = mu_aP_L;
    mu_aF_L   = Phi_approx(mu_p[2]); mu_aF_H = mu_aF_L;
    mu_beta_L = Phi_approx(mu_p[3])*10; mu_beta_H = mu_beta_L;
    mu_aN_H = Phi_approx(mu_common + mu_diff);
    mu_aN_L = Phi_approx(mu_common - mu_diff);
  }else if (cmpParam == 2){
    mu_aN_L   = Phi_approx(mu_p[1]); mu_aN_H = mu_aN_L;
    mu_aF_L   = Phi_approx(mu_p[2]); mu_aF_H = mu_aF_L;
    mu_beta_L = Phi_approx(mu_p[3])*10; mu_beta_H = mu_beta_L;
    mu_aP_H = Phi_approx(mu_common + mu_diff);
    mu_aP_L = Phi_approx(mu_common - mu_diff);
  }else if (cmpParam ==3) {
    mu_aN_L   = Phi_approx(mu_p[1]); mu_aN_H = mu_aN_L;
    mu_aP_L   = Phi_approx(mu_p[2]); mu_aP_H = mu_aP_L;
    mu_beta_L = Phi_approx(mu_p[3])*10; mu_beta_H = mu_beta_L;
    mu_aF_H = Phi_approx(mu_common + mu_diff);
    mu_aF_L = Phi_approx(mu_common - mu_diff);
  }else{
    mu_aN_L   = Phi_approx(mu_p[1]); mu_aN_H = mu_aN_L;
    mu_aP_L   = Phi_approx(mu_p[2]); mu_aP_H = mu_aP_L;
    mu_aF_L   = Phi_approx(mu_p[3]); mu_aF_H = mu_aF_L;
    mu_beta_H = Phi_approx(mu_common + mu_diff)*10;
    mu_beta_L = Phi_approx(mu_common - mu_diff)*10;
  }

  // For log likelihood calculation
  real log_lik[N];
  real log_likMean[N];



  { // local section, this saves time and space
    for (n in 1:N) {
      vector[2] Q; // expected value
      real PE;      // prediction error
      vector[Tsesh[n]] Qdiff;
      vector[2] Q_bias; //biased Q value

      // Initialize values
      Q = initQ;

      log_lik[n] = 0;

      for (t in 1:(Tsesh[n])) {
        if (laser[n,t] == 0){
          // calculate Qdiff
          Qdiff[t] = Q[2] - Q[1];

          // compute log likelihood of current trial
          log_lik[n] = log_lik[n] + bernoulli_logit_lpmf(choice[n, t] | (beta_L[n] * Qdiff[t] + bias[n]));

          // generate posterior prediction for current trial

          //y_pred[n, t] = categorical_rng(softmax(beta[n] * Q_bias));

          if (choice[n,t] == 1) {
            PE = outcome[n, t] - Q[2];
            if (PE < 0){
              Q[2] = Q[2] + aN_L[n] * PE;
            }else{
              Q[2] = Q[2] + aP_L[n] * PE;
            }
            Q[1] = Q[1] * aF_L[n];
          }else{
            PE = outcome[n, t] - Q[1];
            if (PE < 0){
              Q[1] = Q[1] + aN_L[n] * PE;
            }else{
              Q[1] = Q[1] + aP_L[n] * PE;
            }
            Q[2] = Q[2] * aF_L[n];
          }
        }else{
          // calculate Qdiff
          Qdiff[t] = Q[2] - Q[1];

          // compute log likelihood of current trial
          log_lik[n] = log_lik[n] + bernoulli_logit_lpmf(choice[n, t] | (beta_H[n] * Qdiff[t] + bias[n]));

          // generate posterior prediction for current trial

          //y_pred[n, t] = categorical_rng(softmax(beta[n] * Q_bias));

          if (choice[n,t] == 1) {
            PE = outcome[n, t] - Q[2];
            if (PE < 0){
              Q[2] = Q[2] + aN_H[n] * PE;
            }else{
              Q[2] = Q[2] + aP_H[n] * PE;
            }
            Q[1] = Q[1] * aF_H[n];
          }else{
            PE = outcome[n, t] - Q[1];
            if (PE < 0){
              Q[1] = Q[1] + aN_H[n] * PE;
            }else{
              Q[1] = Q[1] + aP_H[n] * PE;
            }
            Q[2] = Q[2] * aF_H[n];
          }

        }
      }
      log_likMean[n] = log_lik[n]/Tsesh[n];

    }
  }
}
