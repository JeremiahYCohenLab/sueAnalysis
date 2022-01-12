data {
  int<lower=1> N; //number of pairs of sessions in group1
  int<lower=1> T; //max number of trials in all session
  int<lower=1, upper=T> Tsesh[2*N];
  int<lower=0, upper=1> choice[2*N, T];
  real outcome[2*N, T];  // no lower and upper bounds
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
  real sigmaDiff;
  real sigmaCommon;

  // Session-level raw parameters
  vector[N] aN_pr;    // learning rate for NPE
  vector[N] aP_pr;    // learning rate for PPE
  vector[N] aF_pr;    // forgetting rate
  vector[N] beta_pr;  // inverse temperature
  vector[N] diff_pr;  //difference between sessions
  vector[2*N] bias;  // side bias
}
transformed parameters {
// Transform session-level raw parameters
  vector<lower=0, upper=1>[2*N] aN;
  vector<lower=0, upper=1>[2*N] aP;
  vector<lower=0, upper=1>[2*N] aF;
  vector<lower=0, upper=10>[2*N] beta;
  vector [N] diff;
  for (n in 1:N){
    diff[n] = mu_diff + sigmaDiff * diff_pr[n];
  }
  for (n in 1:N) {

    if (cmpParam == 1){
      aP[n+N]   = Phi_approx(mu_p[1] + sigma[1] * aP_pr[n]); aP[n] = aP[n+N];
      aF[n+N]   = Phi_approx(mu_p[2] + sigma[2] * aF_pr[n]); aF[n] = aF[n+N];
      beta[n+N] = Phi_approx(mu_p[3] + sigma[3] * beta_pr[n]) * 10; beta[n] = beta[n+N];

      aN[n] = Phi_approx(mu_common + sigmaCommon * aN_pr[n] + diff[n]);
      aN[n+N] = Phi_approx(mu_common + sigmaCommon * aN_pr[n] - diff[n]);
    }else if (cmpParam == 2){
     // print(mu_common);
     // print(aP_pr[n]);
     // print(sigmaCommon);
     // print(diff[n]);
      aN[n+N]   = Phi_approx(mu_p[1] + sigma[1] * aN_pr[n]); aN[n] = aN[n+N];
      aF[n+N]   = Phi_approx(mu_p[2] + sigma[2] * aF_pr[n]); aF[n] = aF[n+N];
      beta[n+N] = Phi_approx(mu_p[3] + sigma[3] * beta_pr[n]) * 10; beta[n] = beta[n+N];

      aP[n] = Phi_approx(mu_common + sigmaCommon * aP_pr[n] + diff[n]);
      aP[n+N] = Phi_approx(mu_common + sigmaCommon * aP_pr[n] - diff[n]);
    }else if (cmpParam == 3){
      aN[n+N]   = Phi_approx(mu_p[1] + sigma[1] * aN_pr[n]); aN[n] = aN[n+N];
      aP[n+N]   = Phi_approx(mu_p[2] + sigma[2] * aP_pr[n]); aP[n] = aP[n+N];
      beta[n+N] = Phi_approx(mu_p[3] + sigma[3] * beta_pr[n]) * 10; beta[n] = beta[n+N];

      aF[n] = Phi_approx(mu_common + sigmaCommon * aF_pr[n] + diff[n]);
      aF[n+N] = Phi_approx(mu_common + sigmaCommon * aF_pr[n] - diff[n]);
    }else{
      aN[n+N]   = Phi_approx(mu_p[1] + sigma[1] * aN_pr[n]); aN[n] = aN[n+N];
      aP[n+N]   = Phi_approx(mu_p[2] + sigma[2] * aP_pr[n]); aP[n] = aP[n+N];
      aF[n+N]   = Phi_approx(mu_p[3] + sigma[3] * aF_pr[n]); aF[n] = aF[n+N];

      beta[n] = Phi_approx(mu_common + sigmaCommon * beta_pr[n] + diff[n])*10;
      beta[n+N] = Phi_approx(mu_common + sigmaCommon * beta_pr[n] - diff[n])*10;
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
  for (n in 1:2*N) {
    vector[2] Q; // expected value
    real PE;      // prediction error
    vector[Tsesh[n]] Qdiff;

    Q = initQ;

    for (t in 1:(Tsesh[n])) {
      Qdiff[t] = Q[2] - Q[1];
      choice[n, t] ~ bernoulli_logit(beta[n] * Qdiff[t] + bias[n]);

      if (choice[n,t] == 1) {
        PE = outcome[n, t] - Q[2];
        if (PE < 0){
          Q[2] = Q[2] + aN[n] * PE;
        }else{
          Q[2] = Q[2] + aP[n] * PE;
        }
        Q[1] = Q[1] * aF[n];
      }else{
        PE = outcome[n, t] - Q[1];
        if (PE < 0){
          Q[1] = Q[1] + aN[n] * PE;
        }else{
          Q[1] = Q[1] + aP[n] * PE;
        }
        Q[2] = Q[2] * aF[n];
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

  real<lower=0, upper=1> mu_aN_diff;
  real<lower=0, upper=1> mu_aP_diff;
  real<lower=0, upper=1> mu_aF_diff;
  real<lower=0, upper=10> mu_beta_diff;

  if (cmpParam == 1) {
    mu_aP   = Phi_approx(mu_p[1]); mu_aP_diff = mu_aP;
    mu_aF   = Phi_approx(mu_p[2]); mu_aF_diff = mu_aF;
    mu_beta = Phi_approx(mu_p[3])*10; mu_beta_diff = mu_beta;
    mu_aN_diff = Phi_approx(mu_common + mu_diff);
    mu_aN = Phi_approx(mu_common - mu_diff);
  }else if (cmpParam == 2){
    mu_aN   = Phi_approx(mu_p[1]); mu_aN_diff = mu_aN;
    mu_aF   = Phi_approx(mu_p[2]); mu_aF_diff = mu_aF;
    mu_beta = Phi_approx(mu_p[3])*10; mu_beta_diff = mu_beta;
    mu_aP_diff = Phi_approx(mu_common + mu_diff);
    mu_aP = Phi_approx(mu_common - mu_diff);
  }else if (cmpParam ==3) {
    mu_aN   = Phi_approx(mu_p[1]); mu_aN_diff = mu_aN;
    mu_aP   = Phi_approx(mu_p[2]); mu_aP_diff = mu_aP;
    mu_beta = Phi_approx(mu_p[3])*10; mu_beta_diff = mu_beta;
    mu_aF_diff = Phi_approx(mu_common + mu_diff);
    mu_aF = Phi_approx(mu_common - mu_diff);
  }else{
    mu_aN   = Phi_approx(mu_p[1]); mu_aN_diff = mu_aN;
    mu_aP   = Phi_approx(mu_p[2]); mu_aP_diff = mu_aP;
    mu_aF   = Phi_approx(mu_p[3]); mu_aF_diff = mu_aF;
    mu_beta_diff = Phi_approx(mu_common + mu_diff)*10;
    mu_beta = Phi_approx(mu_common - mu_diff)*10;
  }

  // For log likelihood calculation
  real log_lik[2*N];
  real log_likMean[2*N];



  { // local section, this saves time and space
    for (n in 1:2*N) {
      vector[2] Q; // expected value
      real PE;      // prediction error
      vector[Tsesh[n]] Qdiff;
      vector[2] Q_bias; //biased Q value

      // Initialize values
      Q = initQ;

      log_lik[n] = 0;

      for (t in 1:(Tsesh[n])) {
        // calculate bias Q 
        Q_bias = Q;
        Q_bias[2] = Q_bias[2] + bias[n]/beta[n];
        // calculate Qdiff
        Qdiff[t] = Q[2] - Q[1];

        // compute log likelihood of current trial
        log_lik[n] = log_lik[n] + bernoulli_logit_lpmf(choice[n, t] | (beta[n] * Qdiff[t] + bias[n]));

        // generate posterior prediction for current trial

        //y_pred[n, t] = categorical_rng(softmax(beta[n] * Q_bias));

        if (choice[n,t] == 1) {
          PE = outcome[n, t] - Q[2];
          if (PE < 0){
            Q[2] = Q[2] + aN[n] * PE;
          }else{
            Q[2] = Q[2] + aP[n] * PE;
          }
          Q[1] = Q[1] * aF[n];
        }else{
          PE = outcome[n, t] - Q[1];
          if (PE < 0){
            Q[1] = Q[1] + aN[n] * PE;
          }else{
            Q[1] = Q[1] + aP[n] * PE;
          }
          Q[2] = Q[2] * aF[n];
        }
      }
      log_likMean[n] = log_lik[n]/Tsesh[n];

    }
  }
}
