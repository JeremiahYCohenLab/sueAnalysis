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
  vector[N] a_pr;    // learning rate for NPE
  vector[N] aF_pr;    // forgetting rate
  vector[N] beta_pr;  // inverse temperature
  vector[N] k_pr; // chocie autoCorr
  vector[N] diff_pr;  //difference between sessions
  vector[2*N] bias;  // side bias
}
transformed parameters {
// Transform session-level raw parameters
  vector<lower=0, upper=1>[2*N] a;
  vector<lower=0, upper=1>[2*N] aF;
  vector<lower=0, upper=10>[2*N] beta;
  vector<lower=-5, upper=5>[2*N] k;
  vector [N] diff;
  for (n in 1:N){
    diff[n] = mu_diff + sigmaDiff * diff_pr[n];
  }
  for (n in 1:N) {

    if (cmpParam == 1){
      aF[n+N]   = Phi_approx(mu_p[1] + sigma[1] * aF_pr[n]); aF[n] = aF[n+N];
      beta[n+N] = Phi_approx(mu_p[2] + sigma[2] * beta_pr[n]) * 10; beta[n] = beta[n+N];
      k[n+N]   = 10*Phi_approx(mu_p[3] + sigma[3] * k_pr[n])-5; k[n] = k[n+N];

      a[n] = Phi_approx(mu_common + sigmaCommon * a_pr[n] + diff[n]);
      a[n+N] = Phi_approx(mu_common + sigmaCommon * a_pr[n] - diff[n]);
    }else if (cmpParam == 2){
      a[n+N]   = Phi_approx(mu_p[1] + sigma[1] * a_pr[n]); a[n] = a[n+N];
      beta[n+N] = Phi_approx(mu_p[2] + sigma[2] * beta_pr[n]) * 10; beta[n] = beta[n+N];
      k[n+N]   = 10*Phi_approx(mu_p[3] + sigma[3] * k_pr[n])-5; k[n] = k[n+N];

      aF[n] = Phi_approx(mu_common + sigmaCommon * aF_pr[n] + diff[n]);
      aF[n+N] = Phi_approx(mu_common + sigmaCommon * aF_pr[n] - diff[n]);
    }else if (cmpParam == 3){
      a[n+N]   = Phi_approx(mu_p[1] + sigma[1] * a_pr[n]); a[n] = a[n+N];
      aF[n+N]   = Phi_approx(mu_p[2] + sigma[2] * aF_pr[n]); aF[n] = aF[n+N];
      k[n+N]   = 10*Phi_approx(mu_p[3] + sigma[3] * k_pr[n])-5; k[n] = k[n+N];

      beta[n] = Phi_approx(mu_common + sigmaCommon * beta_pr[n] + diff[n]) * 10;
      beta[n+N] = Phi_approx(mu_common + sigmaCommon * beta_pr[n] - diff[n]) * 10;
    }else{
      a[n+N]   = Phi_approx(mu_p[1] + sigma[1] * a_pr[n]); a[n] = a[n+N];
      aF[n+N]   = Phi_approx(mu_p[2] + sigma[2] * aF_pr[n]); aF[n] = aF[n+N];
      beta[n+N] = Phi_approx(mu_p[3] + sigma[3] * beta_pr[n]) * 10; beta[n] = beta[n+N];

      k[n] = Phi_approx(mu_common + sigmaCommon * k_pr[n] + diff[n])*10 - 5;
      k[n+N] = Phi_approx(mu_common + sigmaCommon * k_pr[n] - diff[n])*10 - 5;
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
  a_pr   ~ normal(0, 1);
  aF_pr   ~ normal(0, 1);
  beta_pr ~ normal(0, 1);
  k_pr ~ normal(0, 1);
  diff_pr ~ normal(0, 1);
  bias    ~ normal(0, 10);

  // session loop and trial loop
  for (n in 1:2*N) {
    vector[2] Q; // expected value
    real PE;      // prediction error
    vector[Tsesh[n]] Qdiff;
    real prevChoice;

    Q = initQ;
    prevChoice = 0;

    for (t in 1:(Tsesh[n])) {
      Qdiff[t] = Q[2] - Q[1];
      choice[n, t] ~ bernoulli_logit(beta[n] * Qdiff[t] + k[n] * prevChoice + bias[n]);

      if (choice[n,t] == 1) {
        PE = outcome[n, t] - Q[2];
        Q[2] = Q[2] + a[n] * PE;
        Q[1] = Q[1] * aF[n];
        prevChoice = 1;
      }else{
        PE = outcome[n, t] - Q[1];
        Q[1] = Q[1] + a[n] * PE;
        Q[2] = Q[2] * aF[n];
        prevChoice = -1;
      }
    }
  }
}
generated quantities {
  // For group level parameters
  real<lower=0, upper=1> mu_a;
  real<lower=0, upper=1> mu_aF;
  real<lower=0, upper=10> mu_beta;
  real<lower=-5, upper=5> mu_k;

  real<lower=0, upper=1> mu_a_diff;
  real<lower=0, upper=1> mu_aF_diff;
  real<lower=0, upper=10> mu_beta_diff;
  real<lower=-5, upper=5> mu_k_diff;

  if (cmpParam == 1) {
    mu_aF   = Phi_approx(mu_p[1]); mu_aF_diff = mu_aF;
    mu_beta = Phi_approx(mu_p[2])*10; mu_beta_diff = mu_beta;
    mu_k = Phi_approx(mu_p[3])*10-5; mu_k_diff = mu_k;
    mu_a_diff = Phi_approx(mu_common + mu_diff);
    mu_a = Phi_approx(mu_common - mu_diff);
  }else if (cmpParam == 2){
    mu_a   = Phi_approx(mu_p[1]); mu_a_diff = mu_a;
    mu_beta = Phi_approx(mu_p[2])*10; mu_beta_diff = mu_beta;
    mu_k = Phi_approx(mu_p[3])*10-5; mu_k_diff = mu_k;
    mu_aF_diff = Phi_approx(mu_common + mu_diff);
    mu_aF = Phi_approx(mu_common - mu_diff);
  }else if (cmpParam ==3) {
    mu_a   = Phi_approx(mu_p[1]); mu_a_diff = mu_a;
    mu_aF   = Phi_approx(mu_p[2]); mu_aF_diff = mu_aF;
    mu_k = Phi_approx(mu_p[3])*10 - 5; mu_k_diff = mu_k;
    mu_beta_diff = Phi_approx(mu_common + mu_diff)*10;
    mu_beta = Phi_approx(mu_common - mu_diff)*10;
  }else{
    mu_a   = Phi_approx(mu_p[1]); mu_a_diff = mu_a;
    mu_aF   = Phi_approx(mu_p[2]); mu_aF_diff = mu_aF;
    mu_beta = Phi_approx(mu_p[3])*10; mu_beta_diff = mu_beta;
    mu_k_diff = Phi_approx(mu_common + mu_diff)*10-5;
    mu_k = Phi_approx(mu_common - mu_diff)*10-5;
  }

  // For log likelihood calculation
  real log_lik[2*N];
  real log_likMean[2*N];



  { // local section, this saves time and space
    for (n in 1:2*N) {
      vector[2] Q; // expected value
      real PE;      // prediction error
      real prevChoice;
      vector[Tsesh[n]] Qdiff;

      // Initialize values
      Q = initQ;
      prevChoice = 0;

      log_lik[n] = 0;

      for (t in 1:(Tsesh[n])) {

        Qdiff[t] = Q[2] - Q[1];

        // compute log likelihood of current trial
        log_lik[n] = log_lik[n] + bernoulli_logit_lpmf(choice[n, t] | beta[n] * Qdiff[t] + k[n] * prevChoice + bias[n]);

        // generate posterior prediction for current trial
        //y_pred[n, t] = inv_logit(beta[n] * Qdiff[t] + k[n] * prevChoice + bias[n]);

        if (choice[n,t] == 1) {
          PE = outcome[n, t] - Q[2];
          Q[2] = Q[2] + a[n] * PE;
          Q[1] = Q[1] * aF[n];
          prevChoice = 1;


        }else{
          PE = outcome[n, t] - Q[1];
          Q[1] = Q[1] + a[n] * PE;
          Q[2] = Q[2] * aF[n];
          prevChoice = -1;

          }
      log_likMean[n] = log_lik[n]/Tsesh[n];
      }
    }
  }
}
