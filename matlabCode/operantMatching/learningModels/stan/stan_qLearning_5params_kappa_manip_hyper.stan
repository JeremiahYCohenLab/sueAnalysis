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
  vector[2] initCP;
  initQ = rep_vector(0.0, 2);
  initCP = rep_vector(0.0, 2);
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
  vector[N] beta_pr;  // inverse temperature for softmax (decision) function
  vector[N] k_pr;  // choice autocorrelation term

  vector[M] d_aN_pr;    // change in learning rate for NPE
  vector[M] d_aP_pr;    // change in learning rate for PPE
  vector[M] d_aF_pr;    // change in forgetting rate
  vector[M] d_beta_pr;  // change in inverse temperature
  vector[M] d_k_pr;  // change in choice autocorrelation term
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
  vector<lower=0, upper=10>[N] beta;
  vector<lower=-5, upper=5>[N] k;

  vector<lower=0, upper=1>[M] d_aN;
  vector<lower=0, upper=1>[M] d_aP;
  vector<lower=0, upper=1>[M] d_aF;
  vector<lower=0, upper=10>[M] d_beta;
  vector<lower=-5, upper=5>[M] d_k;

  for (n in 1:5){
    mu_p[n] = h_mu_p[1]  + h_sigma[1]  * mu_p_pr[n];
    sigma[n] = h_mu_p[2]  + h_sigma[2]  * sigma_pr[n];
    d_mu_p[n] = h_mu_p[1]  + h_sigma[1]  * d_mu_p_pr[n];
    d_sigma[n] = h_mu_p[2]  + h_sigma[2]  * d_sigma_pr[n];
  }

  for (n in 1:N) {
    aN[n]   = Phi_approx(mu_p[1]  + sigma[1]  * aN_pr[n]);
    aP[n]   = Phi_approx(mu_p[2]  + sigma[2]  * aP_pr[n]);
    aF[n]   = Phi_approx(mu_p[3]  + sigma[3]  * aF_pr[n]);
    beta[n] = Phi_approx(mu_p[4] + sigma[4] * beta_pr[n]) * 10;
    k[n]    = Phi_approx(mu_p[5]  + sigma[5]  * k_pr[n]) * 10 - 5;
    }
  for (n in 1:M){
    d_aN[n]   = Phi_approx(d_mu_p[1]  + d_sigma[1]  * d_aN_pr[n]);
    d_aP[n]   = Phi_approx(d_mu_p[2]  + d_sigma[2]  * d_aP_pr[n]);
    d_aF[n]   = Phi_approx(d_mu_p[3]  + d_sigma[3]  * d_aF_pr[n]);
    d_beta[n] = Phi_approx(d_mu_p[4]  + d_sigma[4]  * d_beta_pr[n]) * 10;
    d_k[n]    = Phi_approx(d_mu_p[5]  + d_sigma[5]  * d_k_pr[n]) * 10 - 5;
  }
}
model {
  //hyper-hyperparameters
  h_mu_p  ~ normal(0, 1);
  h_sigma ~ cauchy(0, 1);

  // Hyperparameters
  mu_p_pr  ~ normal(0, 1);
  sigma_pr ~ cauchy(0, 1);
  d_mu_p_pr  ~ normal(0, 1);
  d_sigma_pr ~ cauchy(0, 1);

  // individual parameters
  aN_pr     ~ normal(0, 1);
  aP_pr     ~ normal(0, 1);
  aF_pr     ~ normal(0, 1);
  beta_pr   ~ normal(0, 1);
  k_pr      ~ normal(0, 1);
  d_aN_pr   ~ normal(0, 1);
  d_aP_pr   ~ normal(0, 1);
  d_aF_pr   ~ normal(0, 1);
  d_beta_pr ~ normal(0, 1);
  d_k_pr    ~ normal(0, 1);



  // session loop and trial loop
  for (n in 1:N) {
    vector[2] Q;        // expected value
    real prevChoice; 
    real PE;           // prediction error
    vector[Tsesh[n]] Qdiff;

    Q = initQ;
    prevChoice = 0;

    for (t in 1:(Tsesh[n])) {
      Qdiff[t] = Q[2] - Q[1];
      choice[n, t] ~ bernoulli_logit(beta[n] * Qdiff[t] + k[n] * prevChoice);

      if (choice[n,t] == 1) {
        PE = outcome[n, t] - Q[2];
        if (PE < 0){
          Q[2] += aN[n] * PE;
        }else{
          Q[2] += aP[n] * PE;
        }
        Q[1] = Q[1] * aF[n];
        prevChoice = 1;
      }else{
        PE = outcome[n, t] - Q[1];
        if (PE < 0){
          Q[1] += aN[n] * PE;
        }else{
          Q[1] += aP[n] * PE;
        }
        Q[2] = Q[2] * aF[n];
        prevChoice = -1;
      }
    }
  }
  for (n in 1:M) {
    vector[2] Q;        // expected value
    real prevChoice; 
    real PE;           // prediction error
    vector[TseshM[n]] Qdiff;

    Q = initQ;
    prevChoice = 0;

    for (t in 1:(TseshM[n])) {
      Qdiff[t] = Q[2] - Q[1];
      choiceM[n, t] ~ bernoulli_logit(d_beta[n] * Qdiff[t] + d_k[n] * prevChoice);

      if (choiceM[n,t] == 1) {
        PE = outcomeM[n, t] - Q[2];
        if (PE < 0){
          Q[2] += d_aN[n] * PE;
        }else{
          Q[2] += d_aP[n] * PE;
        }
        Q[1] = Q[1] * d_aF[n];
        prevChoice = 1;
      }else{
        PE = outcomeM[n, t] - Q[1];
        if (PE < 0){
          Q[1] += d_aN[n] * PE;
        }else{
          Q[1] += d_aP[n] * PE;
        }
        Q[2] = Q[2] * d_aF[n];
        prevChoice = -1;
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
  real<lower=-5, upper=5> mu_k;
  real<lower=0, upper=1> d_mu_aN;
  real<lower=0, upper=1> d_mu_aP;
  real<lower=0, upper=1> d_mu_aF;
  real<lower=0, upper=10> d_mu_beta;
  real<lower=-5, upper=5> d_mu_k;

  // For log likelihood calculation
  real log_lik[N];
  real log_likM[M];

  // For posterior predictive check
  //real y_pred[N, T];
  //real y_predM[M,T];

  // Set all posterior predictions to 0 (avoids NULL values)
 // for (n in 1:N) {
 //   for (t in 1:T) {
 //     y_pred[n, t] = -1;
 //   }
 // }
 // for (n in 1:M) {
 //   for (t in 1:T) {
 //     y_predM[n, t] = -1;
 //   }
 // }

  mu_aN   = Phi_approx(mu_p[1]);
  mu_aP   = Phi_approx(mu_p[2]);
  mu_aF   = Phi_approx(mu_p[3]);
  mu_beta = Phi_approx(mu_p[4]) * 10;
  mu_k    = Phi_approx(mu_p[5]) * 10 - 5;

  d_mu_aN   = Phi_approx(d_mu_p[1]);
  d_mu_aP   = Phi_approx(d_mu_p[2]);
  d_mu_aF   = Phi_approx(d_mu_p[3]);
  d_mu_beta = Phi_approx(d_mu_p[4]) * 10;
  d_mu_k    = Phi_approx(d_mu_p[5]) * 10 - 5;

  { // local section, this saves time and space
    for (n in 1:N) {
      vector[2] Q;        // expected value
      real prevChoice; 
      real PE;           // prediction error
      vector[Tsesh[n]] Qdiff;

      Q = initQ;
      prevChoice = 0;

      log_lik[n] = 0;

      for (t in 1:(Tsesh[n])) {
        Qdiff[t] = Q[2] - Q[1];
        log_lik[n] += bernoulli_logit_lpmf(choice[n, t] | beta[n] * Qdiff[t] + k[n] * prevChoice);

        if (choice[n,t] == 1) {
          PE = outcome[n, t] - Q[2];
          if (PE < 0){
            Q[2] += aN[n] * PE;
          }else{
            Q[2] += aP[n] * PE;
          }
          Q[1] = Q[1] * aF[n];
          prevChoice = 1;
        }else{
          PE = outcome[n, t] - Q[1];
          if (PE < 0){
            Q[1] += aN[n] * PE;
          }else{
            Q[1] += aP[n] * PE;
          }
          Q[2] = Q[2] * aF[n];
          prevChoice = -1;
        }
      }
    }
    for (n in 1:M) {
      vector[2] Q;        // expected value
      real prevChoice; 
      real PE;           // prediction error
      vector[TseshM[n]] Qdiff;

      Q = initQ;
      prevChoice = 0;

      log_likM[n] = 0;

      for (t in 1:(TseshM[n])) {
        Qdiff[t] = Q[2] - Q[1];
        log_likM[n] += bernoulli_logit_lpmf(choiceM[n, t] | d_beta[n] * Qdiff[t] + d_k[n] * prevChoice);

        if (choiceM[n,t] == 1) {
          PE = outcomeM[n, t] - Q[2];
          if (PE < 0){
            Q[2] += d_aN[n] * PE;
          }else{
            Q[2] += d_aP[n] * PE;
          }
          Q[1] = Q[1] * d_aF[n];
          prevChoice = 1;
        }else{
          PE = outcomeM[n, t] - Q[1];
          if (PE < 0){
            Q[1] += d_aN[n] * PE;
          }else{
            Q[1] += d_aP[n] * PE;
          }
          Q[2] = Q[2] * d_aF[n];
          prevChoice = -1;
        }
      }
    }
  }
}
