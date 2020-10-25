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
  vector[6] mu_p_pr;
  vector<lower=0>[6] sigma_pr;
  vector[6] d_mu_p_pr;
  vector<lower=0>[6] d_sigma_pr;

  // Session-level raw parameters
  vector[N] aNscale_pr;    // learning rate scale factor for negative prediction error (NPE)
  vector[N] aNmin_pr;    // minimum learning rate for negative prediction error (NPE)
  vector[N] aP_pr;    // learning rate for positive prediction error (PPE)
  vector[N] aF_pr;    // forgetting rate
  vector[M] aPE_pr;    // learning rate for expected uncertainty
  vector[N] beta_pr;  // inverse temperature for softmax (decision) function

  vector[M] d_aNscale_pr;    // learning rate scale factor for negative prediction error (NPE)
  vector[M] d_aNmin_pr;    // minimum learning rate for negative prediction error (NPE)
  vector[M] d_aP_pr;    // learning rate for PPE
  vector[M] d_aF_pr;    // forgetting rate
  vector[M] d_aPE_pr;    // learning rate for expected uncertainty
  vector[M] d_beta_pr;  // inverse temperature
}
transformed parameters {
  //hyperparameters
  vector[6] mu_p;
  vector[6] sigma;
  vector[6] d_mu_p;
  vector[6] d_sigma;

  // session-level parameters
  vector<lower=0, upper=2>[N] aNscale;
  vector<lower=0, upper=2>[N] aNmin;
  vector<lower=0, upper=1>[N] aP;
  vector<lower=0, upper=1>[N] aF;
  vector<lower=0, upper=1>[N] aPE;
//  vector<lower=0>[N] aPE;
//  vector[N] aPE;
  vector<lower=0, upper=10>[N] beta;

  vector<lower=0, upper=2>[M] d_aNscale;
  vector<lower=0, upper=2>[M] d_aNmin;
  vector<lower=0, upper=1>[M] d_aP;
  vector<lower=0, upper=1>[M] d_aF;
  vector<lower=0, upper=1>[N] d_aPE;
//  vector<lower=0>[M] d_aPE;
//  vector[N] d_aPE;
  vector<lower=0, upper=10>[M] d_beta;

  for (n in 1:6){
    mu_p[n] = h_mu_p[1]  + h_sigma[1]  * mu_p_pr[n];
    sigma[n] = h_mu_p[2]  + h_sigma[2]  * sigma_pr[n];
    d_mu_p[n] = h_mu_p[1]  + h_sigma[1]  * d_mu_p_pr[n];
    d_sigma[n] = h_mu_p[2]  + h_sigma[2]  * d_sigma_pr[n];
  }

  for (n in 1:N) {
    aNscale[n] = Phi_approx(mu_p[1] + sigma[1] * aNscale_pr[n]) * 2;
    aNmin[n]   = Phi_approx(mu_p[2] + sigma[2] * aNmin_pr[n]) * 2;
    aP[n]      = Phi_approx(mu_p[3] + sigma[3] * aP_pr[n]);
    aF[n]      = Phi_approx(mu_p[4] + sigma[4] * aF_pr[n]);
    aPE[n]     = Phi_approx(mu_p[5] + sigma[5] * aPE_pr[n]);
//    aPE[n]     = mu_p[5] + sigma[5] * aPE_pr[n];
    beta[n]    = Phi_approx(mu_p[6] + sigma[6] * beta_pr[n]) * 10;
    }
  for (n in 1:M){
    d_aNscale[n] = Phi_approx(d_mu_p[1] + d_sigma[1] * d_aNscale_pr[n]) * 2;
    d_aNmin[n]   = Phi_approx(d_mu_p[2] + d_sigma[2] * d_aNmin_pr[n]) * 2;
    d_aP[n]      = Phi_approx(d_mu_p[3] + d_sigma[3] * d_aP_pr[n]);
    d_aF[n]      = Phi_approx(d_mu_p[4] + d_sigma[4] * d_aF_pr[n]);
    d_aPE[n]     = Phi_approx(d_mu_p[5] + d_sigma[5] * d_aPE_pr[n]);
//    d_aPE[n]     = d_mu_p[5] + d_sigma[5] * d_aPE_pr[n];
    d_beta[n]    = Phi_approx(d_mu_p[6] + d_sigma[6] * d_beta_pr[n]) * 10;
  }
}
model {
  //hyper-hyperparameters
  h_mu_p  ~ normal(0, 10);
  h_sigma ~ cauchy(0, 10);

  // Hyperparameters
  mu_p_pr    ~ normal(0, 1);
  sigma_pr   ~ cauchy(0, 1);
  d_mu_p_pr  ~ normal(0, 1);
  d_sigma_pr ~ cauchy(0, 1);

  // individual parameters
  aNscale_pr ~ normal(0, 1);
  aNmin_pr   ~ normal(0, 1);
  aP_pr      ~ normal(0, 1);
  aF_pr      ~ normal(0, 1);
  aPE_pr     ~ normal(0, 1);
//  aPE_pr     ~ gamma(2, 0);
  beta_pr    ~ normal(0, 1);

  d_aNscale_pr ~ normal(0, 1);
  d_aNmin_pr   ~ normal(0, 1);
  d_aP_pr      ~ normal(0, 1);
  d_aF_pr      ~ normal(0, 1);
  d_aPE_pr     ~ normal(0, 1);
//  d_aPE_pr     ~ gamma(2, 0);
  d_beta_pr    ~ normal(0, 1);


  // session loop and trial loop
  for (n in 1:N) {
    vector[2] Q;  // expected value
    real PE;      // prediction error
    real pePE;    // prediction error prediction error
    real peBar;   // expected average value
    real aN;      // NPE learning rate
    vector[Tsesh[n]] Qdiff;

    Q = initQ;
    peBar = 0;
    aN = aNmin[n];

    for (t in 1:(Tsesh[n])) {
      Qdiff[t] = Q[2] - Q[1];
      choice[n, t] ~ bernoulli_logit(beta[n] * Qdiff[t]);

      if (choice[n,t] == 1) {
        PE = outcome[n, t] - Q[2];
        pePE = fabs(PE) - peBar;
        aN = pePE * aNscale[n] + aNmin[n];
        if (aN < 0){
          aN = 0;
        }
        if (PE < 0){
          Q[2] += aN * PE;
        }else{
          Q[2] += aP[n] * PE;
        }
        Q[1] = Q[1] * aF[n];
      }else{
        PE = outcome[n, t] - Q[1];
        pePE = fabs(PE) - peBar;
        aN = pePE * aNscale[n] + aNmin[n];
        if (aN < 0){
          aN = 0;
        }
        if (PE < 0){
          Q[1] += aN * PE;
        }else{
          Q[1] += aP[n] * PE;
        }
        Q[2] = Q[2] * aF[n];
      }
      peBar += aPE[n] * (pePE);
    }
  }
  for (n in 1:M) {
    vector[2] Q;  // expected value
    real PE;      // prediction error
    real pePE;    // prediction error prediction error
    real peBar;   // expected average value
    real aN;      // NPE learning rate
    vector[TseshM[n]] Qdiff;

    Q = initQ;
    peBar = 0;
    aN = d_aNmin[n];

    for (t in 1:(TseshM[n])) {
      Qdiff[t] = Q[2] - Q[1];
      choiceM[n, t] ~ bernoulli_logit(d_beta[n] * Qdiff[t]);

      if (choiceM[n,t] == 1) {
        PE = outcomeM[n, t] - Q[2];
        pePE = fabs(PE) - peBar;
        aN = pePE * d_aNscale[n] + d_aNmin[n];
        if (aN < 0){
          aN = 0;
        }
        if (PE < 0){
          Q[2] += aN * PE;
        }else{
          Q[2] += d_aP[n] * PE;
        }
        Q[1] = Q[1] * d_aF[n];
      }else{
        PE = outcomeM[n, t] - Q[1];
        pePE = fabs(PE) - peBar;
        aN = pePE * d_aNscale[n] + d_aNmin[n];
        if (aN < 0){
          aN = 0;
        }
        if (PE < 0){
          Q[1] += aN * PE;
        }else{
          Q[1] += d_aP[n] * PE;
        }
        Q[2] = Q[2] * d_aF[n];
      }
      peBar += d_aPE[n] * (pePE);
    }
  }
}
generated quantities {
  // For group level parameters
  real<lower=0, upper=2> mu_aNscale;
  real<lower=0, upper=2> mu_aNmin;
  real<lower=0, upper=1> mu_aP;
  real<lower=0, upper=1> mu_aF;
  real<lower=0, upper=1> mu_aPE;
//  real<lower=0> mu_aPE;
//  real mu_aPE;
  real<lower=0, upper=10> mu_beta;

  real<lower=0, upper=2> d_mu_aNscale;
  real<lower=0, upper=2> d_mu_aNmin;
  real<lower=0, upper=1> d_mu_aP;
  real<lower=0, upper=1> d_mu_aF;
  real<lower=0, upper=1> d_mu_aPE;
//  real<lower=0> d_mu_aPE;
//  real d_mu_aPE;
  real<lower=0, upper=10> d_mu_beta;

  // For log likelihood calculation
  real log_lik[N];
  real log_likM[M];

  mu_aNscale = Phi_approx(mu_p[1]) * 2;
  mu_aNmin   = Phi_approx(mu_p[2]) * 2;
  mu_aP      = Phi_approx(mu_p[3]);
  mu_aF      = Phi_approx(mu_p[4]);
  mu_aPE     = Phi_approx(mu_p[5]);
//  mu_aPE     = mu_p[5];
  mu_beta    = Phi_approx(mu_p[6]) * 10;

  d_mu_aNscale = Phi_approx(d_mu_p[1]) * 2;
  d_mu_aNmin   = Phi_approx(d_mu_p[2]) * 2;
  d_mu_aP      = Phi_approx(d_mu_p[3]);
  d_mu_aF      = Phi_approx(d_mu_p[4]);
  d_mu_aPE     = Phi_approx(d_mu_p[5]);
//  d_mu_aPE     = d_mu_p[5];
  d_mu_beta    = Phi_approx(d_mu_p[6]) * 10;

  { // local section, this saves time and space
    // session loop and trial loop
    for (n in 1:N) {
      vector[2] Q;  // expected value
      real PE;      // prediction error
      real pePE;    // prediction error prediction error
      real peBar;   // expected average value
      real aN;      // NPE learning rate
      vector[Tsesh[n]] Qdiff;

      Q = initQ;
      peBar = 0;
      aN = aNmin[n];

      for (t in 1:(Tsesh[n])) {
        Qdiff[t] = Q[2] - Q[1];
        // compute log likelihood of current trial
        log_lik[n] += bernoulli_logit_lpmf(choice[n, t] | beta[n] * Qdiff[t]);

        if (choice[n,t] == 1) {
          PE = outcome[n, t] - Q[2];
          pePE = fabs(PE) - peBar;
          aN = pePE * aNscale[n] + aNmin[n];
          if (aN < 0){
            aN = 0;
          }
          if (PE < 0){
            Q[2] += aN * PE;
          }else{
            Q[2] += aP[n] * PE;
          }
          Q[1] = Q[1] * aF[n];
        }else{
          PE = outcome[n, t] - Q[1];
          pePE = fabs(PE) - peBar;
          aN = pePE * aNscale[n] + aNmin[n];
          if (aN < 0){
            aN = 0;
          }
          if (PE < 0){
            Q[1] += aN * PE;
          }else{
            Q[1] += aP[n] * PE;
          }
          Q[2] = Q[2] * aF[n];
        }
        peBar += aPE[n] * (pePE);
      }
    }
    for (n in 1:M) {
      vector[2] Q;  // expected value
      real PE;      // prediction error
      real pePE;    // prediction error prediction error
      real peBar;   // expected average value
      real aN;      // NPE learning rate
      vector[TseshM[n]] Qdiff;

      Q = initQ;
      peBar = 0;
      aN = d_aNmin[n];

      for (t in 1:(TseshM[n])) {
        Qdiff[t] = Q[2] - Q[1];
        // compute log likelihood of current trial
        log_likM[n] += bernoulli_logit_lpmf(choiceM[n, t] | d_beta[n] * Qdiff[t]);

        if (choiceM[n,t] == 1) {
          PE = outcomeM[n, t] - Q[2];
          pePE = fabs(PE) - peBar;
          aN = pePE * d_aNscale[n] + d_aNmin[n];
          if (aN < 0){
            aN = 0;
          }
          if (PE < 0){
            Q[2] += aN * PE;
          }else{
            Q[2] += d_aP[n] * PE;
          }
          Q[1] = Q[1] * d_aF[n];
        }else{
          PE = outcomeM[n, t] - Q[1];
          pePE = fabs(PE) - peBar;
          aN = pePE * d_aNscale[n] + d_aNmin[n];
          if (aN < 0){
            aN = 0;
          }
          if (PE < 0){
            Q[1] += aN * PE;
          }else{
            Q[1] += d_aP[n] * PE;
          }
          Q[2] = Q[2] * d_aF[n];
        }
        peBar += d_aPE[n] * (pePE);
      }
    }
  }
}