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
  vector[N] aNscale_pr;        // learning rate for NPE
  vector[N] aNmin_pr;        // learning rate for NPE
  vector[N] aP_pr;        // learning rate for PPE
  vector[N] aF_pr;        // forgetting rate
  vector[N] aPE_pr;        // learning rate for volatility
  vector[N] beta_pr;   // inverse temp scale val

}
transformed parameters {
  // session-level parameters
  vector<lower=0, upper=1>[N] aNscale;
  vector<lower=0, upper=1>[N] aNmin;
  vector<lower=0, upper=1>[N] aP;
  vector<lower=0, upper=1>[N] aF;
  vector<lower=0, upper=1>[N] aPE;
  vector<lower=0, upper=10>[N] beta;

  for (n in 1:N) {
    aNscale[n] = Phi_approx(mu_p[1] + sigma[1] * aNscale_pr[n]);
    aNmin[n]   = Phi_approx(mu_p[2] + sigma[2] * aNmin_pr[n]);
    aP[n]      = Phi_approx(mu_p[3] + sigma[3] * aP_pr[n]);
    aF[n]      = Phi_approx(mu_p[4] + sigma[4] * aF_pr[n]);
    aPE[n]     = Phi_approx(mu_p[5] + sigma[5] * aPE_pr[n]);
    beta[n]    = Phi_approx(mu_p[6] + sigma[6] * beta_pr[n]) * 10;
  }
}
model {
  // Hyperparameters
  mu_p  ~ normal(0, 0.8);
  sigma ~ cauchy(0, 1);

  // individual parameters
  aNscale_pr ~ normal(0, 1);
  aNmin_pr   ~ normal(0, 1);
  aP_pr      ~ normal(0, 1);
  aF_pr      ~ normal(0, 1);
  aPE_pr     ~ normal(0, 1);
  beta_pr    ~ normal(0, 1);

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
        pePE = PE - peBar;
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
        pePE = PE - peBar;
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
      peBar += aPE[n] * pePE;
    }
  }
}
generated quantities {
  // For group level parameters
  real<lower=0, upper=1> mu_aNscale;
  real<lower=0, upper=1> mu_aNmin;
  real<lower=0, upper=1> mu_aP;
  real<lower=0, upper=1> mu_aF;
  real<lower=0, upper=1> mu_aPE;
  real<lower=0, upper=10> mu_beta;

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

  mu_aNscale = Phi_approx(mu_p[1]);
  mu_aNmin   = Phi_approx(mu_p[2]);
  mu_aP      = Phi_approx(mu_p[3]);
  mu_aF      = Phi_approx(mu_p[4]);
  mu_aPE     = Phi_approx(mu_p[5]); 
  mu_beta    = Phi_approx(mu_p[6]) * 10;

  { // local section, this saves time and space
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

      log_lik[n] = 0;

      for (t in 1:(Tsesh[n])) {
        Qdiff[t] = Q[2] - Q[1];
        // compute log likelihood of current trial
        log_lik[n] += bernoulli_logit_lpmf(choice[n, t] | beta[n] * Qdiff[t]);

        // generate posterior prediction for current trial
  //      y_pred[n, t] = categorical_rng(softmax(beta[n] * Q));

        if (choice[n,t] == 1) {
          PE = outcome[n, t] - Q[2];
          pePE = PE - peBar;
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
          pePE = PE - peBar;
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
        peBar += aPE[n] * pePE;
      }
    }
  }
}
