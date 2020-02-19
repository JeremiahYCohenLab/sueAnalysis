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
  vector[7] mu_p;
  vector<lower=0>[7] sigma;

  // Session-level raw parameters
  vector[N] aNscale_pr;        // learning rate for NPE
  vector[N] aNmin_pr;        // learning rate for NPE
  vector[N] aPscale_pr;        // learning rate for PPE
  vector[N] aPmin_pr;        // learning rate for PPE
  vector[N] aF_pr;        // forgetting rate
  vector[N] v_pr;  // volatility estimate updating rate
  vector[N] beta_pr;   // inverse temp

}
transformed parameters {
  // session-level parameters
  vector<lower=0, upper=1>[N]  aNscale;
  vector<lower=0, upper=1>[N]  aNmin;
  vector<lower=0, upper=1>[N]  aPscale;
  vector<lower=0, upper=1>[N]  aPmin;
  vector<lower=0, upper=1>[N]  aF;
  vector<lower=0, upper=1>[N]  v;
  vector<lower=0, upper=20>[N] beta;

  for (n in 1:N) {
    aNscale[n] = Phi_approx(mu_p[1] + sigma[1] * aNscale_pr[n]);
    aNmin[n]   = Phi_approx(mu_p[2] + sigma[2] * aNmin_pr[n]);
    aPscale[n] = Phi_approx(mu_p[3] + sigma[3] * aPscale_pr[n]);
    aPmin[n]   = Phi_approx(mu_p[4] + sigma[4] * aPmin_pr[n]);
    aF[n]      = Phi_approx(mu_p[5] + sigma[5] * aF_pr[n]);
    v[n]       = Phi_approx(mu_p[6] + sigma[6] * v_pr[n]) ;
    beta[n]    = Phi_approx(mu_p[7] + sigma[7] * beta_pr[n]) * 20;
  }
}
model {
  // Hyperparameters
  mu_p  ~ normal(0, 1);
  sigma ~ cauchy(0, 1);

  // individual parameters
  aNscale_pr  ~ normal(0, 1);
  aNmin_pr    ~ normal(0, 1);
  aPscale_pr  ~ normal(0, 1);
  aPmin_pr    ~ normal(0, 1);
  aF_pr       ~ normal(0, 1);
  v_pr        ~ normal(0, 1);
  beta_pr     ~ normal(0, 1);

  // session loop and trial loop
  for (n in 1:N) {
    vector[2] Q; // expected value
    real PE;      // prediction error
    real peBar; // expected average value
    real aN;
    real aP;
    vector[Tsesh[n]] Qdiff;

    Q = initQ;
    peBar = 0;
    aN = aNmin[n];
    aP = aPmin[n];

    for (t in 1:(Tsesh[n])) {
      Qdiff[t] = Q[2] - Q[1];
      choice[n, t] ~ bernoulli_logit(beta[n] * Qdiff[t]);

      if (choice[n,t] == 1) {
        PE = outcome[n, t] - Q[2];
        peBar = v[n] * abs(PE) + (1-v[n]) * peBar;
        aN = peBar * aNscale[n] + aNmin[n];
        aP = peBar * aPscale[n] + aPmin[n];
        if (PE < 0){
          Q[2] += aN * PE;
        }else{
          Q[2] += aP * PE;
        }
        Q[1] = Q[1] * aF[n];
      }else{
      	PE = outcome[n, t] - Q[1];
        peBar = v[n] * abs(PE) + (1-v[n]) * peBar;
        aN = peBar * aNscale[n] + aNmin[n];
        aP = peBar * aPscale[n] + aPmin[n];
        if (PE < 0){
          Q[1] += aN * PE;
        }else{
          Q[1] += aP * PE;
        }
        Q[2] = Q[2] * aF[n];
      }
    }
  }
}
generated quantities {
  // For group level parameters
  real<lower=0, upper=1>  mu_aNscale;
  real<lower=0, upper=1>  mu_aNmin;
  real<lower=0, upper=1>  mu_aPscale;
  real<lower=0, upper=1>  mu_aPmin;
  real<lower=0, upper=1>  mu_aF;
  real<lower=0, upper=1>  mu_v;
  real<lower=0, upper=20> mu_beta;

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
  mu_aPscale = Phi_approx(mu_p[3]);
  mu_aPmin   = Phi_approx(mu_p[4]);
  mu_aF      = Phi_approx(mu_p[5]);
  mu_v       = Phi_approx(mu_p[6]); 
  mu_beta    = Phi_approx(mu_p[7]) * 20;

  { // local section, this saves time and space
    for (n in 1:N) {
      vector[2] Q; // expected value
      real PE;      // prediction error
      real peBar;
      real aN;
      real aP;
      vector[Tsesh[n]] Qdiff;

      Q = initQ;
      peBar = 0;
      aN = aNmin[n];
      aP = aPmin[n];

      log_lik[n] = 0;

      for (t in 1:(Tsesh[n])) {
        Qdiff[t] = Q[2] - Q[1];
        // compute log likelihood of current trial
        log_lik[n] += bernoulli_logit_lpmf(choice[n, t] | beta[n] * Qdiff[t]);

        // generate posterior prediction for current trial
        y_pred[n, t] = categorical_rng(softmax(beta[n] * Q));

        if (choice[n,t] == 1) {
          PE = outcome[n, t] - Q[2];
          peBar = v[n] * abs(PE) + (1-v[n]) * peBar;
          aN = peBar * aNscale[n] + aNmin[n];
          aP = peBar * aPscale[n] + aPmin[n];
          if (PE < 0){
            Q[2] += aN * PE;
          }else{
            Q[2] += aP * PE;
          }
          Q[1] = Q[1] * aF[n];
        }else{
          PE = outcome[n, t] - Q[1];
          peBar = v[n] * abs(PE) + (1-v[n]) * peBar;
          aN = peBar * aNscale[n] + aNmin[n];
          aP = peBar * aPscale[n] + aPmin[n];
          if (PE < 0){
            Q[1] += aN * PE;
          }else{
            Q[1] += aP * PE;
          }
          Q[2] = Q[2] * aF[n];
        }
      }
    }
  }
}
