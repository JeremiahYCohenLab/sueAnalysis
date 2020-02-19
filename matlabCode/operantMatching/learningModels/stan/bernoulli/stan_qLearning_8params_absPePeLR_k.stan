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
  vector[8] mu_p;
  vector<lower=0>[8] sigma;

  // Session-level raw parameters
  vector[N] aNscale_pr;  // learning rate for NPE
  vector[N] aNmin_pr;    // learning rate for NPE
  vector[N] aPscale_pr;  // learning rate for PPE
  vector[N] aPmin_pr;    // learning rate for PPE
  vector[N] aF_pr;       // forgetting rate
  vector[N] aPE_pr;      // learning rate for volatility
  vector[N] beta_pr;     // inverse temp
  vector[N] k_pr;        // choice autocorrelation

}
transformed parameters {
  // session-level parameters
  vector<lower=0, upper=2>[N] aNscale;
  vector<lower=0, upper=1>[N] aNmin;
  vector<lower=0, upper=2>[N] aPscale;
  vector<lower=0, upper=1>[N] aPmin;
  vector<lower=0, upper=1>[N] aF;
  vector<lower=0, upper=1>[N] aPE;
  vector<lower=0, upper=10>[N] beta;
  vector<lower=-5, upper=5>[N] k;

  for (n in 1:N) {
    aNscale[n] = Phi_approx(mu_p[1] + sigma[1] * aNscale_pr[n]) * 2;
    aNmin[n]   = Phi_approx(mu_p[2] + sigma[2] * aNmin_pr[n]);
    aPscale[n] = Phi_approx(mu_p[3] + sigma[3] * aPscale_pr[n]) * 2;
    aPmin[n]   = Phi_approx(mu_p[4] + sigma[4] * aPmin_pr[n]);
    aF[n]      = Phi_approx(mu_p[5] + sigma[5] * aF_pr[n]);
    aPE[n]     = Phi_approx(mu_p[6] + sigma[6] * aPE_pr[n]);
    beta[n]    = Phi_approx(mu_p[7] + sigma[7] * beta_pr[n]) * 10;
    k[n]       = Phi_approx(mu_p[8] + sigma[8] * k_pr[n]) * 10 - 5;
  }
}
model {
  // Hyperparameters
  mu_p  ~ normal(0, 0.5);
  sigma ~ cauchy(0, 0.5);

  // individual parameters
  aNscale_pr ~ normal(0, 1);
  aNmin_pr   ~ normal(0, 1);
  aPscale_pr ~ normal(0, 1);
  aPmin_pr   ~ normal(0, 1);
  aF_pr      ~ normal(0, 1);
  aPE_pr     ~ normal(0, 1);
  beta_pr    ~ normal(0, 1);
  k_pr       ~ normal(0, 1);

  // session loop and trial loop
  for (n in 1:N) {
    vector[2] Q;  // expected value
    real PE;      // prediction error
    real pePE;    // prediction error prediction error
    real peBar;   // expected average value
    real aN;      // NPE learning rate
    real aP;      // PPE learning rate
    real prevChoice;
    vector[Tsesh[n]] Qdiff;

    Q = initQ;
    prevChoice = 0;
    peBar = 0;
    aN = aNmin[n];
    aP = aPmin[n];

    for (t in 1:(Tsesh[n])) {
      Qdiff[t] = Q[2] - Q[1];
      choice[n, t] ~ bernoulli_logit(beta[n] * Qdiff[t] + k[n] * prevChoice);

      if (choice[n,t] == 1) {
        PE = outcome[n, t] - Q[2];
        pePE = fabs(fabs(PE) - peBar);
        if (PE < 0){
          aN = pePE * aNscale[n] + aNmin[n];
          Q[2] += aN * PE;
        }else{
          aP = pePE * aPscale[n] + aPmin[n];
          Q[2] += aP * PE;
        }
        Q[1] = Q[1] * aF[n];
        prevChoice = 1;
      }else{
        PE = outcome[n, t] - Q[1];
        pePE = fabs(fabs(PE) - peBar);
        if (PE < 0){
          aN = pePE * aNscale[n] + aNmin[n];
          Q[1] += aN * PE;
        }else{
          aP = pePE * aPscale[n] + aPmin[n];
          Q[1] += aP * PE;
        }
        Q[2] = Q[2] * aF[n];
        prevChoice = -1;
      }
      peBar += aPE[n] * (fabs(PE) - peBar);
    }
  }
}
generated quantities {
  // For group level parameters
  real<lower=0, upper=2> mu_aNscale;
  real<lower=0, upper=1> mu_aNmin;
  real<lower=0, upper=2> mu_aPscale;
  real<lower=0, upper=1> mu_aPmin;
  real<lower=0, upper=1> mu_aF;
  real<lower=0, upper=1> mu_aPE;
  real<lower=0, upper=10> mu_beta;
  real<lower=-5, upper=5> mu_k;

  // For log likelihood calculation
  real log_lik[N];

  mu_aNscale = Phi_approx(mu_p[1]) * 2;
  mu_aNmin   = Phi_approx(mu_p[2]);
  mu_aPscale = Phi_approx(mu_p[3]) * 2;
  mu_aPmin   = Phi_approx(mu_p[4]);
  mu_aF      = Phi_approx(mu_p[5]);
  mu_aPE     = Phi_approx(mu_p[6]); 
  mu_beta    = Phi_approx(mu_p[7]) * 10;
  mu_k       = Phi_approx(mu_p[8]) * 10 - 5;

  { // local section, this saves time and space
    for (n in 1:N) {
      vector[2] Q;  // expected value
      real PE;      // prediction error
      real pePE;    // prediction error prediction error
      real peBar;   // expected average value
      real aN;      // NPE learning rate
      real aP;      // PPE learning rate
      real prevChoice;
      vector[Tsesh[n]] Qdiff;

      Q = initQ;
      peBar = 0;
      aN = aNmin[n];
      aP = aPmin[n];
      prevChoice = 0;

      log_lik[n] = 0;

      for (t in 1:(Tsesh[n])) {
        Qdiff[t] = Q[2] - Q[1];
        // compute log likelihood of current trial
        log_lik[n] = log_lik[n] + bernoulli_logit_lpmf(choice[n, t] | beta[n] * Qdiff[t] + k[n] * prevChoice);

        if (choice[n,t] == 1) {
          PE = outcome[n, t] - Q[2];
          pePE = fabs(fabs(PE) - peBar);
          if (PE < 0){
            aN = pePE * aNscale[n] + aNmin[n];
            Q[2] += aN * PE;
          }else{
            aP = pePE * aPscale[n] + aPmin[n];
            Q[2] += aP * PE;
          }
          Q[1] = Q[1] * aF[n];
          prevChoice = 1;
        }else{
          PE = outcome[n, t] - Q[1];
          pePE = fabs(fabs(PE) - peBar);
          if (PE < 0){
            aN = pePE * aNscale[n] + aNmin[n];
            Q[1] += aN * PE;
          }else{
            aP = pePE * aPscale[n] + aPmin[n];
            Q[1] += aP * PE;
          }
          Q[2] = Q[2] * aF[n];
          prevChoice = -1;
        }
        peBar += aPE[n] * (fabs(PE) - peBar);
      }
    }
  }
}
