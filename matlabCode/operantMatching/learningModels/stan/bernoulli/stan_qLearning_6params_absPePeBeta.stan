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
  vector[N] aN_pr;        // learning rate for NPE
  vector[N] aP_pr;        // learning rate for PPE
  vector[N] aF_pr;        // forgetting rate
  vector[N] aPE_pr;  // inverse temperature updating rate
  vector[N] betaScale_pr;   // inverse temp scale val
  vector[N] betaMin_pr;   // inverse temp scale val

}
transformed parameters {
  // session-level parameters
  vector<lower=0, upper=1>[N] aN;
  vector<lower=0, upper=1>[N] aP;
  vector<lower=0, upper=1>[N] aF;
  vector<lower=0, upper=1>[N] aPE;
  vector<lower=0, upper=20>[N] betaScale;
  vector<lower=0, upper=20>[N] betaMin;

  for (n in 1:N) {
    aN[n]        = Phi_approx(mu_p[1] + sigma[1] * aN_pr[n]);
    aP[n]        = Phi_approx(mu_p[2] + sigma[2] * aP_pr[n]);
    aF[n]        = Phi_approx(mu_p[3] + sigma[3] * aF_pr[n]);
    aPE[n]       = Phi_approx(mu_p[4] + sigma[4] * aPE_pr[n]);
    betaScale[n] = Phi_approx(mu_p[5] + sigma[5] * betaScale_pr[n]) * 20;
    betaMin[n]   = Phi_approx(mu_p[6] + sigma[6] * betaMin_pr[n]) * 20;
  }
}
model {
  // Hyperparameters
  mu_p  ~ normal(0, 1);
  sigma ~ cauchy(0, 1);

  // individual parameters
  aN_pr        ~ normal(0, 1);
  aP_pr        ~ normal(0, 1);
  aF_pr        ~ normal(0, 1);
  aPE_pr       ~ normal(0, 1);
  betaScale_pr ~ normal(0, 1);
  betaMin_pr   ~ normal(0, 1);

  // session loop and trial loop
  for (n in 1:N) {
    vector[2] Q; // expected value
    real PE;      // prediction error
    real pePE; // expected average value
    real peBar;
    vector[Tsesh[n]] beta;
    vector[Tsesh[n]] Qdiff;

    Q = initQ;
    beta[1] = 0;
    peBar = 0;

    for (t in 1:(Tsesh[n])) {
      Qdiff[t] = Q[2] - Q[1];
      choice[n, t] ~ bernoulli_logit(beta[t] * Qdiff[t]);

      if (choice[n,t] == 1) {
        PE = outcome[n, t] - Q[2];
        if (PE < 0){
          Q[2] += aN[n] * PE;
        }else{
          Q[2] += aP[n] * PE;
        }
        Q[1] = Q[1] * aF[n];
      }else{
        PE = outcome[n, t] - Q[1];
        if (PE < 0){
          Q[1] += aN[n] * PE;
        }else{
          Q[1] += aP[n] * PE;
        }
        Q[2] = Q[2] * aF[n];
      }
      if (t < Tsesh[n]){
        pePE = fabs(fabs(PE) - peBar);
        beta[t+1] = (1-pePE) * betaScale[n] + betaMin[n];
        peBar += aPE[n] * (fabs(PE) - peBar);
      }
    }
  }
}
generated quantities {
  // For group level parameters
  real<lower=0, upper=1> mu_aN;
  real<lower=0, upper=1> mu_aP;
  real<lower=0, upper=1> mu_aF;
  real<lower=0, upper=1> mu_aPE;
  real<lower=0, upper=20> mu_betaScale;
  real<lower=0, upper=20> mu_betaMin;

  // For log likelihood calculation
  real log_lik[N];

  mu_aN        = Phi_approx(mu_p[1]);
  mu_aP        = Phi_approx(mu_p[2]);
  mu_aF        = Phi_approx(mu_p[3]);
  mu_aPE       = Phi_approx(mu_p[4]); 
  mu_betaScale = Phi_approx(mu_p[5]) * 20;
  mu_betaMin   = Phi_approx(mu_p[6]) * 20;

  { // local section, this saves time and space
    for (n in 1:N) {
      vector[2] Q; // expected value
      real PE;      // prediction error
      real pePE;
      real peBar;
      vector[Tsesh[n]] beta;
      vector[Tsesh[n]] Qdiff;

      // Initialize values
      Q = initQ;
      beta[1] = 0;
      peBar = 0;

      log_lik[n] = 0;

      for (t in 1:(Tsesh[n])) {
        Qdiff[t] = Q[2] - Q[1];

        // compute log likelihood of current trial
        log_lik[n] = log_lik[n] + bernoulli_logit_lpmf(choice[n, t] | beta[t] * Qdiff[t]);

        if (choice[n,t] == 1) {
          PE = outcome[n, t] - Q[2];
          if (PE < 0){
            Q[2] += aN[n] * PE;
          }else{
            Q[2] += aP[n] * PE;
          }
          Q[1] = Q[1] * aF[n];
        }else{
          PE = outcome[n, t] - Q[1];
          if (PE < 0){
            Q[1] += aN[n] * PE;
          }else{
            Q[1] += aP[n] * PE;
          }
          Q[2] = Q[2] * aF[n];
        }
        if (t < Tsesh[n]){
          pePE = fabs(fabs(PE) - peBar);
          beta[t+1] = (1-pePE) * betaScale[n] + betaMin[n];
          peBar += aPE[n] * (fabs(PE) - peBar);
        }
      }
    }
  }
}
