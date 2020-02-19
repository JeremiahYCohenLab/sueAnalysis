data {
  int<lower=1> M;   // number of mice
  int<lower=1> N;   // number of sessions
  int<lower=1> T;   // max number of trials across all sessions
  int<lower=1> MxN;   // number of mice x sessions
  int<lower=1> MxP;   // number of mice x params
  int<lower=1, upper=T> Tsesh[MxN];         // number of trials in each session
  int<lower=0, upper=2> choice[MxN, T];     // choice array
  int<lower=0, upper=1> outcome[MxN, T];    // outcome array
}
transformed data {
  vector[2] initQ;  // initial values for Q
  initQ = rep_vector(0.0, 2);
}
parameters {
// Declare all parameters as vectors for vectorizing
  // Hyper-hyperparameters
  vector[12] h_mu_p;
  vector<lower=0>[12] h_sigma;

  // Hyperparameters
  vector[MxP] mu_p_pr;
  vector<lower=0>[MxP] sigma_pr;


  // Session-level raw parameters
  vector[MxN] aNscale_pr; // learning rate for negative prediction error (NPE)
  vector[MxN] aNmin_pr;   // minimum NPE LR
  vector[MxN] aP_pr;      // learning rate for positive prediction error (PPE)
  vector[MxN] aF_pr;      // forgetting rate
  vector[MxN] aPE_pr;     // expected uncertainty update rate
  vector[MxN] beta_pr;    // inverse temperature for softmax (decision) function

}
transformed parameters {
  //transformed hyperparameters
  vector[MxP] mu_p;
  vector[MxP] sigma;

  // session-level parameters
  vector<lower=0, upper=2>[MxN]  aNscale;
  vector<lower=0, upper=1>[MxN]  aNmin;
  vector<lower=0, upper=1>[MxN]  aP;
  vector<lower=0, upper=1>[MxN]  aF;
  vector<lower=0, upper=1>[MxN]  aPE;
  vector<lower=0, upper=10>[MxN] beta;

  for (m in 1:M){
    for (p in 1:6){
        mu_p[(m - 1)*6 + p]  = h_mu_p[p] + h_sigma[p] * mu_p_pr[(m - 1)*6 + p];
        sigma[(m - 1)*6 + p] = h_mu_p[p+6] + h_sigma[p+6] * sigma_pr[(m - 1)*6 + p];
    }
  }
  for (m in 1:M){
    for (n in 1:MxN) {
      aNscale[n] = Phi_approx(mu_p[(m - 1)*6 + 1] + sigma[(m - 1)*6 + 1] * aNscale_pr[n]);
      aNmin[n]   = Phi_approx(mu_p[(m - 1)*6 + 2] + sigma[(m - 1)*6 + 2] * aNmin_pr[n]);
      aP[n]      = Phi_approx(mu_p[(m - 1)*6 + 3] + sigma[(m - 1)*6 + 3] * aP_pr[n]);
      aF[n]      = Phi_approx(mu_p[(m - 1)*6 + 4] + sigma[(m - 1)*6 + 4] * aF_pr[n]);
      aPE[n]     = Phi_approx(mu_p[(m - 1)*6 + 5] + sigma[(m - 1)*6 + 5] * aPE_pr[n]);
      beta[n]    = Phi_approx(mu_p[(m - 1)*6 + 6] + sigma[(m - 1)*6 + 6] * beta_pr[n]) * 10;
    }
  }
}
model {
  // hyper-hyperparameters
  h_mu_p  ~ normal(0, 1);
  h_sigma ~ cauchy(0, 1);

  // hyperparameters
  mu_p_pr  ~ normal(0, 1);
  sigma_pr ~ cauchy(0, 1);

  // individual parameters
  aNscale_pr ~ normal(0, 1);
  aNmin_pr   ~ normal(0, 1);
  aP_pr      ~ normal(0, 1);
  aF_pr      ~ normal(0, 1);
  aPE_pr     ~ normal(0, 1);
  beta_pr    ~ normal(0, 1);

  // session x mouse loop
  for (n in 1:MxN) {
    vector[2] Q; // expected value
    real PE;      // prediction error
    real pePE;    // prediction error prediction error
    real peBar;   // expected uncertainty
    real aN;      // NPE learning rate
    vector[Tsesh[MxN]] Qdiff; // difference in action values

    Q = initQ;
    peBar = 0;
    aN = aNmin[n];

    //trial loop
    for (t in 1:(Tsesh[MxN])) {
      // compute action probabilities
      Qdiff[t] = Q[2] - Q[1];
      choice[n,t] ~ bernoulli_logit(beta[n] * Qdiff[t]);

      if (choice[n,t] == 1) {
        PE = outcome[n,t] - Q[2];
        pePE = fabs(fabs(PE) - peBar);
        aN = (1 - pePE) * aNscale[n] + aNmin[n];
        if (PE < 0){
          Q[2] += aN * PE;
        }else{
          Q[2] += aP[n] * PE;
        }
        Q[1] = Q[1] * aF[n];
      }else{
        PE = outcome[n,t] - Q[1];
        pePE = fabs(fabs(PE) - peBar);
        aN = (1 - pePE) * aNscale[n] + aNmin[n];
        if (PE < 0){
          Q[1] += aN * PE;
        }else{
          Q[1] += aP[n] * PE;
        }
        Q[2] = Q[2] * aF[n];
      }
      peBar += aPE[n] * (fabs(PE) - peBar);
    }
  }
}
generated quantities {
  // For mouse level parameters
  vector<lower=0, upper=2>[M] mu_aNscale;
  vector<lower=0, upper=1>[M] mu_aNmin;
  vector<lower=0, upper=1>[M] mu_aP;
  vector<lower=0, upper=1>[M] mu_aF;
  vector<lower=0, upper=1>[M] mu_aPE;
  vector<lower=0, upper=10>[M] mu_beta;

  // for population level parameters
  real<lower=0, upper=2> h_mu_aNscale;
  real<lower=0, upper=1> h_mu_aNmin;
  real<lower=0, upper=1> h_mu_aP;
  real<lower=0, upper=1> h_mu_aF;
  real<lower=0, upper=1> h_mu_aPE;
  real<lower=0, upper=10> h_mu_beta;

  // For log likelihood calculation
  real log_lik[MxN];

  for (m in 1:M){
    mu_aNscale[m] = Phi_approx(mu_p[(m - 1)*6 + 1]);
    mu_aNmin[m]   = Phi_approx(mu_p[(m - 1)*6 + 2]);
    mu_aP[m]      = Phi_approx(mu_p[(m - 1)*6 + 3]);
    mu_aF[m]      = Phi_approx(mu_p[(m - 1)*6 + 4]);
    mu_aPE[m]     = Phi_approx(mu_p[(m - 1)*6 + 5]);
    mu_beta[m]    = Phi_approx(mu_p[(m - 1)*6 + 6]) * 10;
  }

  h_mu_aNscale = Phi_approx(h_mu_p[1]);
  h_mu_aNmin   = Phi_approx(h_mu_p[2]);
  h_mu_aP      = Phi_approx(h_mu_p[3]);
  h_mu_aF      = Phi_approx(h_mu_p[4]);
  h_mu_aPE     = Phi_approx(h_mu_p[5]);
  h_mu_beta    = Phi_approx(h_mu_p[6]) * 10;

  { // local section, this saves time and space
    // session x mouse loop
    for (n in 1:MxN) {
      vector[2] Q; // expected value
      real PE;      // prediction error
      real pePE;    // prediction error prediction error
      real peBar;   // expected average value
      real aN;      // NPE learning rate
      vector[Tsesh[MxN]] Qdiff; // difference in action values

      Q = initQ;
      peBar = 0;
      aN = aNmin[n];

      //trial loop
      for (t in 1:(Tsesh[MxN])) {
        // compute action probabilities
        Qdiff[t] = Q[2] - Q[1];
        // compute log likelihood of current trial
        log_lik[n] = log_lik[n] + bernoulli_logit_lpmf(choice[n, t] | beta[n] * Qdiff[t]);

        if (choice[n,t] == 1) {
          PE = outcome[n,t] - Q[2];
          pePE = fabs(fabs(PE) - peBar);
          aN = (1 - pePE) * aNscale[n] + aNmin[n];
          if (PE < 0){
            Q[2] += aN * PE;
          }else{
            Q[2] += aP[n] * PE;
          }
          Q[1] = Q[1] * aF[n];
        }else{
          PE = outcome[n,t] - Q[1];
          pePE = fabs(fabs(PE) - peBar);
          aN = (1 - pePE) * aNscale[n] + aNmin[n];
          if (PE < 0){
            Q[1] += aN * PE;
          }else{
            Q[1] += aP[n] * PE;
          }
          Q[2] = Q[2] * aF[n];
        }
        peBar += aPE[n] * (fabs(PE) - peBar);
      }
    }
  }
}
