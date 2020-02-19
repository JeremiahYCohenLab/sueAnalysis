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
  vector[N] aN_pr;        // learning rate for NPE
  vector[N] aP_pr;        // learning rate for PPE
  vector[N] aF_pr;        // forgetting rate
  vector[N] v_pr;  // inverse temperature updating rate
  vector[N] betaScale_pr;   // inverse temp min val
  vector[N] betaMin_pr;   // inverse temp min val
  vector[N] k_pr;   // choice autocorrelation

}
transformed parameters {
  // session-level parameters
  vector<lower=0, upper=1>[N] aN;
  vector<lower=0, upper=1>[N] aP;
  vector<lower=0, upper=1>[N] aF;
  vector<lower=0, upper=1>[N] v;
  vector<lower=0, upper=20>[N] betaScale;
  vector<lower=0, upper=20>[N] betaMin;
  vector<lower=-5, upper=5>[N] k;

  for (n in 1:N) {
    aN[n]        = Phi_approx(mu_p[1] + sigma[1] * aN_pr[n]);
    aP[n]        = Phi_approx(mu_p[2] + sigma[2] * aP_pr[n]);
    aF[n]        = Phi_approx(mu_p[3] + sigma[3] * aF_pr[n]);
    v[n]         = Phi_approx(mu_p[4] + sigma[4] * v_pr[n]);
    betaScale[n] = Phi_approx(mu_p[5] + sigma[5] * betaScale_pr[n]) * 20;
    betaMin[n]   = Phi_approx(mu_p[6] + sigma[6] * betaMin_pr[n]) * 20;
    k[n]         = Phi_approx(mu_p[7] + sigma[7] * k_pr[n]) * 10 - 5;
  }
}
model {
  // Hyperparameters
  mu_p  ~ normal(0, 1);
  sigma ~ cauchy(0, 5);

  // individual parameters
  aN_pr        ~ normal(0, 1);
  aP_pr        ~ normal(0, 1);
  aF_pr        ~ normal(0, 1);
  v_pr         ~ normal(0, 1);
  betaScale_pr ~ normal(0, 1);
  betaMin_pr   ~ normal(0, 1);
  k_pr         ~ normal(0, 1);

  // session loop and trial loop
  for (n in 1:N) {
    vector[2] Q; // expected value
    real prevChoice;
    real PE;      // prediction error
    real peBar; // expected average value
    vector[Tsesh[n]] beta;
    vector[Tsesh[n]] Qdiff;

    Q = initQ;
    prevChoice = 0;
    peBar = 0;
    beta[1] = betaMin[n];

    for (t in 1:(Tsesh[n])) {
      // compute action probabilities
      Qdiff[t] = Q[2] - Q[1];
      choice[n, t] ~ bernoulli_logit(beta[t] * Qdiff[t] + k[n] * prevChoice);

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
      peBar = v[n] * abs(PE) + (1-v[n]) * peBar;
      if (t < Tsesh[n]){
        beta[t+1] = betaMin[n] + (1-peBar) * betaScale[n];
      }
    }
  }
}
generated quantities {
  // For group level parameters
  real<lower=0, upper=1> mu_aN;
  real<lower=0, upper=1> mu_aP;
  real<lower=0, upper=1> mu_aF;
  real<lower=0, upper=1> mu_v;
  real<lower=0, upper=20> mu_betaScale;
  real<lower=0, upper=20> mu_betaMin;
  real<lower=-5, upper=5> mu_k;

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

  mu_aN        = Phi_approx(mu_p[1]);
  mu_aP        = Phi_approx(mu_p[2]);
  mu_aF        = Phi_approx(mu_p[3]);
  mu_v         = Phi_approx(mu_p[4]); 
  mu_betaScale = Phi_approx(mu_p[5]) * 20;
  mu_betaMin   = Phi_approx(mu_p[6]) * 20;
  mu_k         = Phi_approx(mu_p[7]) * 10 - 5; 

  { // local section, this saves time and space
    for (n in 1:N) {
      vector[2] Q; // expected value
      real prevChoice;
      real PE;      // prediction error
      real peBar;
      vector[Tsesh[n]] beta;
      vector[Tsesh[n]] Qdiff;

      // Initialize values
      Q = initQ;
      prevChoice = 0;
      peBar = 0;
      beta[1] = betaMin[n];

      log_lik[n] = 0;

      for (t in 1:(Tsesh[n])) {
        Qdiff[t] = Q[2] - Q[1];

        // compute log likelihood of current trial
        log_lik[n] = log_lik[n] + bernoulli_logit_lpmf(choice[n, t] | beta[t] * Qdiff[t] + k[n] * prevChoice);

        // generate posterior prediction for current trial
     //   y_pred[n, t] = bernoulli_rng(softmax(beta[t] * (Qdiff[t] + k[n]*prevChoice)));

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
        peBar = v[n] * abs(PE) + (1-v[n]) * peBar;
        if (t < Tsesh[n]){
          beta[t+1] = betaMin[n] + (1-peBar) * betaScale[n];
        }
      }
    }
  }
}
