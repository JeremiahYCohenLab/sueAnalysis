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
  vector[N] v_pr;  // inverse temperature updating rate
  vector[N] betaScale_pr;   // inverse temp scale val
  vector[N] betaInit_pr;    // inverse temp init val

}
transformed parameters {
  // session-level parameters
  vector<lower=0, upper=1>[N] aN;
  vector<lower=0, upper=1>[N] aP;
  vector<lower=0, upper=1>[N] aF;
  vector<lower=0, upper=1>[N] v;
  vector<lower=0, upper=20>[N] betaScale;
  vector<lower=0, upper=20>[N] betaInit;

  for (i in 1:N) {
    aN[i]   = Phi_approx(mu_p[1]  + sigma[1]  * aN_pr[i]);
    aP[i]   = Phi_approx(mu_p[2]  + sigma[2]  * aP_pr[i]);
    aF[i]   = Phi_approx(mu_p[3]  + sigma[3]  * aF_pr[i]);
    v[i] = Phi_approx(mu_p[4] + sigma[4] * v_pr[i]) ;
    betaScale[i] = Phi_approx(mu_p[5] + sigma[5] * betaScale_pr[i]) * 20;
    betaInit[i] = Phi_approx(mu_p[6] + sigma[6] * betaInit_pr[i]) * 20;
  }
}
model {
  // Hyperparameters
  mu_p  ~ normal(0, 1);
  sigma ~ cauchy(0, 1);

  // individual parameters
  aN_pr   ~ normal(0, 1);
  aP_pr   ~ normal(0, 1);
  aF_pr   ~ normal(0, 1);
  v_pr ~ normal(0, 1);
  betaScale_pr ~ normal(0, 1);
  betaInit_pr ~ normal(0, 1);

  // session loop and trial loop
  for (i in 1:N) {
    vector[2] Q; // expected value
    real PE;      // prediction error
    real R; // expected average value
    real beta;

    Q = initQ;
    R = 0;
    beta = betaInit[i];

    for (t in 1:(Tsesh[i])) {
      // compute action probabilities
      choice[i, t] ~ categorical_logit(beta * Q);

      // prediction error
      PE = outcome[i, t] - Q[choice[i, t]];

      // value updating (learning)
      if (PE < 0){
        Q[choice[i, t]] = Q[choice[i, t]] + aN[i] * PE;
      }
      else{
        Q[choice[i, t]] = Q[choice[i, t]] + aP[i] * PE;
      }
      if (choice[i, t] == 1){
        Q[2] = Q[2] * aF[i];
      }else{
        Q[1] = Q[1] * aF[i];
      }
      R = v[i] * outcome[i, t] + (1-v[i]) * R;
      beta = R * betaScale[i];
    }
  }
}
generated quantities {
  // For group level parameters
  real<lower=0, upper=1> mu_aN;
  real<lower=0, upper=1> mu_aP;
  real<lower=0, upper=1> mu_aF;
  real<lower=0, upper=20> mu_v;
  real<lower=0, upper=20> mu_betaScale;
  real<lower=0, upper=20> mu_betaInit;

  // For log likelihood calculation
  real log_lik[N];

  // For posterior predictive check
  real y_pred[N, T];

  // Set all posterior predictions to 0 (avoids NULL values)
  for (i in 1:N) {
    for (t in 1:T) {
      y_pred[i, t] = -1;
    }
  }

  mu_aN   = Phi_approx(mu_p[1]);
  mu_aP   = Phi_approx(mu_p[2]);
  mu_aF   = Phi_approx(mu_p[3]);
  mu_v = Phi_approx(mu_p[4]); 
  mu_betaScale = Phi_approx(mu_p[5]) * 20;
  mu_betaInit = Phi_approx(mu_p[6]) * 20;

  { // local section, this saves time and space
    for (i in 1:N) {
      vector[2] Q; // expected value
      real PE;      // prediction error
      real R;
      real beta;

      // Initialize values
      Q = initQ;
      R = 0;
      beta = betaInit[i];

      log_lik[i] = 0;

      for (t in 1:(Tsesh[i])) {
        // compute log likelihood of current trial
        log_lik[i] = log_lik[i] + categorical_logit_lpmf(choice[i, t] | beta * Q);

        // generate posterior prediction for current trial
        y_pred[i, t] = categorical_rng(softmax(beta * Q));

        // prediction error
        PE = outcome[i, t] - Q[choice[i, t]];

        // value updating (learning)
        if (PE < 0){
          Q[choice[i, t]] = Q[choice[i, t]] + aN[i] * PE;
        }
        else{
          Q[choice[i, t]] = Q[choice[i, t]] + aP[i] * PE;
        }
        if (choice[i, t] == 1){
          Q[2] = Q[2] * aF[i];
        }else{
          Q[1] = Q[1] * aF[i];
        }
        R = v[i] * outcome[i, t] + (1-v[i]) * R;
        beta = R * betaScale[i];
      }
    }
  }
}
