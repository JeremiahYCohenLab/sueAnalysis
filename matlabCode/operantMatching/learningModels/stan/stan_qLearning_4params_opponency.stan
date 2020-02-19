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
  vector[4] mu_p;
  vector<lower=0>[4] sigma;

  // Session-level raw parameters
  vector[N] a_pr;    // learning rate for NPE
  vector[N] aF_pr;    // forgetting rate
  vector[N] beta_pr;  // inverse temperature
  vector[N] v_pr;     // expected average value learning rate
}
transformed parameters {
  // session-level parameters
  vector<lower=0, upper=1>[N] a;
  vector<lower=0, upper=1>[N] aF;
  vector<lower=0, upper=20>[N] beta;
  vector<lower=0, upper=1>[N] v;

  for (i in 1:N) {
    a[i]   = Phi_approx(mu_p[1]  + sigma[1]  * a_pr[i]);
    aF[i]   = Phi_approx(mu_p[2]  + sigma[2]  * aF_pr[i]);
    beta[i] = Phi_approx(mu_p[3] + sigma[3] * beta_pr[i]) * 20;
    v[i]   = Phi_approx(mu_p[4]  + sigma[4]  * v_pr[i]);
  }
}
model {
  // Hyperparameters
  mu_p  ~ normal(0, 1);
  sigma ~ cauchy(0, 1);

  // individual parameters
  a_pr   ~ normal(0, 1);
  aF_pr   ~ normal(0, 1);
  beta_pr ~ normal(0, 1);
  v_pr    ~ normal(0, 1);

  // session loop and trial loop
  for (i in 1:N) {
    vector[2] Q; // expected value
    real PE;      // prediction error
    real rBar; // expected average value

    Q = initQ;
    rBar = 0.4;

    for (t in 1:(Tsesh[i])) {
      // compute action probabilities
      choice[i, t] ~ categorical_logit(beta[i] * Q);

      // prediction error
      PE = outcome[i, t] - Q[choice[i, t]] - rBar;

      // value updating (learning)
      Q[choice[i, t]] = Q[choice[i, t]] + a[i] * PE;
      if (choice[i, t] == 1){
        Q[2] = Q[2] * aF[i];
      }else{
        Q[1] = Q[1] * aF[i];
      }
      rBar = v[i] * outcome[i, t] + (1-v[i]) * rBar;
    }
  }
}
generated quantities {
  // For group level parameters
  real<lower=0, upper=1> mu_a;
  real<lower=0, upper=1> mu_aF;
  real<lower=0, upper=20> mu_beta;
  real<lower=0, upper=1> mu_v;

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

  mu_a   = Phi_approx(mu_p[1]);
  mu_aF   = Phi_approx(mu_p[2]);
  mu_beta = Phi_approx(mu_p[3]) * 20;
  mu_v   = Phi_approx(mu_p[4]);

  { // local section, this saves time and space
    for (i in 1:N) {
      vector[2] Q; // expected value
      real PE;      // prediction error
      real rBar; // expected average value

      // Initialize values
      Q = initQ;
      rBar = 0.4;

      log_lik[i] = 0;

      for (t in 1:(Tsesh[i])) {
        // compute log likelihood of current trial
        log_lik[i] = log_lik[i] + categorical_logit_lpmf(choice[i, t] | beta[i] * Q);

        // generate posterior prediction for current trial
        y_pred[i, t] = categorical_rng(softmax(beta[i] * Q));

        // prediction error
        PE = outcome[i, t] - Q[choice[i, t]] - rBar;

        // value updating (learning)
        Q[choice[i, t]] = Q[choice[i, t]] + a[i] * PE;
        if (choice[i, t] == 1){
          Q[2] = Q[2] * aF[i];
        }else{
          Q[1] = Q[1] * aF[i];
        }
        rBar = v[i] * outcome[i, t] + (1-v[i]) * rBar;
      }
    }
  }
}
