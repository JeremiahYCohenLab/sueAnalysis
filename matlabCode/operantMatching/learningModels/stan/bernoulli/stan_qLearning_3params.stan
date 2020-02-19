data {
  int<lower=1> N;
  int<lower=1> T;
  int<lower=1, upper=T> Tsesh[N];
  int<lower=0, upper=2> choice[N, T];
  real outcome[N, T];  // no lower and upper bounds
}
transformed data {
  vector[2] initQ;  // initial values for Q
  initQ = rep_vector(0.0, 2);
}
parameters {
// Declare all parameters as vectors for vectorizing
  // Hyper(animal)-parameters
  vector[3] mu_p;
  vector<lower=0>[3] sigma;

  // Session-level raw parameters
  vector[N] a_pr;    // learning rate
  vector[N] aF_pr;    // forgetting rate
  vector[N] beta_pr;  // inverse temperature
}
transformed parameters {
// Transform session-level raw parameters
  vector<lower=0, upper=1>[N] a;
  vector<lower=0, upper=1>[N] aF;
  vector<lower=0, upper=10>[N] beta;

  for (n in 1:N) {
    a[n]    = Phi_approx(mu_p[1] + sigma[1] * a_pr[n]);
    aF[n]   = Phi_approx(mu_p[2] + sigma[2] * aF_pr[n]);
    beta[n] = Phi_approx(mu_p[3] + sigma[3] * beta_pr[n]) * 10;
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

  // session loop and trial loop
  for (n in 1:N) {
    vector[2] Q; // expected value
    real PE;      // prediction error
    vector[Tsesh[n]] Qdiff;

    Q = initQ;

    for (t in 1:(Tsesh[n])) {
      Qdiff[t] = Q[2] - Q[1];
      choice[n, t] ~ bernoulli_logit(beta[n] * Qdiff[t]);

      if (choice[n,t] == 1) {
        PE = outcome[n, t] - Q[2];
        Q[2] += a[n] * PE;
        Q[1] = Q[1] * aF[n];
      }else{
        PE = outcome[n, t] - Q[1];
        Q[1] += a[n] * PE;
        Q[2] = Q[2] * aF[n];
      }
    }
  }
}
generated quantities {
  // For group level parameters
  real<lower=0, upper=1> mu_a;
  real<lower=0, upper=1> mu_aF;
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

  mu_a    = Phi_approx(mu_p[1]);
  mu_aF   = Phi_approx(mu_p[2]);
  mu_beta = Phi_approx(mu_p[3])*10;

  { // local section, this saves time and space
    for (n in 1:N) {
      vector[2] Q; // expected value
      real PE;      // prediction error
      vector[Tsesh[n]] Qdiff;

      // Initialize values
      Q = initQ;

      log_lik[n] = 0;

      for (t in 1:(Tsesh[n])) {
        Qdiff[t] = Q[2] - Q[1];

        // compute log likelihood of current trial
        log_lik[n] = log_lik[n] + bernoulli_logit_lpmf(choice[n, t] | beta[n] * Qdiff[t]);

        // generate posterior prediction for current trial
        y_pred[n, t] = categorical_rng(softmax(beta[n] * Q));

        if (choice[n,t] == 1) {
          PE = outcome[n, t] - Q[2];
          Q[2] += a[n] * PE;
          Q[1] = Q[1] * aF[n];
        }else{
          PE = outcome[n, t] - Q[1];
          Q[1] += a[n] * PE;
          Q[2] = Q[2] * aF[n];
        }
      }
    }
  }
}
