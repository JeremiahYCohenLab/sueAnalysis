data {
  int<lower=1> N;
  int<lower=1> T;
  int<lower=1, upper=T> Tsesh[N];
  int<lower=0, upper=1> choice[N, T];
  real outcome[N, T];  // no lower and upper bounds
  real<lower=0, upper=1> aN;
  real<lower=0, upper=1> aF;
  real<lower=0, upper=10> beta;
}
transformed data {
  vector[2] initQ;  // initial values for Q
  initQ = rep_vector(0.0, 2);
}
parameters {
// Declare all parameters as vectors for vectorizing
  // Hyper(animal)-parameters
  real mu_p;
  real<lower=0> sigma;

  // Session-level raw parameters
  vector[N] aP_pr;    // learning rate for PPE
}
transformed parameters {
// Transform session-level raw parameters
  vector<lower=0, upper=1>[N] aP;

  for (n in 1:N) {
    aP[n]   = Phi_approx(mu_p + sigma * aP_pr[n]);
  }
}
model {
  // Hyperparameters
  mu_p  ~ normal(0, 1);
  sigma ~ cauchy(0, 1);

  // individual parameters
  aP_pr   ~ normal(0, 1);

  // session loop and trial loop
  for (n in 1:N) {
    vector[2] Q; // expected value
    real PE;      // prediction error
    vector[Tsesh[n]] Qdiff;

    Q = initQ;

    for (t in 1:(Tsesh[n])) {
      Qdiff[t] = Q[2] - Q[1];
      choice[n, t] ~ bernoulli_logit(beta * Qdiff[t]);

      if (choice[n,t] == 1) {
        PE = outcome[n, t] - Q[2];
        if (PE < 0){
          Q[2] = Q[2] + aN * PE;
        }else{
          Q[2] = Q[2] + aP[n] * PE;
        }
        Q[1] = Q[1] * aF;
      }else{
        PE = outcome[n, t] - Q[1];
        if (PE < 0){
          Q[1] = Q[1] + aN * PE;
        }else{
          Q[1] = Q[1] + aP[n] * PE;
        }
        Q[2] = Q[2] * aF;
      }
    }
  }
}
generated quantities {
  // For group level parameters
  real<lower=0, upper=1> mu_aP;

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

  mu_aP   = Phi_approx(mu_p);

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
        log_lik[n] = log_lik[n] + bernoulli_logit_lpmf(choice[n, t] | beta * Qdiff[t]);

        // generate posterior prediction for current trial
        y_pred[n, t] = categorical_rng(softmax(beta * Q));

        if (choice[n,t] == 1) {
          PE = outcome[n, t] - Q[2];
          if (PE < 0){
            Q[2] = Q[2] + aN * PE;
          }else{
            Q[2] = Q[2] + aP[n] * PE;
          }
          Q[1] = Q[1] * aF;
        }else{
          PE = outcome[n, t] - Q[1];
          if (PE < 0){
            Q[1] = Q[1] + aN * PE;
          }else{
            Q[1] = Q[1] + aP[n] * PE;
          }
          Q[2] = Q[2] * aF;
        }
      }
    }
  }
}
