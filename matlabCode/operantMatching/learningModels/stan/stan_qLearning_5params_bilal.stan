data {
  int<lower=1> N;
  int<lower=1> T;
  int<lower=1, upper=T> Tsesh[N];
  int<lower=0, upper=2> choice[N, T];
  real outcome[N, T];  // no lower and upper bounds
}
transformed data {
  vector[2] initQ;  // initial values for Q
  vector[2] initCP;
  initQ = rep_vector(0.0, 2);
  initCP = rep_vector(0.0, 2);
}
parameters {
// Declare all parameters as vectors for vectorizing
  // Hyper(animal)-parameters
  vector[4] mu_p;
  vector<lower=0>[4] sigma;

  // Session-level raw parameters
  vector[N] alpha_pr; // learning rate 
  vector[N] aF_pr;    // forgetting rate
  vector[N] beta_pr;  // inverse temperature
  vector[N] k_pr;     // choice autocorrelation
  vector[N] bias_pr;  // bias
}
transformed parameters {
// Transform session-level raw parameters
  vector<lower=0, upper=1>[N] alpha;
  vector<lower=0, upper=1>[N] aF;
  vector<lower=0, upper=40>[N] beta;
  vector<lower=-5, upper=5>[N] k;
  vector<lower=-5, upper=5>[N] bias;

  for (n in 1:N) {
    alpha[n] = Phi_approx(mu_p[1] + sigma[1] * alpha_pr[n]);
    aF[n]    = Phi_approx(mu_p[2] + sigma[2] * aF_pr[n]);
    beta[n]  = Phi_approx(mu_p[3] + sigma[3] * beta_pr[n]) * 40;
    k[n]     = Phi_approx(mu_p[4] + sigma[4] * k_pr[n]) * 10 - 5;
    bias[n]  = Phi_approx(bias_pr[n]) * 10 - 5;
  }
}
model {
  // Hyperparameters
  mu_p  ~ normal(0, 0.5);
  sigma ~ cauchy(0, 0.5);

  // individual parameters
  alpha_pr ~ normal(0, 1);
  aF_pr    ~ normal(0, 1);
  beta_pr  ~ normal(0, 0.0001);
  k_pr     ~ normal(0, 1);
  bias_pr  ~ normal(0,1);

  // session loop and trial loop
  for (n in 1:N) {
    vector[2] Q; // expected value
    vector[2] cP;
    real PE;      // prediction error
    real prevChoice;

    prevChoice = 0;
    Q = initQ;
    cP = initCP;

    for (t in 1:(Tsesh[n])) {
      // compute action probabilities
      cP[2] = 1 / (1 + exp(-(beta[n] * (Q[2] - Q[1]) + k[n] * prevChoice +  bias[n])));
      cP[1] = 1 - cP[2];
      choice[n, t] ~ categorical(cP);

      // prediction error
      PE = outcome[n, t] - Q[choice[n, t]];
      // value updating (learning)
      Q[choice[n, t]] = Q[choice[n, t]] * aF[n] + alpha[n] * PE;
      if (choice[n, t] == 2){
        Q[1] = Q[1] * aF[n];
        prevChoice = 1;
      }else{
        Q[2] = Q[2] * aF[n];
        prevChoice = -1;
      }
    }
  }
}
generated quantities {
  // For group level parameters
  real<lower=0, upper=1> mu_alpha;
  real<lower=0, upper=1> mu_aF;
  real<lower=0, upper=40> mu_beta;
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

  mu_alpha = Phi_approx(mu_p[1]);
  mu_aF    = Phi_approx(mu_p[2]);
  mu_beta  = Phi_approx(mu_p[3]) * 40;
  mu_k     = Phi_approx(mu_p[4]) * 10 - 5;
  
  { // local section, this saves time and space
    for (n in 1:N) {
      vector[2] Q; // expected value
      vector[2] cP;
      real PE;      // prediction error
      real prevChoice;

      // Initialize values
      prevChoice = 0;
      Q = initQ;
      cP = initCP;
      log_lik[n] = 0;

      for (t in 1:(Tsesh[n])) {
        // compute log likelihood of current trial
        cP[2] = 1 / (1 + exp(-(beta[n] * (Q[2] - Q[1]) + k[n] * prevChoice + bias[n])));
        cP[1] = 1 - cP[2];
        log_lik[n] += categorical_lpmf(choice[n, t] | cP);

        // generate posterior prediction for current trial
        y_pred[n, t] = categorical_rng(cP);

        // prediction error
        PE = outcome[n, t] - Q[choice[n, t]];
        // value updating (learning)
        Q[choice[n, t]] = Q[choice[n, t]] * aF[n] + alpha[n] * PE;
        if (choice[n, t] == 2){
          Q[1] = Q[1] * aF[n];
          prevChoice = 1;
        }else{
          Q[2] = Q[2] * aF[n];
          prevChoice = -1;
        }
      }
    }
  }
}
