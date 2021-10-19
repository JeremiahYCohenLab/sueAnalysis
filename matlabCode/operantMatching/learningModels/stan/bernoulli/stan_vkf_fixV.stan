data {
  int<lower=1> N;
  int<lower=1> T;
  int<lower=1, upper=T> Tsesh[N];
  int<lower=0, upper=1> choice[N, T];
  int<lower=0, upper=1> outcome[N, T]; 
}
transformed data {
  row_vector[2] initm;  // initial values for m
  initm = rep_row_vector(0.0, 2);
}
parameters {
// Declare all parameters as vectors for vectorizing
  // Hyper(animal)-parameters
  vector[4] mu_p;
  vector<lower=0>[4] sigma;

  // Session-level raw parameters
  vector[N] lambda_pr;    // learning rate for volatility
  vector[N] omega_pr;    // 'noise'
  vector[N] beta_pr;  // inverse temperature
  //vector[N] aF_pr;       //forget rate
}
transformed parameters {
// Transform session-level raw parameters
  vector<lower=0, upper=1>[N] lambda;
  real<lower=0, upper=10> vInit;
  vector<lower=0, upper=10>[N] omega;
  vector<lower=0, upper=1>[N] beta;

  for (n in 1:N) {
    lambda[n]   = inv_logit(mu_p[1] + sigma[1] * lambda_pr[n]);
    omega[n]   = inv_logit(mu_p[3] + sigma[3] * omega_pr[n])*10;
    beta[n] = inv_logit(mu_p[4] + sigma[4] * beta_pr[n]);
  }

  vInit = inv_logit(mu_p[2])*10;
}
model {
  // Hyperparameters
  mu_p  ~ normal(0, 1);
  sigma[1] ~ cauchy(0, 2);
  sigma[2] ~ cauchy(0, 2);
  sigma[3] ~ cauchy(0, 1);
  sigma[4] ~ cauchy(0, 2);

  // individual parameters
  lambda_pr   ~ normal(0, 1);
  omega_pr   ~ normal(0, 1);
  beta_pr ~ normal(0, 1);

  // session loop and trial loop
  for (n in 1:N) {
    matrix[Tsesh[n]+1, 2] m; // expected value
    matrix[Tsesh[n]+1, 2] w; // estimation uncertainty
    matrix[Tsesh[n]+1, 2] v; // estimated volatility
    real wcov; //estimation covariance
    real k; //kalman filter
    real alpha; //learning rate
    real PE;      // prediction error
    vector[Tsesh[n]] mdiff;

    m[1] = initm;
    w[1] = rep_row_vector(omega[n], 2);
    v[1] = rep_row_vector(vInit, 2);

    for (t in 1:(Tsesh[n])) {
      mdiff[t] = m[t,2] - m[t,1];
      //print(mdiff[t]);
      choice[n, t] ~ bernoulli_logit(beta[n] * mdiff[t]);

      if (choice[n,t] == 1) {
        PE = outcome[n, t] - inv_logit(m[t,2]);
        k = (w[t,2] + v[t,2])/(w[t,2] + v[t,2] + omega[n]);
        alpha = sqrt(w[t,2] + v[t,2]);
        m[t+1,2] = m[t,2] + alpha * PE;
        w[t+1,2] = (1 - k) * (w[t,2] + v[t,2]);
        wcov = (1 - k) * w[t, 2];
        v[t+1, 2] = v[t,2] + lambda[n] * ( (m[t+1, 2] - m[t, 2])^2 + w[t, 2] + w[t+1, 2] - 2 * wcov - v[t, 2] );

        m[t+1, 1] = m[t, 1]; 
        w[t+1, 1] = w[t, 1]; 
        v[t+1, 1] = v[t, 1]; 
      }else{
        PE = outcome[n, t] - inv_logit(m[t,1]);
        k = (w[t,1] + v[t,1])/(w[t,1] + v[t,1] + omega[n]);
        alpha = sqrt(w[t,1] + v[t,1]);
        m[t+1,1] = m[t,1] + alpha * PE;
        w[t+1,1] = (1 - k) * (w[t,1] + v[t,1]);
        wcov = (1 - k) * w[t, 1];
        v[t+1, 1] = v[t,1] + lambda[n] * ( (m[t+1, 1] - m[t, 1])^2 + w[t, 1] + w[t+1, 1] - 2 * wcov - v[t, 1] );

        m[t+1, 2] = m[t, 2]; 
        w[t+1, 2] = w[t, 2]; 
        v[t+1, 2] = v[t, 2]; 

      // print("mdiff = ", mdiff[t]);
      // print("PE =", PE);
      }
    }
  }
}


generated quantities {
  // For animal level parameters

  real<lower=0, upper=1> mu_lambda;
  real<lower=0, upper=10> mu_vInit;
  real<lower=0, upper=10> mu_omega;
  real<lower=0, upper=1> mu_beta;

  // For log likelihood calculation
  real log_lik[N];

  mu_lambda  = inv_logit(mu_p[1]);
  mu_vInit   = inv_logit(mu_p[2])*10;
  mu_omega   = inv_logit(mu_p[3])*10;
  mu_beta    = inv_logit(mu_p[4]);
  

  { // local section, this saves time and space
    for (n in 1:N) {
      matrix[Tsesh[n]+1, 2] m; // expected value
      matrix[Tsesh[n]+1, 2] w; // estimation uncertainty
      matrix[Tsesh[n]+1, 2] v; // estimated volatility
      real wcov; //estimation covariance
      real k; //kalman filter
      real alpha; //learning rate
      real PE;      // prediction error
      vector[Tsesh[n]] mdiff;

      m[1] = initm;
      w[1] = rep_row_vector(omega[n], 2);
      v[1] = rep_row_vector(vInit, 2);

      log_lik[n] = 0;

      for (t in 1:(Tsesh[n])) {
        mdiff[t] = m[t,2] - m[t,1];
        log_lik[n] += bernoulli_logit_lpmf(choice[n, t] | beta[n] * mdiff[t]);

        if (choice[n,t] == 1) {
          PE = outcome[n, t] - inv_logit(m[t,2]);
          k = (w[t,2] + v[t,2])/(w[t,2] + v[t,2] + omega[n]);
          alpha = sqrt(w[t,2] + v[t,2]);
          m[t+1,2] = m[t,2] + alpha * PE;
          w[t+1,2] = (1 - k) * (w[t,2] + v[t,2]);
          wcov = (1 - k) * w[t, 2];
          v[t+1, 2] = v[t,2] + lambda[n] * ( (m[t+1, 2] - m[t, 2])^2 + w[t, 2] + w[t+1, 2] - 2 * wcov - v[t, 2] );

          m[t+1, 1] = m[t, 1]; 
          w[t+1, 1] = w[t, 1]; 
          v[t+1, 1] = v[t, 1]; 
        }else{
          PE = outcome[n, t] - inv_logit(m[t,1]);
          k = (w[t,1] + v[t,1])/(w[t,1] + v[t,1] + omega[n]);
          alpha = sqrt(w[t,1] + v[t,1]);
          m[t+1,1] = m[t,1] + alpha * PE;
          w[t+1,1] = (1 - k) * (w[t,1] + v[t,1]);
          wcov = (1 - k) * w[t, 1];
          v[t+1, 1] = v[t,1] + lambda[n] * ( (m[t+1, 1] - m[t, 1])^2 + w[t, 1] + w[t+1, 1] - 2 * wcov - v[t, 1] );

          m[t+1, 2] = m[t, 2]; 
          w[t+1, 2] = w[t, 2]; 
          v[t+1, 2] = v[t, 2]; 
        }
      }


    }
  }
}