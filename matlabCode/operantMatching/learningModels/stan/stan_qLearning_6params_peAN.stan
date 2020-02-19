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
  vector[N] aNscale_pr;        // learning rate for NPE
  vector[N] aNmin_pr;        // learning rate for NPE
  vector[N] aP_pr;        // learning rate for PPE
  vector[N] aF_pr;        // forgetting rate
  vector[N] v_pr;  // inverse temperature updating rate
  vector[N] beta_pr;   // inverse temp scale val

}
transformed parameters {
  // session-level parameters
  vector<lower=0, upper=1>[N] aNscale;
  vector<lower=0, upper=1>[N] aNmin;
  vector<lower=0, upper=1>[N] aP;
  vector<lower=0, upper=1>[N] aF;
  vector<lower=0, upper=1>[N] v;
  vector<lower=0, upper=20>[N] beta;

  for (i in 1:N) {
    aNscale[i]   = Phi_approx(mu_p[1]  + sigma[1]  * aNscale_pr[i]);
    aNmin[i]   = Phi_approx(mu_p[2]  + sigma[2]  * aNmin_pr[i]);
    aP[i]   = Phi_approx(mu_p[3]  + sigma[3]  * aP_pr[i]);
    aF[i]   = Phi_approx(mu_p[4]  + sigma[4]  * aF_pr[i]);
    v[i] = Phi_approx(mu_p[5] + sigma[5] * v_pr[i]) ;
    beta[i] = Phi_approx(mu_p[6] + sigma[6] * beta_pr[i]) * 20;
  }
}
model {
  // Hyperparameters
  mu_p  ~ normal(0, 1);
  sigma ~ cauchy(0, 1);

  // individual parameters
  aNscale_pr   ~ normal(0, 1);
  aNmin_pr   ~ normal(0, 1);
  aP_pr   ~ normal(0, 1);
  aF_pr   ~ normal(0, 1);
  v_pr ~ normal(0, 1);
  beta ~ normal(0, 1);

  // session loop and trial loop
  for (i in 1:N) {
    vector[2] Q; // expected value
    real PE;      // prediction error
    real R; // expected average value
    real aN;

    Q = initQ;
    R = 0;
    aN = aNmin[i];

    for (t in 1:(Tsesh[i])) {
      // compute action probabilities
      choice[i, t] ~ categorical_logit(beta[i] * Q);

      // prediction error
      PE = outcome[i, t] - Q[choice[i, t]];

      // value updating (learning)
      if (PE < 0){
        Q[choice[i, t]] += aN * PE;
      }
      else{
        Q[choice[i, t]] += aP[i] * PE;
      }
      if (choice[i, t] == 1){
        Q[2] = Q[2] * aF[i];
      }else{
        Q[1] = Q[1] * aF[i];
      }
      R = v[i] * abs(PE) + (1-v[i]) * R;
      aN = R * aNscale[i] + aNmin[i];
    }
  }
}
generated quantities {
  // For group level parameters
  real<lower=0, upper=1> mu_aNscale;
  real<lower=0, upper=1> mu_aNmin;
  real<lower=0, upper=1> mu_aP;
  real<lower=0, upper=1> mu_aF;
  real<lower=0, upper=20> mu_v;
  real<lower=0, upper=20> mu_beta;

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

  mu_aNscale = Phi_approx(mu_p[1]);
  mu_aNmin   = Phi_approx(mu_p[2]);
  mu_aP      = Phi_approx(mu_p[3]);
  mu_aF      = Phi_approx(mu_p[4]);
  mu_v       = Phi_approx(mu_p[5]); 
  mu_beta    = Phi_approx(mu_p[6]) * 20;

  { // local section, this saves time and space
    for (i in 1:N) {
      vector[2] Q; // expected value
      real PE;      // prediction error
      real R;
      real aN;

      Q = initQ;
      R = 0;
      aN = aNmin[i];

      log_lik[i] = 0;

      for (t in 1:(Tsesh[i])) {
        // compute log likelihood of current trial
        log_lik[i] += categorical_logit_lpmf(choice[i, t] | beta[i] * Q);

        // generate posterior prediction for current trial
        y_pred[i, t] = categorical_rng(softmax(beta[i] * Q));

        // prediction error
        PE = outcome[i, t] - Q[choice[i, t]];

        // value updating (learning)
        if (PE < 0){
          Q[choice[i, t]] += aN * PE;
        }
        else{
          Q[choice[i, t]] += aP[i] * PE;
        }
        if (choice[i, t] == 1){
          Q[2] = Q[2] * aF[i];
        }else{
          Q[1] = Q[1] * aF[i];
        }
        R = v[i] * abs(PE) + (1-v[i]) * R;
        aN = R * aNscale[i] + aNmin[i];
      }
    }
  }
}
