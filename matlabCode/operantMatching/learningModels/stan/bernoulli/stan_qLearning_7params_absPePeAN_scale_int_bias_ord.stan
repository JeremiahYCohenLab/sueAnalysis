data {
  int<lower=1> N;                      // number of sessions
  int<lower=1> T;                      // max number of trials per session
  int<lower=1, upper=T> Tsesh[N];      // number of trials in each session
  int<lower=0, upper=2> choice[N, T];  // matrix of choices (0 for left, 1 for right)
  int<lower=0, upper=1> outcome[N, T]; // matrix of outcomes (0 for no reward, 1 for reward)
}
transformed data {
  vector[2] initQ;  // initial values for Q
  initQ = rep_vector(0.0, 2);
}
parameters {
// declare all parameters as vectors for vectorizing
  // hyper(animal)-parameters
  vector[6] mu_p;
  vector<lower=0>[4] sigma;

  // session-level raw parameters
  vector[N] aNmin_pr;   // median value for alphaNPE (negative reward prediction error learning rate)
  vector[N] aP_pr;      // alphaPPE (positive reward prediction error learning rate)
  vector[N] aF_pr;      // forgetting rate
  vector[N] v_pr;       // integration rate for alphaNPE
  vector[N] aPE_pr;     // updating rate for expected uncertainty
  vector[N] beta_pr;    // inverse temp for softmax decision function
  vector[N] bias;       // bias term

}
transformed parameters {
  // session-level parameters
  vector<lower=0, upper=1>[N] aNmin;
  vector<lower=0, upper=1>[N] aP;
  vector<lower=0, upper=1>[N] aF;
  vector<lower=0, upper=1>[N] aPE;
  vector<lower=0, upper=1>[N] v;
  vector<lower=0, upper=10>[N] beta;

  for (n in 1:N) {    // transform parameters to [0 1] range or [0 1*x]
    aNmin[n]   = Phi_approx(mu_p[1] + sigma[1] * aNmin_pr[n]);
    aP[n]      = Phi_approx(mu_p[2] + sigma[2] * aP_pr[n]);
    aF[n]      = Phi_approx(mu_p[3] + sigma[3] * aF_pr[n]);
    beta[n]    = Phi_approx(mu_p[4] + sigma[4] * beta_pr[n]) * 10;
    v[n]       = inv_logit(mu_p[5] + v_pr[n]);
    aPE[n]     = v[n] * inv_logit(mu_p[6] + aPE_pr[n]);
  }
}
model {
  // hyperparameters
  mu_p  ~ normal(0, 1);
  sigma ~ cauchy(0, 3);

  // individual parameters
  aNmin_pr ~ normal(0, 1);
  aP_pr    ~ normal(0, 1);
  aF_pr    ~ normal(0, 1);
  beta_pr  ~ normal(0, 1);
  bias     ~ normal(0, 10);
  v_pr     ~ normal(0, 5);
  aPE_pr   ~ normal(0, 5);

  // session loop and trial loop
  for (n in 1:N) {
    vector[2] Q;            // expected value for each action
    real PE;                // prediction error
    real pePE;              // prediction error prediction error (surprise)
    real peBar;             // expected uncertainty
    real aN;                // NPE learning rate
    vector[Tsesh[n]] Qdiff; // relative value of actions (for softmax decision function)

    Q = initQ;
    peBar = 0;
    aN = aNmin[n];

    for (t in 1:(Tsesh[n])) {
      Qdiff[t] = Q[2] - Q[1];
      choice[n, t] ~ bernoulli_logit(beta[n] * Qdiff[t] + bias[n]);
      if (choice[n,t] == 1) {
        PE = outcome[n, t] - Q[2];
        pePE = fabs(PE) - peBar;
        if (PE < 0){
          aN = v[n] * (pePE + aNmin[n]) + (1 - v[n]) * aN;
          if (aN < 0){
            aN = 0;
          }
          if (aN > 1){
            aN = 1;
          }
          Q[2] += aN * PE * (1-peBar);
        }else{
          Q[2] += aP[n] * PE * (1-peBar);
        }
        Q[1] = Q[1] * aF[n];
      }else{
        PE = outcome[n, t] - Q[1];
        pePE = fabs(PE) - peBar;
        if (PE < 0){
          aN = v[n] * (pePE + aNmin[n]) + (1 - v[n]) * aN;
          if (aN < 0){
            aN = 0;
          }
          if (aN > 1){
            aN = 1;
          }
          Q[1] += aN * PE * (1-peBar);
        }else{
          Q[1] += aP[n] * PE * (1-peBar);
        }
        Q[2] = Q[2] * aF[n];
      }
      peBar += aPE[n] * pePE;
      //print("Qdiff = ", Qdiff[t]);
      //print("PE =", PE);
      //print("pePE = ", pePE);
    }
  }
}
generated quantities {
  // For group level parameters
  real<lower=0, upper=1> mu_aNmin;
  real<lower=0, upper=1> mu_aP;
  real<lower=0, upper=1> mu_aF;
  real<lower=0, upper=1> mu_aPE;
  real<lower=0, upper=1> mu_v;
  real<lower=0, upper=10> mu_beta;

  // For log likelihood calculation
  real log_lik[N];

  // For posterior predictive check
//  real y_pred[N, T];

  // Set all posterior predictions to 0 (avoids NULL values)
//  for (n in 1:N) {
//    for (t in 1:T) {
//      y_pred[n, t] = -1;
//    }
//  }

  mu_aNmin   = Phi_approx(mu_p[1]);
  mu_aP      = Phi_approx(mu_p[2]);
  mu_aF      = Phi_approx(mu_p[3]);
  mu_beta    = Phi_approx(mu_p[4]) * 10;
  mu_v       = inv_logit(mu_p[5]);
  mu_aPE     = mu_v * inv_logit(mu_p[6]);
  

  { // local section, this saves time and space
    for (n in 1:N) {
      vector[2] Q;            // expected value for each action
      real PE;                // prediction error
      real pePE;              // prediction error prediction error (surprise)
      real peBar;             // expected average value
      real aN;                // NPE learning rate
      vector[Tsesh[n]] Qdiff; // relative value of actions (for softmax decision function)
  //    vector[2] cP;

      Q = initQ;
  //    cP = initQ;
      peBar = 0;
      aN = aNmin[n];

      log_lik[n] = 0;

      for (t in 1:(Tsesh[n])) {
        Qdiff[t] = Q[2] - Q[1];

        // compute log likelihood of current trial
        log_lik[n] += bernoulli_logit_lpmf(choice[n, t] | beta[n] * Qdiff[t] + bias[n]);

        // generate posterior prediction for current trial
    //    cP[2] = 1 / (1 + exp(-(beta[n] * Qdiff[t] + bias[n])));
    //    cP[1] = 1 - cP[2];
    //    y_pred[n, t] = categorical_rng(cP);

        if (choice[n,t] == 1) {
          PE = outcome[n, t] - Q[2];
          pePE = fabs(PE) - peBar;
          if (PE < 0){
            aN = v[n] * (pePE + aNmin[n]) + (1 - v[n]) * aN;
            if (aN < 0){
              aN = 0;
            }
            if (aN > 1){
              aN = 1;
            }
            Q[2] += aN * PE * (1-peBar);
          }else{
            Q[2] += aP[n] * PE * (1-peBar);
          }
          Q[1] = Q[1] * aF[n];
        }else{
          PE = outcome[n, t] - Q[1];
          pePE = fabs(PE) - peBar;
          if (PE < 0){
            aN = v[n] * (pePE + aNmin[n]) + (1 - v[n]) * aN;
            if (aN < 0){
              aN = 0;
            }
          if (aN > 1){
            aN = 1;
          }
            Q[1] += aN * PE * (1-peBar);
          }else{
            Q[1] += aP[n] * PE * (1-peBar);
          }
          Q[2] = Q[2] * aF[n];
        }
        peBar += aPE[n] * pePE;
      }
    }
  }
}
