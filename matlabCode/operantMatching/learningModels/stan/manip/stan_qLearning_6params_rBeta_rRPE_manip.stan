data {
  int<lower=1> N;
  int<lower=1> M;
  int<lower=1> T;
  int<lower=1, upper=T> Tsesh[N];
  int<lower=0, upper=2> choice[N, T];
  int<lower=0, upper=1> outcome[N, T];
  int<lower=1, upper=T> TseshM[M];
  int<lower=0, upper=2> choiceM[M, T];
  int<lower=0, upper=1> outcomeM[M, T];
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
  vector[6] d_mu_p;
  vector<lower=0>[6] d_sigma;

  // Session-level raw parameters
  vector[N] aN_pr;    // learning rate for NPE
  vector[N] aP_pr;    // learning rate for PPE
  vector[N] aF_pr;    // forgetting rate
  vector[N] rRate_pr;     // expected average value learning rate
  vector[N] betaMin_pr;  // inverse temperature min value
  vector[N] betaDiff_pr;  // inverse temperature max-min value
  

  vector[M] d_aN_pr;    // change in learning rate for NPE
  vector[M] d_aP_pr;    // change in learning rate for PPE
  vector[M] d_aF_pr;    // change in forgetting rate
  vector[M] d_rRate_pr;     // change in expected average value learning rate
  vector[N] d_betaMin_pr;  // change in inverse temperature min value
  vector[N] d_betaDiff_pr;  // change in inverse temperature max-min value
}
transformed parameters {
  // session-level parameters
  vector<lower=0, upper=1>[N] aN;
  vector<lower=0, upper=1>[N] aP;
  vector<lower=0, upper=1>[N] aF;
  vector<lower=0, upper=1>[N] rRate;
  vector<lower=0, upper=20>[N] betaMin;
  vector<lower=0, upper=20>[N] betaDiff;

  vector<lower=-1, upper=1>[M] d_aN;
  vector<lower=-1, upper=1>[M] d_aP;
  vector<lower=-1, upper=1>[M] d_aF;
  vector<lower=-1, upper=1>[M] d_rRate;
  vector<lower=-10, upper=10>[N] d_betaMin;
  vector<lower=-10, upper=10>[N] d_betaDiff;  


  for (i in 1:N) {
    aN[i]   = Phi_approx(mu_p[1]  + sigma[1]  * aN_pr[i]);
    aP[i]   = Phi_approx(mu_p[2]  + sigma[2]  * aP_pr[i]);
    aF[i]   = Phi_approx(mu_p[3]  + sigma[3]  * aF_pr[i]);
    rRate[i]   = Phi_approx(mu_p[4]  + sigma[4]  * rRate_pr[i]);
    betaMin[i] = Phi_approx(mu_p[5] + sigma[5] * betaMin_pr[i]) * 10;
    betaDiff[i] = Phi_approx(mu_p[6] + sigma[6] * betaDiff_pr[i]) * 10;
    }
  for (i in 1:M){
    d_aN[i]   = -1 + 2*Phi_approx(d_mu_p[1]  + d_sigma[1]  * d_aN_pr[i]);
    d_aP[i]   = -1 + 2*Phi_approx(d_mu_p[2]  + d_sigma[2]  * d_aP_pr[i]);
    d_aF[i]   = -1 + 2*Phi_approx(d_mu_p[3]  + d_sigma[3]  * d_aF_pr[i]);
    d_rRate[i]    = -1 + 2*Phi_approx(d_mu_p[4] + d_sigma[4] * d_rRate_pr[i]);
    d_betaMin[i] = -10 + 20*Phi_approx(d_mu_p[5]  + d_sigma[5]  * d_betaMin_pr[i]);
    d_betaDiff[i] = -10 + 20*Phi_approx(d_mu_p[6]  + d_sigma[6]  * d_betaDiff_pr[i]);    
  }
}
model {
  // Hyperparameters
  mu_p  ~ normal(0, 1);
  sigma ~ cauchy(0, 1);
  d_mu_p ~ normal(0,1);
  d_sigma ~ cauchy(0,1);

  // individual parameters
  aN_pr   ~ normal(0, 1);
  aP_pr   ~ normal(0, 1);
  aF_pr   ~ normal(0, 1);
  rRate_pr    ~ normal(0, 1);
  betaMin_pr ~ normal(0, 1);
  betaDiff_pr ~ normal(0, 1);    
  d_aN_pr   ~ normal(0, 1);
  d_aP_pr   ~ normal(0, 1);
  d_aF_pr   ~ normal(0, 1);
  d_rRate_pr    ~ normal(0, 1);
  d_betaMin_pr ~ normal(0, 1);
  d_betaDiff_pr ~ normal(0, 1);



  // session loop and trial loop
  for (i in 1:N) {
    vector[2] Q; // expected value
    real PE;      // prediction error
    real R; // expected average value
    real beta; //inverse temperature

    Q = initQ;
    R = 0;
    beta = 0;

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
      R = R + rRate[i] * (outcome[i, t] - R);
      if (R < 0){
        R = 0;
      }
      beta = betaMin[i] + betaDiff[i] * R;
    }
  }
  for (i in 1:M) {
    vector[2] Q; // expected value
    real PE;      // prediction error
    real R; // expected average value
    real beta; // inverse temperature

    Q = initQ;
    R = 0.4;
    beta = 0;

    for (t in 1:(TseshM[i])) {
      // compute action probabilities
      choiceM[i, t] ~ categorical_logit(beta * Q);

      // prediction error
      PE = outcomeM[i, t] - Q[choiceM[i, t]];

      // value updating (learning)
      if (PE < 0){
        Q[choiceM[i, t]] = Q[choiceM[i, t]] + (aN[i] + d_aN[i]) * PE;
      }
      else{
        Q[choiceM[i, t]] = Q[choiceM[i, t]] + (aP[i] + d_aP[i]) * PE;
      }
      if (choiceM[i, t] == 1){
        Q[2] = Q[2] * (aF[i] + d_aF[i]);
      }else{
        Q[1] = Q[1] * (aF[i] + d_aF[i]);
      }
      R = R + (rRate[i] + d_rRate[i]) * (outcomeM[i, t] - R);
      if (R < 0){
        R = 0;
      }
      beta = betaMin[i] + d_betaMin[i] + (betaDiff[i] + d_betaDiff[i]) * R;
    }
  }
}
generated quantities {
  // For group level parameters
  real<lower=0, upper=1> mu_aN;
  real<lower=0, upper=1> mu_aP;
  real<lower=0, upper=1> mu_aF;
  real<lower=0, upper=10> mu_beta;
  real<lower=0, upper=1> mu_rRate;
  real<lower=-1, upper=1> d_mu_aN;
  real<lower=-1, upper=1> d_mu_aP;
  real<lower=-1, upper=1> d_mu_aF;
  real<lower=-10, upper=10> d_mu_beta;
  real<lower=-1, upper=1> d_mu_rRate;

  // For log likelihood calculation
  real log_lik[N];
  real log_likM[M];

  // For posterior predictive check
  real y_pred[N, T];
  real y_predM[M,T];

  // Set all posterior predictions to 0 (avoids NULL values)
  for (i in 1:N) {
    for (t in 1:T) {
      y_pred[i, t] = -1;
    }
  }
  for (i in 1:M) {
    for (t in 1:T) {
      y_predM[i, t] = -1;
    }
  }


  mu_aN   = Phi_approx(mu_p[1]);
  mu_aP   = Phi_approx(mu_p[2]);
  mu_aF   = Phi_approx(mu_p[3]);
  mu_beta = Phi_approx(mu_p[4]) * 10;
  mu_rRate  = Phi_approx(mu_p[5]);

  d_mu_aN   = -1 + 2*(Phi_approx(d_mu_p[1]));
  d_mu_aP   = -1 + 2*(Phi_approx(d_mu_p[2]));
  d_mu_aF   = -1 + 2*(Phi_approx(d_mu_p[3]));
  d_mu_beta = -10 + 20*(Phi_approx(d_mu_p[4]));
  d_mu_rRate   = -1 + 2*(Phi_approx(d_mu_p[5]));

  { // local section, this saves time and space
    for (i in 1:N) {
      vector[2] Q; // expected value
      real PE;      // prediction error
      real R; // expected average value
      real beta; // inverse temperature

      // Initialize values
      Q = initQ;
      R = 0;
      beta = 0;

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
        R = R + rRate[i] * (outcome[i, t] - R);
        if (R < 0){
          R = 0;
        }
        beta = betaMin[i] + betaDiff[i] * R;
      }
    }
    for (i in 1:M) {
      vector[2] Q; // expected value
      real PE;      // prediction error
      real R; // expected average value
      real beta; // inverse temperature

      // Initialize values
      Q = initQ;
      R = 0;
      beta = 0;

      log_likM[i] = 0;

      for (t in 1:(TseshM[i])) {
        // compute log likelihood of current trial
        log_likM[i] = log_likM[i] + categorical_logit_lpmf(choiceM[i, t] | beta * Q);

        // generate posterior prediction for current trial
        y_predM[i, t] = categorical_rng(softmax(beta * Q));

        // prediction error
        PE = outcomeM[i, t] - Q[choiceM[i, t]];

        // value updating (learning)
        if (PE < 0){
          Q[choiceM[i, t]] = Q[choiceM[i, t]] + (aN[i] + d_aN[i]) * PE;
        }
        else{
          Q[choiceM[i, t]] = Q[choiceM[i, t]] + (aP[i] + d_aP[i]) * PE;
        }
        if (choiceM[i, t] == 1){
          Q[2] = Q[2] * (aF[i] + d_aF[i]);
        }else{
          Q[1] = Q[1] * (aF[i] + d_aF[i]);
        }
        R = R + (rRate[i] + d_rRate[i]) * (outcomeM[i, t] - R);
        if (R < 0){
          R = 0;
        }
        beta = betaMin[i] + d_betaMin[i] + (betaDiff[i] + d_betaDiff[i]) * R;
      }
    }
  }
}
