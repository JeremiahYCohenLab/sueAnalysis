data {
  int<lower=1> N;
  int<lower=1> T;
  int<lower=1, upper=T> Tsesh[N];
  int<lower=0, upper=2> choice[N, T];
  int<lower=0, upper=1> outcome[N, T];
}
transformed data {
  vector[2] initQ;  // initial values for Q
  vector[2] initCP;  // initial values for choice probability
  initQ  = rep_vector(0.0, 2);
  initCP = rep_vector(0.0, 2);
}
parameters {
// Declare all parameters as vectors for vectorizing
  // Hyper(animal)-parameters
  vector[9] mu_p;
  vector<lower=0>[9] sigma;

  // Session-level raw parameters
  vector[N] aNscale_pr;    // learning rate for NPE
  vector[N] aNmin_pr;      // learning rate for NPE
  vector[N] aP_pr;         // learning rate for PPE
  vector[N] aF_pr;         // forgetting rate
  vector[N] aPE_pr;        // learning rate for volatility
  vector[N] v_pr;          // inverse temperature updating rate
  vector[N] betaScale_pr;  // inverse temp scale val
  vector[N] betaMin_pr;    // inverse temp scale val
  vector[N] k_pr;          // choice autocorrelation
  vector[N] b_pr;          // choice bias

}
transformed parameters {
  // session-level parameters
  vector<lower=0, upper=1>[N] aNscale;
  vector<lower=0, upper=1>[N] aNmin;
  vector<lower=0, upper=1>[N] aP;
  vector<lower=0, upper=1>[N] aF;
  vector<lower=0, upper=1>[N] aPE;
  vector<lower=0, upper=1>[N] v;
  vector<lower=0, upper=20>[N] betaScale;
  vector<lower=0, upper=20>[N] betaMin;
  vector<lower=-5, upper=5>[N] k;
  vector<lower=-5, upper=5>[N] b;

  for (i in 1:N) {
    aNscale[i]    = Phi_approx(mu_p[1]  + sigma[1]  * aNscale_pr[i]);
    aNmin[i]      = Phi_approx(mu_p[2]  + sigma[2]  * aNmin_pr[i]);
    aP[i]         = Phi_approx(mu_p[3]  + sigma[3]  * aP_pr[i]);
    aF[i]         = Phi_approx(mu_p[4]  + sigma[4]  * aF_pr[i]);
    aPE[i]        = Phi_approx(mu_p[5]  + sigma[5]  * aPE_pr[i]) ;
    v[i]          = Phi_approx(mu_p[6]  + sigma[6]  * v_pr[i]) ;
    betaScale[i]  = Phi_approx(mu_p[7]  + sigma[7]  * betaScale_pr[i]) * 20;
    betaMin[i]    = Phi_approx(mu_p[8]  + sigma[8]  * betaMin_pr[i]) * 20;
    k[i]          = Phi_approx(mu_p[9]  + sigma[9]  * k_pr[i]) * 10 - 5;
    b[i]          = Phi_approx(b_pr[i]) * 10 - 5;
  }
}
model {
  // Hyperparameters
  mu_p  ~ normal(0, 1);
  sigma ~ cauchy(0, 1);

  // individual parameters
  aNscale_pr    ~ normal(0, 1);
  aNmin_pr      ~ normal(0, 1);
  aP_pr         ~ normal(0, 1);
  aF_pr         ~ normal(0, 1);
  aPE_pr        ~ normal(0, 1);
  v_pr          ~ normal(0, 1);
  betaScale_pr  ~ normal(0, 1);
  betaMin_pr    ~ normal(0, 1);
  k_pr          ~ normal(0, 1);
  b_pr          ~ normal(0, 1);

  // session loop and trial loop
  for (i in 1:N) {
    vector[2] Q;        // expected value
    vector[2] cP;       // choice probability
    real prevChoice;
    real PE;            // prediction error
    real peBar;         // expected average value
    real aN;
    real R;
    real beta; 
    vector[Tsesh[i]] Qdiff;

    prevChoice = 0;
    Q = initQ;
    cP = initCP;
    peBar = 0;
    aN = aNmin[i];
    R = 0;
    beta = betaMin[i];

    for (t in 1:(Tsesh[i])) {
      cP[1] = 1 / (1 + exp(-beta * ((Q[1] - Q[2]) + k[i] * prevChoice) + b[i]));
      cP[2] = 1 - cP[1];
      choice[i, t] ~ bernoulli(cP);

      if (choice[i,t] == 1) {
        PE = outcome[i, t] - Q[2];
        peBar += aPE[i] * (abs(PE) - peBar);
        aN = peBar * aNscale[i] + aNmin[i];
        if (PE < 0){
          Q[2] += aN * PE;
        }else{
          Q[2] += aP[i] * PE;
        }
        Q[1] = Q[1] * aF[i];
        prevChoice = 1;
      }else{
        PE = outcome[i, t] - Q[1];
        peBar += aPE[i] * (abs(PE) - peBar);
        aN = peBar * aNscale[i] + aNmin[i];
        if (PE < 0){
          Q[1] += aN * PE;
        }else{
          Q[1] += aP[i] * PE;
        }
        Q[2] = Q[2] * aF[i];
        prevChoice = -1;
      }
      R = v[i] * outcome[i, t] + (1-v[i]) * R;
      beta = R * betaScale[i] + betaMin[i];
    }
  }
}
generated quantities {
  // For group level parameters
  real<lower=0, upper=1>  mu_aNscale;
  real<lower=0, upper=1>  mu_aNmin;
  real<lower=0, upper=1>  mu_aP;
  real<lower=0, upper=1>  mu_aF;
  real<lower=0, upper=1>  mu_aPE;
  real<lower=0, upper=1>  mu_v;
  real<lower=0, upper=20> mu_betaScale;
  real<lower=0, upper=20> mu_betaMin;
  real<lower=-5, upper=5>  mu_k;
  real<lower=-5, upper=5>  mu_b;

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
  mu_aNmin     = Phi_approx(mu_p[2]);
  mu_aP        = Phi_approx(mu_p[3]);
  mu_aF        = Phi_approx(mu_p[4]);
  mu_aPE       = Phi_approx(mu_p[5]);
  mu_v         = Phi_approx(mu_p[6]); 
  mu_betaScale = Phi_approx(mu_p[7]) * 20;
  mu_betaMin   = Phi_approx(mu_p[8]) * 20;
  mu_k         = Phi_approx(mu_p[9]);
  mu_b         = Phi_approx(b_pr) * 10 - 5; 

  { // local section, this saves time and space
    for (i in 1:N) {
      vector[2] Q;        // expected value
      vector[2] cP;       // choice probability
      real prevChoice;
      real PE;            // prediction error
      real peBar;         // expected average value
      real aN;
      real R;
      real beta; 
      vector[Tsesh[i]] Qdiff;

      prevChoice = 0;
      Q = initQ;
      cP = initCP;
      peBar = 0;
      aN = aNmin[i];
      R = 0;
      beta = betaMin[i];

      log_lik[i] = 0;

      for (t in 1:(Tsesh[i])) {
        Qdiff[t] = Q[2] - Q[1];
        // compute log likelihood of current trial
        cP[1] = 1 / (1 + exp(-beta * ((Q[1] - Q[2]) + k[i] * prevChoice) + b[i]));
        cP[2] = 1 - cP[1];
        log_lik[i] += bernoulli_lpmf(choice[i, t] | cP);

        // generate posterior prediction for current trial
        y_pred[i, t] = categorical_rng(cP);

        if (choice[i,t] == 1) {
          PE = outcome[i, t] - Q[2];
          peBar += aPE[i] * (abs(PE) - peBar);
          aN = peBar * aNscale[i] + aNmin[i];
          if (PE < 0){
            Q[2] += aN * PE;
          }else{
            Q[2] += aP[i] * PE;
          }
          Q[1] = Q[1] * aF[i];
          prevChoice = -1;
        }else{
          PE = outcome[i, t] - Q[1];
          peBar += aPE[i] * (abs(PE) - peBar);
          aN = peBar * aNscale[i] + aNmin[i];
          if (PE < 0){
            Q[1] += aN * PE;
          }else{
            Q[1] += aP[i] * PE;
          }
          Q[2] = Q[2] * aF[i];
          prevChoice = -1;
        }
        R = v[i] * outcome[i, t] + (1-v[i]) * R;
        beta = R * betaScale[i] + betaMin[i];
      }
    }
  }
}
