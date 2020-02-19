data {
  int<lower=1>  N; // number of sessions
  int<lower=1>  T; // longest session
  int<lower=1, upper=T> Tsesh[N]; // vector of trial count for each _session_
  int<lower=0>  choice_R[N, T]; // choice data; 0 for left, 1 for right
  real<lower=0> outcome[N, T]; // 0 for no reward, 1 for reward
}
parameters {
  // hierarchical parameters
  real<lower=0> a_alpha; // beta distributed
  real<lower=0> b_alpha;
  real<lower=0> a_aF;
  real<lower=0> b_aF;
  real<lower=0> a_beta; // gamma distributed
  real<lower=0> b_beta;

  // session-level raw parameters
  vector<lower=0,upper=1>[N]  alpha;
  vector<lower=0,upper=1>[N]  aF;
  vector<lower=0,upper=30>[N] betaTemp; // inverse temp
  vector<lower=-5,upper=5>[N] bias;
  vector<lower=-5,upper=5>[N] k;
}
model {
  // priors over hyperparameters
  a_alpha ~ cauchy(0, 5);
  b_alpha ~ cauchy(0, 5);
  a_aF    ~ cauchy(0, 5);
  b_aF    ~ cauchy(0, 5);
  a_beta  ~ cauchy(0, 5);
  b_beta  ~ cauchy(0, 5);

  // distribution for parameters
  alpha ~    beta(a_alpha, b_alpha);
  aF ~   beta(a_aF, b_aF);
  betaTemp ~      gamma(a_beta, 1/b_beta);

  // session loop and trial loop
  for (n in 1:N) {
    vector[2]        Q; // Q-value; [Left, Right]
    real             RPE; // reward prediction error; trials - 1
    vector[Tsesh[n]] Qdiff;
    vector[Tsesh[n] + 1] prevChoice;

    Q[1] = 0;
    Q[2] = 0;
    prevChoice[1] = 0;

    for (t in 1:Tsesh[n]) {
      Qdiff[t] = Q[2] - Q[1];
      if (choice_R[n, t] == 1) { // right choice
        prevChoice[t + 1] = 1;
        RPE = outcome[n, t] - Q[2];
        Q[2] = aF[n]*Q[2] + alpha[n]*RPE;
        Q[1] = aF[n]*Q[1];
      } else if (choice_R[n, t] == 0) { // left choice
        prevChoice[t + 1] = -1;
        RPE = outcome[n, t] - Q[1];
        Q[1] = aF[n]*Q[1] + alpha[n]*RPE;
        Q[2] = aF[n]*Q[2];
      }
    }
    choice_R[n, 1:Tsesh[n]] ~ bernoulli_logit(betaTemp[n]*Qdiff + 
                                              bias[n] + 
                                              k[n]*prevChoice[1:Tsesh[n]]);
  }
}
generated quantities{
  real Q[N, T + 1, 2];
  real Qdiff[N, T];
  real RPE[N, T];
  real prevChoice[N, T + 1];
  real choiceR_predicted[N, T];
  real log_lik[N]; // for log likelihood comparison

  for (n in 1:N) {
    Q[n,1,1] = 0;
    Q[n,1,2] = 0;
    prevChoice[n, 1] = 0;
    log_lik[n] = 0; // initialize log likelihood calculation
    for (t in 1:Tsesh[n]) {
      Qdiff[n,t] = Q[n,t,2] - Q[n,t,1];

      choiceR_predicted[n, t] = inv_logit(betaTemp[n]*Qdiff[n,t] +
                                                            bias[n] + 
                                                            k[n]*prevChoice[n,t]);
      log_lik[n] = log_lik[n] + bernoulli_logit_lpmf(choice_R[n,t] | betaTemp[n]*Qdiff[n,t] +
                                          bias[n] + 
                                          k[n]*prevChoice[n,t]);
      if (choice_R[n,t] == 1) { // right choice
        prevChoice[n,t+1] = 1;
        RPE[n,t] = outcome[n,t] - Q[n,t,2];
        Q[n,t+1,2] = aF[n]*Q[n,t,2] + alpha[n]*RPE[n,t];
        Q[n,t+1,1] = aF[n]*Q[n,t,1];
      } else if (choice_R[n,t] == 0) { // left choice
        prevChoice[n,t+1] = -1;
        RPE[n,t] = outcome[n,t] - Q[n,t,1];
        Q[n,t+1,1] = aF[n]*Q[n,t,1] + alpha[n]*RPE[n,t];
        Q[n,t+1,2] = aF[n]*Q[n,t,2];
      }
    }
  }
}