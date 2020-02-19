data {
  int<lower=1, upper=1000> T;
  int<lower=0, upper=2> choice[T];
  int<lower=0, upper=1> outcome[T];
}
parameters {
  // Session-level raw parameters
  vector<lower=0, upper=1>[1] aN;    // learning rate for NPE
  vector<lower=0, upper=1>[1] aP;    // learning rate for PPE
  vector<lower=0, upper=1>[1] aF;    // forgetting rate
  vector<lower=0, upper=20>[1] beta;  // inverse temperature
  vector<lower=0, upper=1>[1] v;     // expected average value learning rate
}

model {
  vector[2] Q;  // initial values for Q
  real PE;
  real rBar;
  real aF;

  // individual parameters
  aN   ~ normal(0.5, 1);
  aP   ~ normal(0.5, 1);
  aF   ~ normal(0.7, 0.01);
  beta ~ normal(5, 5);
  v    ~ normal(0.005, 0.01);



  Q = rep_vector(0.0, 2);
  rBar = 0.4;


  for (t in 1:(T)) {
    // compute action probabilities
    choice[t] ~ categorical_logit(beta[1] * Q);

    // prediction error
    PE = outcome[t] - Q[choice[t]] - rBar;

    // value updating (learning)
    if (PE < 0){
      Q[choice[t]] = Q[choice[t]] + aN[1] * PE;
    }
    else{
      Q[choice[t]] = Q[choice[t]] + aP[1] * PE;
    }
    if (choice[t] == 1){
      Q[2] = Q[2] * aF[1];
    }else{
      Q[1] = Q[1] * aF[1];
    }
    rBar = v[1] * outcome[t] + (1-v[1]) * rBar;
  }
}