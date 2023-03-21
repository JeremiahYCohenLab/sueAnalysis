data{
  int<lower=1> N;
  int<lower=1> T;
  int<lower=1, upper=T> Tsesh[N];
  matrix[N, T] rpeMat;
  matrix[N, T] signalMat;
  matrix[N, T] nRwdMat;
}
parameters{
  real<lower=0> sigma[N];
  real<lower=0> maxF[N];
  real half_a;
  real slope_a;
  // real w_a;

  
  row_vector[N] half_pr;
  row_vector[N] slope_pr;
  row_vector[N] intercept_pr;
  row_vector[N] w_pr;

  real<lower=0> sigma_half;
  real<lower=0> sigma_slope;
  // real<lower=0> sigma_w;
}

transformed parameters{
  row_vector<lower=-2, upper=2>[N] half;
  row_vector[N] slope;
  row_vector[N] intercept;
  row_vector[N] w;

  half = 2 * Phi_approx(half_a + half_pr * sigma_half)-1;
  slope = 4 * Phi_approx(slope_a + slope_pr * sigma_slope);
  intercept = 2 * Phi_approx(intercept_pr);
  w = Phi_approx(w_pr);
}

model{
  sigma ~ cauchy(0, 1);
  maxF ~ normal(0, 1);
  intercept_pr ~ normal(0, 1);

  half_a ~ normal(0, 1);
  slope_a ~ normal(0, 1);
  // w_a ~ normal(0, 1);

  half_pr ~ normal(0, 1);
  slope_pr ~ normal(0 ,1);
  w_pr ~ normal(0, 1);

  sigma_half ~ cauchy(0, 1);
  sigma_slope ~ cauchy(0, 1);
  // sigma_w ~ cauchy(0, 1);


  for (n in 1:N) {
    row_vector[Tsesh[n]] est;
    row_vector[Tsesh[n]] currRpe;
    row_vector[Tsesh[n]] currNrwd; 
    row_vector[Tsesh[n]] currSignal;

    currRpe = rpeMat[n, 1:Tsesh[n]];
    currNrwd = nRwdMat[n, 1:Tsesh[n]];
    currSignal = signalMat[n, 1:Tsesh[n]];

    est = 1+exp(-slope[n]*(w[n]*currRpe + (1-w[n])*currNrwd - half[n]));
    est = maxF[n]./est + intercept[n]; 

    currSignal ~ normal(est, sigma[n]);
  }

}