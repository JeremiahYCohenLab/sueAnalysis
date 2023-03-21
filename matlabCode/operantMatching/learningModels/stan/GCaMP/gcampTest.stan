data{
  int<lower=1> N;
  int<lower=1> T;
  int<lower=1, upper=T> Tsesh[N];
  matrix[N, T] rpeMat;
  matrix[N, T] signalMat;
}
parameters{
  real<lower=0> sigma[N];
  real<lower=0> maxF[N];
  real half_a;
  real slope_a;

  
  row_vector[N] half_pr;
  row_vector[N] slope_pr;
  row_vector[N] intercept_pr;

  real<lower=0> sigma_half;
  real<lower=0> sigma_slope;
}

transformed parameters{
  row_vector<lower=-2, upper=2>[N] half;
  row_vector[N] slope;
  row_vector[N] intercept;

  half = 2 * Phi_approx(half_a + half_pr * sigma_half)-1;
  slope = 4 * Phi_approx(slope_a + slope_pr * sigma_slope);
  intercept = 2 * Phi_approx(intercept_pr);
}

model{
  sigma ~ cauchy(0, 1);
  maxF ~ normal(0, 1);
  intercept_pr ~ normal(0, 1);

  half_a ~ normal(0, 1);
  slope_a ~ normal(0, 1);

  half_pr ~ normal(0, 1);
  slope_pr ~ normal(0 ,1);

  sigma_half ~ cauchy(0, 1);
  sigma_slope ~ cauchy(0, 1);


  for (n in 1:N) {
    row_vector[Tsesh[n]] est;
    row_vector[Tsesh[n]] currRpe;
    row_vector[Tsesh[n]] currSignal;

    currRpe = rpeMat[n, 1:Tsesh[n]];
    currSignal = signalMat[n, 1:Tsesh[n]];

    est = 1+exp(-slope[n]*(currRpe - half[n]));
    est = maxF[n]./est + intercept[n]; 

    currSignal ~ normal(est, sigma[n]);
  }

}