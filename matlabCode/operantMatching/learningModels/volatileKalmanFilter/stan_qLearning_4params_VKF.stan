data {
  int<lower=1> N;
  int<lower=1> T;
  int<lower=1, upper=T> Tsesh[N];
  int<lower=0, upper=2> choice[N, T];
  int<lower=0, upper=1> outcome[N, T];
}
transformed data {
  vector[2] initM;  // initial values for M
  initM = rep_vector(0.0, 2);
  vector[2] initW;  // initial values for W
  initW = rep_vector(0.0, 2);  
}
parameters {
// Declare all parameters as vectors for vectorizing
  // Hyper(animal)-parameters
  vector[4] mu_p;
  vector<lower=0>[4] sigma;

  // Session-level raw parameters
  vector[N] lambda_pr;      // volatility update rate
  vector[N] nu0_pr;        	// initial volatility
  vector[N] omega_pr;       // observation noise
  vector[N] beta_pr;        // inverse temp for action choice
}
transformed parameters {
  // session-level parameters
  vector<lower=0, upper=2>[N] lambda;
  vector<lower=0, upper=2>[N] nu0;
  vector<lower=0, upper=2>[N] omega;
  vector<lower=0, upper=10>[N] beta;

  for (n in 1:N) {
    lambda[n] 	= Phi_approx(mu_p[1] + sigma[1] * lambda_pr[n]) * 2;
    nu0[n]   	  = Phi_approx(mu_p[2] + sigma[2] * nu0_pr[n]) * 2;
    omega[n]    = Phi_approx(mu_p[3] + sigma[3] * omega_pr[n]) * 2;
    beta[n]    	= Phi_approx(mu_p[4] + sigma[4] * beta_pr[n]) * 10;
  }
}
model {
  // Hyperparameters
  mu_p  ~ normal(0, 0.5);
  sigma ~ cauchy(0, 0.5);

  // individual parameters
  lambda_pr 	~ normal(0, 1);
  nu0_pr   		~ normal(0, 1);
  omega_pr    ~ normal(0, 1);
  beta_pr    	~ normal(0, 1);

  // session loop and trial loop
  for (n in 1:N) {
    vector[2] M;  		// estimated mean
    vector[2] M_pre;  	// estimated mean t-1
    vector[2] W;		// estimated variance
    vector[2] W_pre;	// estimated variance t-1
    vector[2] W_auto;	// variance autocorrelation
    real k;
    vector[2] nu;
    vector[2] nu_pre;
    vector[Tsesh[n]] Mdiff;

    int c;

    M_pre = initM;
    W_pre = initW;
    nu_pre[1] = nu0[n];
    nu_pre[2] = nu0[n];

    for (t in 1:(Tsesh[n])) {
      Mdiff[t] = M[2] - M[1];
      choice[n, t] ~ bernoulli_logit(beta[n] * Mdiff[t]);
      c = choice[n, t] + 1;

	  k = (W_pre[c] + nu_pre[c]) / (W_pre[c] + nu_pre[c] + omega[n]);
	  M[c] = M_pre[c] + k * (outcome[n, t] - M_pre[c]);
	  W[c] = (1 - k) * (W_pre[c] + nu_pre[c]);
	  W_auto[c] = (1 - k) * W_pre[c];
	  nu[c] = nu_pre[c] + lambda[n] * ((M[c] - M_pre[c]) * (M[c] - M_pre[c]) + W_pre[c] + W[c] - 2 * W_auto[c] - nu_pre[c]);
	  // below, a slightly shorter way to calculate nu
	  // nu[c] = nu_pre[c] + lambda[n] * ((M[c] - M_pre[c]) * (M[c] - M_pre[c]) + k * (W_pre[c] - nu_pre[c]));

	  M_pre[c] = M[c];
	  W_pre[c] = W[c];
	  nu_pre[c] = nu[c];
    }
  }
}
generated quantities {
  // For group level parameters
  real<lower=0, upper=2> mu_lambda;
  real<lower=0, upper=2> mu_nu0;
  real<lower=0, upper=2> mu_omega;
  real<lower=0, upper=10> mu_beta;

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

  mu_lambda = Phi_approx(mu_p[1]) * 2;
  mu_nu0	  = Phi_approx(mu_p[2]) * 2;
  mu_omega  = Phi_approx(mu_p[3]) * 2;
  mu_beta   = Phi_approx(mu_p[4]) * 10;

  { // local section, this saves time and space
    for (n in 1:N) {
	  vector[2] M;  		// estimated mean
   	vector[2] M_pre;  	// estimated mean t-1
    vector[2] W;		// estimated variance
	  vector[2] W_pre;	// estimated variance t-1
	  vector[2] W_auto;	// variance autocorrelation
	  real k;
	  vector[2] nu;
	  vector[2] nu_pre;
	  vector[Tsesh[n]] Mdiff;

	  int c;

	  M_pre = initM;
	  W_pre = initW;
 	  nu_pre[1] = nu0[n];
	  nu_pre[2] = nu0[n];

	  log_lik[n] = 0;

      for (t in 1:(Tsesh[n])) {
        Mdiff[t] = M[2] - M[1];
        // compute log likelihood of current trial
        log_lik[n] += bernoulli_logit_lpmf(choice[n, t] | beta[n] * Mdiff[t]);

        // generate posterior prediction for current trial
  //      y_pred[n, t] = categorical_rng(softmax(beta[n] * Q));
 	    c = choice[n, t] + 1;

	  	k = (W_pre[c] + nu_pre[c]) / (W_pre[c] + nu_pre[c] + omega[n]);
	  	M[c] = M_pre[c] + k * (outcome[n, t] - M_pre[c]);
	  	W[c] = (1 - k) * (W_pre[c] + nu_pre[c]);
	  	W_auto[c] = (1 - k) * W_pre[c];
	  	nu[c] = nu_pre[c] + lambda[n] * ((M[c] - M_pre[c]) * (M[c] - M_pre[c]) + W_pre[c] + W[c] - 2 * W_auto[c] - nu_pre[c]);
	  	// below, a slightly shorter way to calculate nu
	  	// nu[c] = nu_pre[c] + lambda[n] * ((M[c] - M_pre[c]) * (M[c] - M_pre[c]) + k * (W_pre[c] - nu_pre[c]));

	  	M_pre[c] = M[c];
	  	W_pre[c] = W[c];
	  	nu_pre[c] = nu[c];
	  }
    }
  }
}
