data {
  int<lower=1> T;                   // number of observations (length)
  int<lower=1> K;                   // number of hidden states
  int<lower=1> N;                   // number of possible observations
  int y[T];                        // observations

}

parameters {
  // Discrete state model
  simplex[K] pi1;                   // initial state probabilities
  simplex[K] expP;                  // exploration transition probabilities
                                    // A[i][j] = p(z_t = j | z_{t-1} = i)
  simplex[N] expO;                  // exploration emission
  vector<lower = 0, upper = 1>[K-1] pStay;                  // exploitation transition probabilities

}

transformed parameters{
  matrix[K,K] A;
  matrix[K,N] B;

  for (i in 1:K){
    A[1,i] = expP[i]; 
  }

  for (i in 1:N){
    B[1,i] = expO[i];
  }

  for (i in 2:K){
    for (j in 2:K){
      if (i != j){
        A[i,j] = 0.000000001;
      } else{
        A[i,j] = pStay[i-1];
        A[i,1] = 1-pStay[i-1];
      }
    }
    for (j in 1:N){
      if (j != i-1){
        B[i,j] = 0.000000001;
      }else{
          B[i,j] = 1-0.000000001;
      }
    }
  }
}

model {
  pStay ~ beta(8, 2);
  // print(pStay[1]);
  matrix[T,K] logalpha;
// Forward algorithm log p(z_t = j | x_{1:t})
    real accumulator[K];

    for (j in 1:K){
      logalpha[1,j] = log(pi1[j]) + log(B[j,y[1]]);
    }

    for (t in 2:T) {
      for (j in 1:K) { // j = current (t)
        for (i in 1:K) { // i = previous (t-1)
                         // Murphy (2012) Eq. 17.48
                         // belief state      + transition prob + local evidence at t
          accumulator[i] = logalpha[t-1, i] + log(A[i, j]) + log(B[j,y[t]]);
        }
        logalpha[t, j] = log_sum_exp(accumulator);
      }
    }// Forward

  target += log_sum_exp(logalpha[T]); // Note: update based only on last logalpha


}

// generated quantities {
//   vector[K] logbeta[T];
//   vector[K] loggamma[T];

//   vector[K] alpha[T];
//   vector[K] beta[T];
//   vector[K] gamma[T];

//   int<lower=1, upper=K> zstar[T];
//   real logp_zstar;

//   { // Forward algortihm
//     for (t in 1:T)
//       alpha[t] = softmax(logalpha[t]);
//   } // Forward

//   { // Backward algorithm log p(x_{t+1:T} | z_t = j)
//     real accumulator[K];

//     for (j in 1:K)
//       logbeta[T, j] = 1;

//     for (tforward in 0:(T-2)) {
//       int t;
//       t = T - tforward;

//       for (j in 1:K) { // j = previous (t-1)
//         for (i in 1:K) { // i = next (t)
//                          // Murphy (2012) Eq. 17.58
//                          // backwards t    + transition prob + local evidence at t
//           accumulator[i] = logbeta[t, i] + log(A[j, i]) + normal_lpdf(y[t] | mu[i], sigma[i]);
//           }
//         logbeta[t-1, j] = log_sum_exp(accumulator);
//       }
//     }

//     for (t in 1:T)
//       beta[t] = softmax(logbeta[t]);
//   } // Backward

//   { // Forwards-backwards algorithm log p(z_t = j | x_{1:T})
//     for(t in 1:T) {
//         loggamma[t] = alpha[t] .* beta[t];
//     }

//     for(t in 1:T)
//       gamma[t] = normalize(loggamma[t]);
//   } // Forwards-backwards

//   { // Viterbi algorithm
//     int bpointer[T, K];             // backpointer to the most likely previous state on the most probable path
//     real delta[T, K];               // max prob for the seq up to t
//                                     // with final output from state k for time t

//     for (j in 1:K)
//       delta[1, K] = normal_lpdf(y[1] | mu[j], sigma[j]);

//     for (t in 2:T) {
//       for (j in 1:K) { // j = current (t)
//         delta[t, j] = negative_infinity();
//         for (i in 1:K) { // i = previous (t-1)
//           real logp;
//           logp = delta[t-1, i] + log(A[i, j]) + normal_lpdf(y[t] | mu[j], sigma[j]);
//           if (logp > delta[t, j]) {
//             bpointer[t, j] = i;
//             delta[t, j] = logp;
//           }
//         }
//       }
//     }

//     logp_zstar = max(delta[T]);

//     for (j in 1:K)
//       if (delta[T, j] == logp_zstar)
//         zstar[T] = j;

//     for (t in 1:(T - 1)) {
//       zstar[T - t] = bpointer[T - t + 1, zstar[T - t + 1]];
//     }
//   }
// }