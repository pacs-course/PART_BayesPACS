
///////////////////////// DATA /////////////////////////////////////
  data {
    int<lower = 0> N;                 // number of users
    int<lower = 0> Q;                 // number of user covariates
    int<lower = 0, upper = 1> Y[N];
    matrix[Q,N] X;   		              // design matrix of user covariates
    vector[N] f;   		              // design matrix of books
    matrix[Q,Q] L0;               // empirical variance-covariance matrix of beta
    int<lower = 0> L1;          // empirical variance of alpha 
    vector[Q] m0;                 // mean of beta
    int<lower = 0> m1;                 // mean of alpha
  }

//////////////////// PARAMETERS /////////////////////////////////
  parameters {
    row_vector[Q] beta;        // coefficients related to users
    real alpha;       // coefficients related to books
    real mu;          // user-level random effect 
    real<lower = 0> tau;
    real<lower = 0> eta;
    
  }


////////////////// MODEL ////////////////////////
  model {
    
    // Likelihood
    for (i in 1:N){
      target += bernoulli_logit_lpmf(Y | beta * X[,i] + alpha * f[i] + mu);
    }
      

    
    // Priors
    
    mu ~ normal( 0, sqrt(100));
    beta ~ multi_normal( m0, inv_sqrt(tau)*L0);
    alpha ~ normal( m1, inv_sqrt(eta)*L1);
    tau ~ gamma(0.01, 0.01);
    eta ~ gamma(0.01, 0.01);
  
  }


