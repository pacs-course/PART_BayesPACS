///////////////////////// DATA /////////////////////////////////////
  data {
    int<lower = 0> n;       		  // number of data
    int<lower = 0> d;       		  // number of coefficients of regression
    int<lower = 0, upper = 1> Y[n];   // response vector
    matrix[d,n] X;   		  		  // design matrix
    
  }

//////////////////// PARAMETERS /////////////////////////////////
  parameters {
    row_vector[d] beta;        // coefficients of regression
  }


////////////////// MODEL ////////////////////////
  model {
 
  // Likelihood     
    for (s in 1:n)
    {
      Y[s] ~ bernoulli(inv_logit(-3 + beta*X[,s]));  
    } 
    
    // Prior
    // beta
    for (j in 1:(d)) 
    {
      beta[j] ~ normal(0, 5);
    }

  }