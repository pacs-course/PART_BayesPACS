##### Modify the following path according to your system!!
setwd('PART_BayesPACS-master/Interface/Rdir')

library(MASS)
library(rstan)

subX <- as.matrix(read.table("subdata.csv", header = F))

# number of data
n <- dim(subX)[1] 

# number of covariates
d <- dim(subX)[2]

# response vector
b0 <- -3
b  <-  rnorm(d, 0, 25)

odd  <- exp(cbind(b0*rep(1, n),b*subX))
prob <- odd/(odd +1)

Y <- NULL
for (j in seq(n)){
  Y[j] <- rbinom(1,1,prob[j])
}

# design matrix
X <- t(subX)

stan_rdump(c('n','d', 'Y','X'), "subdata.data.R")


