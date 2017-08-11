##### Modify the following path according to your system!!
setwd('PART_BayesPACS-master/Interface/Rdir')

library(MASS)
library(rstan)

subX <- as.matrix(read.table("subdata.csv", header = F))
d    <- dim(subX)[2]
beta <- rep(0.1,d)

stan_rdump(c('beta'), "inits.init.R")






