rm ( list = ls() )
library("mvtnorm")
library(gjam)

#### INPUT OF THE DATA ####

#Simulated data are used to check that the algorithm can recover true
#parameter values and predict data, including underlying latent variables.  
#To illustrate I simulate a sample of size $n = 500$ for $S = 10$ species and
#$Q = 4$ predictors.  To indicate that all species are *continuous abundance*
#data I specify `typeNames` as `'CA'`.  `CA` data are continuous above zero, 
#with point mass at zero. 

f <- gjamSimData(n = 500, S = 10, Q = 4, typeNames = 'CA')
# The object `f` includes elements needed to analyze the simulated data set.  
# `f$typeNames` is a length-$S$ `character vector`. The `formula` follows 
# standard R syntax. It does not start with `y ~`, because gjam is multivariate.
# The multivariate response is supplied as a $n \times S$ `matrix` or
# `data.frame ydata`.
summary(f)
ml  <- list(ng = 1000, burnin = 100, typeNames = f$typeNames)
out <- gjam(f$formula, f$xdata, f$ydata, modelList = ml)
summary(out)

#### CONSTRUCTION OF THE COVARIATES ####

#### HYPERPARAMETERS ####

alpha0 = 1 # Mass parameter of the Dirichlet process

#### SOME MATRICES TO BE PRECOMPUTED ####

#### FIT OF THE MODEL ####

# Number of iterations
posteriorDraws = 1
burnInIterations = 0

for ( niter in 1:(posteriorDraws + burnInIterations) ) { # MCMC loop
  
  
  
  
} # End MCMC loop