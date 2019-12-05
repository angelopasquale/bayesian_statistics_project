rm ( list = ls() )
library("mvtnorm")
library("matlib")
library("cIRT")
library(gjam)

#### INPUT OF THE DATA ####

#Simulated data are used to check that the algorithm can recover true
#parameter values and predict data, including underlying latent variables.  
#To illustrate I simulate a sample of size $n = 500$ for $S = 10$ species and
#$Q = 4$ predictors.  To indicate that all species are *continuous abundance*
#data I specify `typeNames` as `'CA'`.  `CA` data are continuous above zero, 
#with point mass at zero. 

n = 500
S = 10
k = 4

f <- gjamSimData(n = n, S = S, Q = k, typeNames = 'CA')
# The object `f` includes elements needed to analyze the simulated data set.  
# `f$typeNames` is a length-$S$ `character vector`. The `formula` follows 
# standard R syntax. It does not start with `y ~`, because gjam is multivariate.
# The multivariate response is supplied as a $n \times S$ `matrix` or
# `data.frame ydata`.
summary(f)
#ml  <- list(ng = 1000, burnin = 100, typeNames = f$typeNames)
#out <- gjam(f$formula, f$xdata, f$ydata, modelList = ml)
#summary(out)

x = as.matrix(f$xdata) # matrix of measured (environmental) covariates (n * k)
v = as.matrix(f$ydata) # matrix of species presence/absence data (in a continuous framework) (n * S)

#### CONSTRUCTION OF THE COVARIATES ####

#### HYPERPARAMETERS ####

alpha0 = 1 # Mass parameter of the Dirichlet process

#### SOME MATRICES TO BE PRECOMPUTED ####
B <- matrix(rnorm(S * k, 0, 100), S) # true covariates
R <- riwishart(S + 1, diag(n))
e <- matrix(rnorm(n * S, 0, 100), n)

zsim <- x %*% t(B) + e

#### FIT OF THE MODEL ####

# Number of iterations
posteriorDraws = 1
burnInIterations = 0

for ( niter in 1:(posteriorDraws + burnInIterations) ) { # MCMC loop
  
  
  
  
} # End MCMC loop