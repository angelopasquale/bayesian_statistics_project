rm ( list = ls() )
library("mvtnorm")
library("matlib")
library("devtools")
#install_github("cran/MCMCpack")#, for inverse wishart, but just for the moment
library("MCMCpack")
library(gjam)

#### INPUT OF THE DATA ####

#Simulated data are used to check that the algorithm can recover true
#parameter values and predict data, including underlying latent variables.  
#To illustrate I simulate a sample of size $n = 500$ for $S = 10$ species and
#$Q = 4$ predictors.  To indicate that all species are *continuous abundance*
#data I specify `typeNames` as `'CA'`.  `CA` data are continuous above zero, 
#with point mass at zero. 

n = 5
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
muBeta = rep ( 0, times=k ) # Prior mean of the beta coefficients
sigmaBeta = diag ( k ) # Prior variance-covariance matrix of the beta coefficients
B <- matrix(data = 0, nrow = S, ncol = k) # true covariates

for (j in seq(1,S,1)) {
  B[j,] <- rmvnorm ( n = 1, mean = muBeta, sigma = 100*sigmaBeta )
}

mue = rep ( 0, times=S ) # Prior mean of e
R <- riwish(S + 1, diag(S)) # Prior variance-covariance matrix of the e 

e <- matrix(data = 0, nrow = n, ncol = S) # error

for (i in seq(1,n,1)) {
  e[i,] <- rmvnorm ( n = 1, mean = mue, sigma = R )
}

vsim <- x %*% t(B) + e
# How to interpret these values? They do not seem to have any sense, moreover there are any hyperparameters
# in the Core Model, how to write a first basic Gibbs Sampler fot it?

# for inverse wishart and related computations see gjam library (for instance gjamSimData function), 
# it is done using Rcpp 
