rm ( list = ls() )
library("mvtnorm")
library("matlib")
library("devtools")
#install_github("cran/MCMCpack")#, for inverse wishart, but just for the moment
library("MCMCpack")
library("invgamma")
#library("gdirmn")
library(gjam)

#### INPUT OF THE DATA ####

n = 5 # number of sites
S = 10 # number of species
k = 4 # number of covariates (no intercept)

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
V = as.matrix(f$ydata) # matrix of n_species presence/absence data (in a continuous framework) (n * S)

#### CONSTRUCTION OF THE COVARIATES ####

#### HYPERPARAMETERS ####

#### SOME MATRICES TO BE PRECOMPUTED ####

mu_e = rep ( 0, times=S ) # Prior mean of e
R <- riwish(S + 1, diag(S)) # Prior variance-covariance matrix of the e 

# for inverse wishart and related computations see gjam library (for instance gjamSimData function), 
# it is done using Rcpp 

# r = number of latent factors, user chosen
#eta_h = rep ( 0, times=r )

r = 10
eta_h <- rinvgamma(r, 1/2, 1/10^4)

Dz <- riwish(2 + r - 1, 4 * 1/eta_h * diag(r)) 

# N_stick = truncation level for stick breaking factors, user chosen
N_stick = 30 
#Z <- rmvnorm ( n = N_stick, mean = rep(0, times = r), sigma = Dz ) 

#mu_Zj <- rep(0, times = r)
#for (j in seq(1, N)) {
#  Z[j,] <- rmvnorm ( n = 1, mean = mu_Zj, sigma = Dz ) 
#}

W <- rmvnorm ( n = n, mean = rep(0, times = r), sigma = diag(r) ) 
# need to transpose for the computation of Bx_i + Q(k)Zw_i'


nu = 1
G = 1
sigmaeps2 <- 1e4

muBeta = rep ( 0, times=k ) # Prior mean of the beta coefficients
sigmaBeta = diag ( k ) # Prior variance-covariance matrix of the beta coefficients
B <- matrix(data = 0, nrow = S, ncol = k) # coefficient matrix

pl = matrix(0,ncol = N_stick,nrow = n_species)

for (j in seq(1,S,1)) {
  B[j,] <- rmvnorm ( n = 1, mean = muBeta, sigma = 100*sigmaBeta )
}

# Dirichlet process

n_species = ncol ( V )
n_sites = n

# each label belongs to {1,...,N}
k = c();

p <- numeric(N_stick)
pl <- matrix(0, nrow = n_species, ncol = N_stick)

# Number of iterations
posteriorDraws = 3
burnInIterations = 0

alpha0 = 1e2 # Mass parameter of the Dirichlet process
compute_GD_prior <- function(N_stick, alpha0){
  xi <- numeric(N_stick)
  for (j in seq(1, N_stick-1)) {
    shape1 = alpha0/N_stick + sum(k == j)
    shape2 = (N_stick - j)/N_stick * alpha0
    for (s in seq(j + 1, N_stick)) {
      shape2 = shape2 + sum(k == s)
    }
    xi[j] <- rbeta(1, shape1 = shape1, shape2 = shape2)
  }
  return(xi)
}

norm_vec <- function(x) sqrt(sum(x^2))

for ( niter in 1:(posteriorDraws + burnInIterations) ) { # MCMC loop
  # Step 1 : resample the cluster assignments
  #for each species j = 1,...,n_species do
  Z <- matrix(data = 0, nrow = N_stick, ncol = r) 
  Q <- matrix(data = 0, nrow = n_species, ncol = N_stick) 
  
  xi <- compute_GD_prior(N_stick,alpha0)
  # c <- 2
  # given the vector k = [1,..., 1] of length = S
  p[1] <- xi[1]
  p[2:N_stick] <- sapply(2:N_stick, function(j) xi[j] * prod(1 - xi[1:(j-1)]))
  p[N_stick] <- 1 - sum( p[1:N_stick-1] )
  
  for ( j in 1:N_stick ) { 
    
    # Step 1
    if (!(j %in% k)) {
      Z[j,] = rmvnorm ( n = 1, mean = rep(0, times = r), sigma = Dz )
    } 
    else {
      cardinality_Sj = sum( j == k )
      Sigma_Zj = solve(cardinality_Sj/sigmaeps2 * t(W) %*% W + solve(Dz) )
      mu_Zj = c(0)
      for (l in 1:n_species) {
        if (k[l] == j){
          mu_Zj = mu_Zj + 1/sigmaeps2 * Sigma_Zj %*% t(W) %*% ( t(t(V[,l])) - x %*% B[l,] )
        }
      }
      Z[j,] = rmvnorm ( n = 1, mean = mu_Zj, sigma = Sigma_Zj )
      Q[j,k[l]] = 1
    }
  }
  
   
  A = Q %*% Z
  
  # Step 2
  for ( i in 1:n_sites ) {
    Sigma_W = solve(1/sigmaeps2 * t(A) %*% A + diag(r))
    mu_W = 1/sigmaeps2 * Sigma_W %*% A %*% (t(t(V[i,])) - B %*% x[i,])
    W[i,] = rmvnorm ( n = 1, mean = mu_W, sigma = Sigma_W )
  }
  
  # Step 3
  for ( l in 1:n_species ) { # Loop on species
    for (j in 1:N_stick) {
      pl[l,j] <- p[j] * exp( -1/(2*sigmaeps2) * norm_vec(V[,l] - x %*% B[l,] - W %*% Z[j,])^2)
    }
    k[l] = sample(1:N_stick, size = 1, replace=TRUE, prob = pl[l,])
  }
  
  # Step 4 : maybe once and for all outside all
  
  # Step 5
  dx = c(0)
  for (i in 1:n_sites) {
    dx = dx + norm_vec(V[i,] - B %*% x[i,] - A %*% W[i,])^2
  }
  #sigmaeps2 = invgamma::rinvgamma((n * S + nu)/2 + 1, dx/2 + nu/G^2)
  sigmaeps2 = 1
  
  # Step 6
  Dz = riwish(2 + r + N_stick - 1, t(Z) %*% Z + 4 * 1/eta_h * diag(r))
} 