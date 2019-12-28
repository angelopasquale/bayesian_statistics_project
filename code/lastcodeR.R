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
S = 15 # number of species
n_cov = 4 # number of covariates (no intercept)

f <- gjamSimData(n = n, S = S, Q = n_cov, typeNames = 'CA')
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

r = 5
eta_h <- rinvgamma(r, 1/2, 1/10^4)

Dz <- riwish(2 + r - 1, 4 * 1/eta_h * diag(r)) 

# N_stick = truncation level for stick breaking factors, user chosen
N_stick = 10 
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

muBeta = rep ( 0, times=n_cov ) # Prior mean of the beta coefficients
sigmaBeta = diag ( k ) # Prior variance-covariance matrix of the beta coefficients
B <- matrix(data = 0, nrow = S, ncol = n_cov) # coefficient matrix

pl = matrix(0,ncol = N_stick,nrow = n_species)

for (j in seq(1,S,1)) {
  B[j,] <- rmvnorm ( n = 1, mean = muBeta, sigma = 100*sigmaBeta )
}

cardinality_S = rep ( 0, times=N_stick )

# Dirichlet process

n_species = ncol ( V )
n_sites = n

# each label belongs to {1,...,N_stick}
k = c();

p <- numeric(N_stick)
pl <- matrix(0, nrow = n_species, ncol = N_stick)

# Number of iterations
posteriorDraws = 3
burnInIterations = 0

alpha0 = 1e1 # Mass parameter of the Dirichlet process

###############################################
#STA FUNZIONE è FATTA NELL'hfunction DA .sampleP <- CHIAMATA IN .getPars <-, HO PROVATO A DARCI UN OCCHIO (LA CHIAMA IN 
# tmp <- .getPars, MA CI VUOLE UN BEL PO DI TEMPO. SE CI CAPITE QUALCOSA ABBIAMO RISOLTO IL PROBLEMA DEI p,
#MA COME ABBIAMO FATTO NOI MI SEMBRA TORNI A LIVELLO DI NOTAZIONI CON QUELLO CHE HA FATTO LUI IN hfunction
###############################################

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

xi <- compute_GD_prior(N_stick,alpha0)
# c <- 2
# given the vector k = [1,..., 1] of length = S
p[1] <- xi[1]
p[2:N_stick-1] <- sapply(2:N_stick, function(j) xi[j] * prod(1 - xi[1:(j-1)]))
p[N_stick] <- 1 - sum( p[1:N_stick-1] )

k = sample(1:N_stick, size = n_species, replace=TRUE, prob = p)

names(k) <- seq(1,n_species)

for ( niter in 1:(posteriorDraws + burnInIterations) ) { # MCMC loop
  
  # Step 1 : resample the cluster assignments
  #for each species j = 1,...,n_species do
  Z <- matrix(data = 0, nrow = N_stick, ncol = r) 
  Q <- matrix(data = 0, nrow = n_species, ncol = N_stick) 
  
  # Step 3
  for ( l in 1:n_species ) { # Loop on species
    for (j in 1:N_stick) {
      pl[l,j] <- p[j] * exp( -1/(2*sigmaeps2) * norm_vec(V[,l] - x %*% B[l,] - W %*% Z[j,])^2)
    }
    k[l] = sample(N_stick, size = 1, replace=TRUE, prob = pl[l,])
  }
  
  cardinality_S = rep(0, times = N_stick)
                      
  for ( j in 1:N_stick ) { 
    # Step 1
      print(j)
    if (!(j %in% k)) {
      print('entered if')
      Z[j,] = rmvnorm ( n = 1, mean = rep(0, times = r), sigma = Dz )
      cardinality_S[j] = 1
    } 
    else { # if j in k
      print('entered else')
      cardinality_S[j] = cardinality_S[j] + 1
      print(cardinality_S[j])
      Sigma_Zj = solve(cardinality_S[j]/sigmaeps2 * t(W) %*% W + solve(Dz) )
      mu_Zj = c(0)
      for (l in 1:n_species) {
        if (k[l] == j){
          print('entered else-if')
          mu_Zj = mu_Zj + 1/sigmaeps2 * Sigma_Zj %*% t(W) %*% ( t(t(V[,l])) - x %*% B[l,] )
        }
        Q[l,k[l]] = 1
      }
      Z[j,] = rmvnorm ( n = 1, mean = mu_Zj, sigma = Sigma_Zj )
    }
  }
  
  A = Q %*% Z

} 

print(k)
print(A)
