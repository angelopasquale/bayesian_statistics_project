rm ( list = ls() )
library("mvtnorm")
library("gdirmn")
library("matlib")
library("devtools")
#install_github("cran/MCMCpack")#, for inverse wishart, but just for the moment
library("MCMCpack")
library("invgamma")
#library("gdirmn")
library(gjam)

#### INPUT OF THE DATA ####

#Simulated data are used to check that the algorithm can recover true
#parameter values and predict data, including underlying latent variables.  
#To illustrate I simulate a sample of size $n = 500$ for $S = 10$ n_species and
#$Q = 4$ predictors.  To indicate that all n_species are *continuous abundance*
#data I specify `typeNames` as `'CA'`.  `CA` data are continuous above zero, 
#with point mass at zero. 

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
v = as.matrix(f$ydata) # matrix of n_species presence/absence data (in a continuous framework) (n * S)

#### CONSTRUCTION OF THE COVARIATES ####

#### HYPERPARAMETERS ####

#### SOME MATRICES TO BE PRECOMPUTED ####

mu_e = rep ( 0, times=S ) # Prior mean of e
R <- riwish(S + 1, diag(S)) # Prior variance-covariance matrix of the e 

e <- matrix(data = 0, nrow = n, ncol = S) # error

for (i in seq(1,n,1)) {
  e[i,] <- rmvnorm ( n = 1, mean = mu_e, sigma = R )
}

vsim <- x %*% t(B) + e
# How to interpret these values? They do not seem to have any sense, moreover there are any hyperparameters
# in the Core Model, how to write a first basic Gibbs Sampler fot it?

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

mu_Zj <- c(0)
#for (j in seq(1, N)) {
#  Z[j,] <- rmvnorm ( n = 1, mean = mu_Zj, sigma = Dz ) 
#}

W <- rmvnorm ( n = n, mean = rep(0, times = r), sigma = diag(r) ) 
# need to transpose for the computation of Bx_i + Q(k)Zw_i'

muBeta = rep ( 0, times=k ) # Prior mean of the beta coefficients
sigmaBeta = diag ( k ) # Prior variance-covariance matrix of the beta coefficients
B <- matrix(data = 0, nrow = S, ncol = k) # coefficient matrix

for (j in seq(1,S,1)) {
  B[j,] <- rmvnorm ( n = 1, mean = muBeta, sigma = 100*sigmaBeta )
}

sigma_eps2 = 1e6

# Dirichlet process

n_species = ncol ( v )

# each label belongs to {1,...,N}
#labels = rep ( 1, times = n_species )
labels = seq(1,n_species)
clusterSize = rep(0,N_stick)
names(clusterSize) <- seq(1,N_stick)
cluster_dx = matrix(0,nrow = k+1, ncol = N_stick)
colnames(cluster_dx) <- seq(0,N_stick)
cluster_dx2 = matrix(0,nrow = k+1, ncol = N_stick)
colnames(cluster_dx2) <- seq(0,N_stick)
clusterMatrix = matrix ( rep(1, times=n_species*n_species), nrow = n_species, ncol = n_species )

p <- rep ( 1/n_species, times = N_stick)
pl <- rep ( 0, times = N_stick)

# Number of iterations
posteriorDraws = 1
burnInIterations = 0

alpha0 = 1e-1 # Mass parameter of the Dirichlet process

norm_vec <- function(x) sqrt(sum(x^2))

compute_GD_prior <- function(N_stick, alpha0){
  xi <- numeric(N_stick)
  for (j in seq(1, N_stick-1)) {
    shape1 = alpha0/N_stick + sum(labels == j)
    shape2 = (N_stick - j)/N_stick * alpha0
    for (s in seq(j + 1, N_stick)) {
      shape2 = shape2 + sum(labels == s)
    }
    xi[j] <- rbeta(1, shape1 = shape1, shape2 = shape2)
  }
  return(xi)
}
Z <- matrix(data = 0, nrow = N_stick, ncol = r) 

# UPDATE of p
xi <- compute_GD_prior(N_stick,alpha0)
# c <- 2
# given the vector labels = [1,..., 1] of length = S
p[1] <- xi[1]
p[2:N_stick] <- sapply(2:N_stick, function(j) xi[j] * prod(1 - xi[1:(j-1)]))
p[N_stick] <- abs(1 - sum( p[1:N_stick-1] ))

for ( niter in 1:(posteriorDraws + burnInIterations) ) { # MCMC loop
  # Step 1 : resample the cluster assignments
  #for each species l = 1,...,n_species do
  mu_Zj <- c(0)
  for ( l in 1:n_species ) { # Loop on species
    print(labels)
    
    cat ( "iter: ", niter, " species: ", l, "\n", sep = "" )
    
    for (j in seq(1,N_stick)) {
      #print(j)
      cluster_dx[,labels[l]] = cluster_dx[,labels[l]] + (v[,l] - x %*% B[l,])
      cluster_dx2[,labels[l]] = cluster_dx2[,labels[l]] + (v[,l] - x %*% B[l,] - W %*% Z[j, ])
      
      if ( !(j %in% labels) ){ 
        # if j is not a label (e.g. not a stick position), we generate it 
        Z[j, ] = rmvnorm ( n = 1, mean = rep(0, times = r), sigma = Dz )
        print('entered if')
        clusterSize[labels[l]] = 1
      }
      else {
        # if j is a label, we count its occurencies, i.e. the cluster size
        # we work on unique labels vector in order to avoid useless loops  
        uniqueLabels = unique(labels)
        for ( c in 1:length(uniqueLabels) ) {
          #print(which(labels == j))
          if ( uniqueLabels[c] == labels[l] ) {
            clusterSize[labels[l]] = clusterSize[labels[l]] + 1
            #print(uniqueLabels[c])
            #print(l)
            print('entered ELSE')
            
            Sigma_Zj = solve(clusterSize[labels[l]]/sigma_eps2 * t(W) %*% W + solve(Dz))
            mu_Zj = (Sigma_Zj %*% t(W) / sigma_eps2) %*% cluster_dx[,labels[l]]
            Z[j, ] = rmvnorm ( n = 1, mean = mu_Zj, sigma = Sigma_Zj ) 
            
          }
        }
        
        #clusterSize[labels[l]] = clusterSize[labels[l]] + sum(labels[l] == j)
      }
      
      # update labels according to full conditional 3 
      
      pl[1:N_stick] <- p[1:N_stick] * exp( - 1 / (2*sigma_eps2) * norm_vec(cluster_dx2)^2 )
      
      pl[1:N_stick] <- pl[1:N_stick] / sum(pl[1:N_stick])
      
      #cat ( "\014" ) # clear the screen
      #print(p)
      #print(pl)
      labels[l] = sample(1:N_stick, size = 1, replace=TRUE, prob = pl)
      
      # UPDATE of p
      xi <- compute_GD_prior(N_stick,alpha0)
      # c <- 2
      # given the vector labels = [1,..., 1] of length = S
      p[1] <- xi[1]
      p[2:N_stick] <- sapply(2:N_stick, function(j) xi[j] * prod(1 - xi[1:(j-1)]))
      p[N_stick] <- abs(1 - sum( p[1:N_stick-1] ))
      
      # xi <- compute_GD_prior(N_stick,alpha0)
      # # c <- 2
      # # given the vector labels = [1,..., 1] of length = S
      # p[1] <- xi[1]
      # p[2:N_stick] <- sapply(2:N_stick, function(j) xi[j] * prod(1 - xi[1:(j-1)]))
      # p[N_stick] <- 1 - sum( p[1:N_stick-1] )
      #print(labels)
    }
  }
  
  #   uniqueLabels = unique(labels)
  #   prob.old = rep ( 0, times = length(uniqueLabels) ) # Probability of assigning to an existing cluster
  #   twopisigma = (2 * pi * sigma2) ^ (-length(times) / 2)
  #   for ( c in 1:length(uniqueLabels) ) {
  #     if ( uniqueLabels[c] != labels[i] ) {
  #       dx = y[i,] - h %*% beta[c,]
  #       prob.old[c] = clusterSize[c] * twopisigma * exp ( -t(dx) %*% dx / (2*sigma2) )
  #     }
  #   }
  #   
  #   # The i-th cell is removed from its cluster
  #   clusterMatrix[,i] = rep ( 0, times = cells )
  #   clusterMatrix[i,] = rep ( 0, times = cells )
  #   clusterSize[labels[i]] = clusterSize[labels[i]] - 1
  #   
  #   if ( clusterSize[labels[i]] == 0 )
  #     beta[labels[i],] = rep ( NA, times = ncol(h) )
  #   
  #   # Then, pick the new cluster for the i-th cell
  #   a = runif ( 1, 0, sum(prob.old) + prob.new )
  #   if ( a < prob.new ) {
  #     labels[i] = max(labels) + 1
  #     clusterMatrix[i,i] = 1
  #     clusterSize[labels[i]] = 1
  #   }
  #   else {
  #     a = a - prob.new
  #     for ( j in (1:cells)[-i] ) {
  #       if ( a < prob.old[j] ) {
  #         labels[i] = labels[j];
  #         clusterMatrix[j,i] = 1;
  #         clusterMatrix[i,] = clusterMatrix[j,]
  #         clusterMatrix[,i] = clusterMatrix[i,]
  #         clusterSize[labels[i]] = clusterSize[labels[i]] + 1
  #         break
  #       }
  #       else a = a - prob.old[j]
  #     }
  #   }
  #   
  #   # If the cluster of i contains only i, we need to resample its parameters
  #   if ( clusterSize[labels[i]] == 1 ) {
  #     while ( nrow(beta) < labels[i] )
  #       beta = rbind ( beta, rep ( NA, times = ncol(h) ) )
  #     
  #     beta[labels[i],] = rmvnorm ( 1, mean = betaBar, sigma = sigma )
  #   }
}
#labels
#Z

# TODO:
# need to update cluster_size in the loop
