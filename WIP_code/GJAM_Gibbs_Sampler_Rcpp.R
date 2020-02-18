setwd("/Users/angelopasquale/Documents/University/LM/YEAR2/SEM1/BS/Project/implementation/gjam/")
Rcpp::sourceCpp('src/cppFns.cpp') #in gjam sources
source(file = "R/gjamHfunctions.R")

GJAM_Gibbs_Sampler_Rcpp <- function(x, Y, r, N_stick, alpha0, posteriorDraws, burnInIterations){
  tic("Data generation:")
  # INPUT PARAMETERS
  n = nrow(Y) # number of sites
  S = ncol(Y) # number of species
  n_cov = ncol(x) # number of covariates (no intercept)
  # r = number of latent factors
  # N_stick = truncation level for stick breaking factors, user chosen
  # alpha0 = Mass parameter of the Dirichlet process
  
  # OUTPUT PARAMETERS
  list_B <- list()
  list_A <- list()
  list_Z <- list()
  list_k <- list()
  list_R <- list()
  list_sigmaeps2 <- list()
  
  # 3.A.2 for simulation of latent variable V
  #V = Y
  
  #### CONSTRUCTION OF THE COVARIATES ####
  
  #### HYPERPARAMETERS ####
  
  #### SOME MATRICES TO BE PRECOMPUTED ####
  
  mu_e = rep ( 0, times=S ) # Prior mean of e
  R <- .riwish(S + 1, diag(S)) # Prior variance-covariance matrix of the e 
  
  # for inverse wishart and related computations see gjam library (for instance gjamSimData function), 
  # it is done using Rcpp 
  
  # r = number of latent factors, user chosen
  #eta_h = rep ( 0, times=r )
  
  eta_h <- 1/rgamma(r, shape = 1/2, rate = 1/(10^4)) # R
  
  Dz <- .riwish(2 + r - 1, 4 * 1/eta_h * diag(r)) # RCPP
  
  #Z <- rmvnorm ( n = N_stick, mean = rep(0, times = r), sigma = Dz ) 
  
  #mu_Zj <- rep(0, times = r)
  #for (j in seq(1, N)) {
  #  Z[j,] <- rmvnorm ( n = 1, mean = mu_Zj, sigma = Dz ) 
  #}
  
  W <- rmvnormRcpp ( n = n, mu = rep(0, times = r), sigma = diag(r) ) 
  # need to transpose for the computation of Bx_i + Q(k)Zw_i'
  
  nu = 0
  G = 1
  sigmaeps2 <- 1e-2
  
  muBeta = rep ( 0, times=n_cov ) # Prior mean of the beta coefficients
  sigmaBeta = diag ( n_cov ) # Prior variance-covariance matrix of the beta coefficients
  B <- matrix(data = 0, nrow = S, ncol = n_cov) # coefficient matrix
  
  for (j in seq(1,S,1)) {
    B[j,] <- rmvnormRcpp ( n = 1, mu = muBeta, sigma = 10*sigmaBeta )
  }
  sigmaB = 10
  
  cardinality_S = rep ( 0, times=N_stick )
  
  # Dirichlet process
  
  # each label belongs to {1,...,N_stick}
  k = c();
  
  p <- numeric(N_stick)
  pl <- matrix(0, nrow = S, ncol = N_stick)
  logpl <- matrix(0, nrow = S, ncol = N_stick)
  
  xi <- compute_GD_prior(N_stick,alpha0,k)
  # c <- 2
  # given the vector k = [1,..., 1] of length = S
  p[1] <- xi[1]
  p[2:N_stick-1] <- sapply(2:N_stick, function(j) xi[j] * prod(1 - xi[1:(j-1)]))
  p[N_stick] <- 1 - sum( p[1:N_stick-1] )
  
  k = sample(1:N_stick, size = S, replace=TRUE, prob = p)
  
  names(k) <- seq(1,S)
  
  V = Y
  V_star = matrix(0,nrow = nrow(Y),ncol = ncol(Y))
  
  S_Dz = matrix(0, nrow = r, ncol = r)
  toc()
  
  tic("GJAM R Gibbs Sampler:")
  # start of the Gibbs sampler
  for ( niter in 1:(posteriorDraws + burnInIterations) ) { # MCMC loop
    
    # Step 1 : resample the cluster assignments
    #for each species j = 1,...,S do
    Z <- matrix(data = 0, nrow = N_stick, ncol = r) 
    Q <- matrix(data = 0, nrow = S, ncol = N_stick) 
    
    # Step 3 -> RCPP
    
    # OUR CODE 
    # for ( l in 1:S ) { # Loop on species
    #   for (j in 1:N_stick) {
    #     pl[l,j] <- p[j] * exp( -1/(2*sigmaeps2) * norm_vec(V[,l] - x %*% B[l,] - W %*% Z[j,])^2)
    #   }
    #   
    #   logpl[l,] = log(pl[l,])
    #   logpl[l,] = logpl[l,] - max(logpl[l,])
    #   pl[l,] = exp(pl[l,])
    #   pl[l,] = pl[l,] / sum(pl[l,])
    #   
    #   k[l] = sample(N_stick, size = 1, replace=TRUE, prob = pl[l,])
    # }
    # BUT this is done in Rcpp by the function getPmatKRcpp in RcppExports, which computes the matrix pl
    pl <- getPmatKRcpp(pveck = p,Yk = V, Zk = Z,
                       Xk = x, Bk = B, Wk = W,
                       sigmasqk = sigmaeps2)
    # and then the vector of labels k
    k <- unlist( apply(pl, 1, function(x)sample(1:N_stick, size=1, prob=x)) )
    
    # Step 4 -> R
    # RESAMPLE p|k
    # OUR CODE
    # xi <- compute_GD_prior(N_stick,alpha0,k)
    # # c <- 2
    # # given the vector k = [1,..., 1] of length = S
    # p[1] <- xi[1]
    # p[2:N_stick-1] <- sapply(2:N_stick, function(j) xi[j] * prod(1 - xi[1:(j-1)]))
    # p[N_stick] <- 1 - sum( p[1:N_stick-1] )
    # BUT in Rcpp:
    p <- .sampleP(N = N_stick, avec = rep(alpha0/N_stick,(N_stick-1)),
                  bvec = ((N_stick-1):1)*alpha0/N_stick, K = k) 
    
    #cardinality_S = rep(0, times = N_stick)
    
    #names(cardinality_S) = seq(1,N_stick)
    
    # Step 1
    
    nn   <- nrow(x)
    pp    <- ncol(x)
    SS    <- ncol(V)
    ntot <- nrow(V)
    
    covR <- solveRcpp( (1/sigmaeps2)*crossprod(Z[k,]) + diag(r) ) # Sigma_W
    z1   <- crossprod( Z[k,]/sigmaeps2,t(Y - x %*% t(B) ))
    RR   <- rmvnormRcpp(ntot, mu = rep(0,r), sigma = covR ) + t(crossprod( covR,z1))
    
    D    <- .riwish(df = (2 + r + N_stick - 1), S = (crossprod(Z) + 2*2*diag(1/eta_h)))
    Z    <- fnZRcpp(kk=k, Yk=V, Xk=x, Dk=Dz, Bk=B, 
                    Wk=RR, sigmasqk=sigmaeps2, Nz=N_stick)
    
    A = Z[k,]
    
    # Our CODE in R 
    # for ( j in 1:N_stick ) { 
    #   
    #   # Step 1 : TO CHECK, and CAREFULLY !!!!!!!!!
    #   if (!(j %in% k)) {
    #     #print('entered if')
    #     Z[j,] = rmvnormRcpp ( n = 1, mu = rep(0, times = r), sigma = Dz )
    #     cardinality_S[j] = 1
    #   } 
    #   else { # if j in k
    #     #print('entered else')
    #     cardinality_S[j] = sum(j==k)
    #     #print(paste('cardinality',cardinality_S[j]))
    #     Sigma_Zj = solveRcpp(cardinality_S[j]/sigmaeps2 * t(W) %*% W + solveRcpp(Dz) )
    #     mu_Zj = c(0)
    #     for (l in 1:S) {
    #       if (k[l] == j){
    #         #print('entered else-if')
    #         mu_Zj = mu_Zj + 1/sigmaeps2 * Sigma_Zj %*% t(W) %*% ( t(t(V[,l])) - x %*% B[l,] )
    #       }
    #       Q[l,k[l]] = 1
    #     }
    #     # Sigma_Zj <- Sigma_Zj + diag(ncol(Sigma_Zj))*0.01
    #     Sigma_Zj <- make.positive.definite(Sigma_Zj, tol=1e-3)
    #     # because of machine precision the matrix could seem to be not positive definite:
    #     # we overcome this by making sure your det(Sigma_Zj) returns a positive 
    #     # value. One way is to add some variance in all directions:
    #     Sigma_Zj[lower.tri(Sigma_Zj)] = t(Sigma_Zj)[lower.tri(Sigma_Zj)]
    #     # because of machine precision the matrix doesn't seem symmetric
    #     Z[j,] = rmvnormRcpp ( n = 1, mu = mu_Zj, sigma = Sigma_Zj )
    #   }
    # }
    # A = Q %*% Z # in R, todo in Rcpp?
    
    # Step 2 : TO CHECK, and CAREFULLY !!!!!!!!!
    for ( i in 1:n ) {  # cycle in R, todo in Rcpp?
      Sigma_W = solveRcpp(1/sigmaeps2 * t(A) %*% A + diag(r))
      mu_W = 1/sigmaeps2 * Sigma_W %*% t(A) %*% (t(t(V[i,])) - B %*% x[i,])
      #Sigma_W <- Sigma_W + diag(ncol(Sigma_W))*0.01
      Sigma_W <- make.positive.definite(Sigma_W, tol=1e-3)
      # make.positive.definite finds the closest positive definite matrix using the
      # algorithm of NJ Higham
      # we can find an implemented version in Rcpp here:
      # https://stackoverflow.com/questions/51490499/results-for-calculating-nearest-positive-definite-matrix-are-different-in-r-func
      Sigma_W[lower.tri(Sigma_W)] = t(Sigma_W)[lower.tri(Sigma_W)]
      W[i,] = rmvnormRcpp ( n = 1, mu = mu_W, sigma = Sigma_W )
    }
    
    # Step 5 in R, todo in Rcpp
    
    # nn   <- nrow(x)
    # pp    <- ncol(x)
    # SS    <- ncol(V)
    # ntot <- nrow(V)
    # 
    # covR <- solveRcpp( (1/sigmaeps2)*crossprod(Z[k,]) + diag(r) ) # Sigma_W
    # z1   <- crossprod( Z[k,]/sigmaeps2,t(Y - x %*% t(B) ))        
    # RR   <- rmvnormRcpp(ntot, mu = rep(0,r), sigma = covR ) + t(crossprod( covR,z1))
    # 
    # rndEff <- RR%*%t(Z[k,])
    # 
    # res        <- sum((V - x %*% t(B) - rndEff )^2)
    # sigmaeps2 <- 1/rgamma(1,shape=(SS*nn + 1)/2, rate=res/2)  
    # print(sigmaeps2)
    
    dx = c(0)
    for (i in 1:n) {
      dx = dx + norm_vec(V[i,] - B %*% x[i,] - A %*% W[i,])^2
    }
    sigmaeps2 = 1/rgamma(1,shape = (n * S + nu)/2, rate = dx/2 + nu/G^2)
    #sigmaeps2 = 1
    
    # Step 6 in R
    eta_h <- 1/rgamma(r, shape = (2 + r )/2, 
                      rate = ((1/1000000) + 2*diag(solveRcpp(Dz)) ) )
    
    S_Dz <- crossprod(Z) + 4 * diag(1/eta_h)
    # S_Dz <- S_Dz + diag(ncol(S_Dz))*0.01
    #S_Dz <- make.positive.definite(S_Dz, tol=1e-3)
    #S_Dz[lower.tri(S_Dz)] = t(S_Dz)[lower.tri(S_Dz)]
    Dz = .riwish(df = 2 + r + N_stick - 1, S = S_Dz)
    
    # 
    # Step 7 : TO CHECK, and REALLY CAREFULLY !!!!!!!!! in R, todo Rcpp?
    Sigma_star = A %*% t(A) + sigmaeps2 * diag(S)
    D = diag(diag(Sigma_star))
    R = solveRcpp(D)^(1/2) %*% Sigma_star %*% solveRcpp(D)^(1/2)
    B_star = solveRcpp(D)^(1/2) %*% B
    
    for (i in seq(1,n)) {
      for (j in seq(1,S)) {
        
        mean1 = as.vector(B_star[j,] %*% x[i,] + A[j,] %*% W[i,])
        
        if (Y[i,j] == 1){
          V_star[i,j] <- .tnorm(1, 0, Inf, mean1, sigmaeps2)
        } else {
          V_star[i,j] <- .tnorm(1, -Inf, 0, mean1, sigmaeps2)
        }
      }
    }
    
    V = V_star %*% solveRcpp(D)^(1/2)
    
    
    # L<-x%*%t(B) #We create the mean by multiplying B with the design matrix X
    # 
    # Sigma_star <-A%*%t(A)+sigmaeps2*diag(S) #We obtain Sigma. Here sigma_epsilon^2 is 0.1
    # 
    # 
    # #We obtain the correlation matrix R
    # B_star =cov2cor(Sigma_star) 
    # V<-L+rmvnormRcpp(n = n, mu=rep(0,S), sigma=B_star)
    
    # Step 8 : TO CHECK, and REALLY CAREFULLY !!!!!!!!! to do in Rcpp?
    for (j in seq(1,S,1)) {
      muBetaj = solveRcpp(1/(sigmaB)^2 * diag(n_cov) + 1/sigmaeps2 * t(x) %*% x) %*% t(x) %*% (V_star[,j] - W %*% A[j,]) * 1/sigmaeps2
      sigmaBetaj = solveRcpp(1/(sigmaB^2) * diag(n_cov) + 1/sigmaeps2 * t(x) %*% x)
      # sigmaBetaj <- sigmaBetaj + diag(ncol(sigmaBetaj))*0.01
      sigmaBetaj <- make.positive.definite(sigmaBetaj, tol=1e-3)
      sigmaBetaj[lower.tri(sigmaBetaj)] = t(sigmaBetaj)[lower.tri(sigmaBetaj)]
      B_star[j,] <- rmvnormRcpp ( n = 1, mu = muBetaj, sigma = sigmaBetaj )
    }
    
    B = solveRcpp(D)^(1/2) %*% B_star; 
    
    list_B[[niter]] <- B
    list_A[[niter]] <- A
    list_Z[[niter]] <- Z
    list_k[[niter]] <- k
    list_R[[niter]] <- R
    list_sigmaeps2[[niter]] <- sigmaeps2
  } 
  # end of Gibbs sampler
  toc()
  return(list(B = list_B, A = list_A, Z = list_Z, k = list_k, R = list_R, sigmaeps2 = list_sigmaeps2))
}

.tnorm <- function(n,lo,hi,mu,sig){   
  
  #normal truncated lo and hi
  
  tiny <- 10e-6
  
  if(length(lo) == 1 & length(mu) > 1)lo <- rep(lo,length(mu))
  if(length(hi) == 1 & length(mu) > 1)hi <- rep(hi,length(mu))
  
  q1 <- pnorm(lo,mu,sig)
  q2 <- pnorm(hi,mu,sig) 
  
  z <- runif(n,q1,q2)
  z <- qnorm(z,mu,sig)
  
  z[z == Inf]  <- lo[z == Inf] + tiny
  z[z == -Inf] <- hi[z == -Inf] - tiny
  z
}

compute_GD_prior <- function(N_stick, alpha0, k){
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