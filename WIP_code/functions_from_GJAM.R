# USEFUL FUNCTIONS FROM gjamHfunctions.R
.sampleP <- function(N, avec, bvec, K){
  
  a    <- avec + vapply(1:(N-1), function(k)sum(K == k), 0)
  b    <- bvec + vapply(1:(N-1), function(k)sum(K > k), 0)
  V    <- rbeta((N - 1), a, b)
  p    <- vector("numeric",length=N)
  p[1] <- V[1]
  for(l in 2:(N - 1))p[l] <- prod(1 - V[1:(l - 1)])*V[l]
  p[N] <- prod(1 - V)   
  p
}

.getPars <- function(CLUST, x, N, r, Y, B, D, Z, sigmaerror, K, pvec,
                     alpha.DP, inSamples,...){      
  
  # Y includes all terms but x%*%beta
  
  nn   <- length(inSamples)
  p    <- ncol(x)
  S    <- ncol(Y)
  ntot <- nrow(Y)
  nn   <- length(inSamples)
  
  covR <- solveRcpp( (1/sigmaerror)*crossprod(Z[K,]) + diag(r) ) # Sigma_W
  z1   <- crossprod( Z[K,]/sigmaerror,t(Y - x%*%t(B)) )        
  RR   <- rmvnormRcpp(ntot, mu = rep(0,r), sigma = covR ) + t(crossprod( covR,z1))
  if(nn < ntot)RR[-inSamples,] <- rmvnormRcpp(ntot-nn,mu=rep(0,r), sigma=diag(r))
  rndEff <- RR%*%t(Z[K,])
  
  res        <- sum((Y[inSamples,] - x[inSamples,]%*%t(B) - rndEff[inSamples,] )^2)
  sigmaerror <- 1/rgamma(1,shape=(S*nn + 1)/2, rate=res/2)  
  
  if(CLUST){   #only until convergence
    avec <- 1/rgamma(r, shape = (2 + r )/2, 
                     rate = ((1/1000000) + 2*diag(solveRcpp(D)) ) )  
    
    D    <- .riwish(df = (2 + r + N - 1), S = (crossprod(Z) + 2*2*diag(1/avec)))
    Z    <- fnZRcpp(kk=K, Yk=Y[inSamples,], Xk=x[inSamples,], Dk=D, Bk=B, 
                    Wk=RR[inSamples,], sigmasqk=sigmaerror, Nz=N)
    
    pmat <- getPmatKRcpp(pveck = pvec,Yk = Y[inSamples,], Zk = Z,
                         Xk = x[inSamples,], Bk = B, Wk = RR[inSamples,],
                         sigmasqk = sigmaerror)
    K    <- unlist( apply(pmat, 1, function(x)sample(1:N, size=1, prob=x)) )
    pvec <- .sampleP(N = N, avec = rep(alpha.DP/N,(N-1)),
                     bvec = ((N-1):1)*alpha.DP/N, K = K)  
  }
  
  list(A = Z[K,], D = D, Z = Z, K = K, pvec = pvec, 
       sigmaerror = sigmaerror, rndEff = rndEff)
} 

.riwish <- function(df,S){
  solveRcpp(.rwish(df,solveRcpp(S)))
}

#IMPORTANT COMMENT ON IMPORTANT LINES
#line 4006 - gjam()
#line 4780 - starts loop for gibbs sampling that ends at line 5323
#line 5217 - if (g>burnin) ->here Taylor does important things (I think)
rmvnormRcpp

updateBeta (that calls) betaWrapper

.wWrapper

tnormMVNmatrix

.xpredSetup

.setupFactors #function that update cluster size (I think), used in line 4615 in gjamHfunctions.R

solveRcpp

.updateW (that calls) .wWrapper

.contrastCoeff

.dMVN #mvn density for mean 0 
dmvnormRcpp

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
.tnormMVNmatrix <- function(avec, muvec, smat, 
                            lo=matrix(-1000,nrow(muvec),ncol(muvec)), 
                            hi=matrix(1000,nrow(muvec),ncol(muvec)),
                            whichSample = c(1:nrow(smat))){
  
  #lo, hi must be same dimensions as muvec,avec
  
  lo[lo < -1000] <- -1000
  hi[hi > 1000]  <- 1000
  
  if(max(whichSample) > length(muvec))
    stop('whichSample outside length(muvec)')
  
  r <- avec
  a <- trMVNmatrixRcpp(avec, muvec, smat, lo, hi, whichSample, 
                       idxALL = c(0:(nrow(smat)-1)) )  
  r[,whichSample] <- a[,whichSample]
  r
}
