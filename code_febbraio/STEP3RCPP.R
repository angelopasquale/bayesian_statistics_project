STEP3RCPP<-function(p,V,Z,x,B,W,sigmaeps2,N_stick){
  ###STEP 3###
  #This function computes vector of labels a posteriori
  pl <- getPmatKRcpp(pveck = p,Yk = V, Zk = Z,
                     Xk = x, Bk = B, Wk = W,
                     sigmasqk = sigmaeps2)
  # and then the vector of labels k
  k <- unlist( apply(pl, 1, function(x)sample(1:N_stick, size=1, prob=x)) )
  k
}