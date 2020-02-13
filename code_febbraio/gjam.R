#source(file = "gjamHfunctions.R")
gjamReduct <- function(formula, xdata, ydata, modelList){
  
  holdoutN      <-  0
  holdoutIndex  <- numeric(0)
  modelSummary  <- betaPrior  <- traitList <- effort <- NULL
  specByTrait   <- traitTypes <- breakList <- notStandard <- NULL
  censor <- censorCA <- censorDA <- CCgroups <- FCgroups <- intMat <- NULL
  reductList <- y0 <- N  <- r <- otherpar <- pg <- NULL
  ng     <- 2000
  burnin <- 500
  REDUCT <- TRAITS <- FULL <- F
  PREDICTX <- T
  lambdaPrior <- betaPrior <- NULL
  
  RANDOM <- F              # random group intercepts
  
  TIME <- F
  
  ematAlpha <- .5
  
  alpha.DP <- ncol(ydata)          # large values give more variation
  
  if(alpha.DP == 1)
    stop('multivariate model: at least 2 columns needed in ydata')
  
  for(k in 1:length(modelList))assign( names(modelList)[k], modelList[[k]] )
  
  if('CCgroups' %in% names(modelList))attr(typeNames,'CCgroups')  <- CCgroups
  if('FCgroups' %in% names(modelList))attr(typeNames,'FCgroups')  <- FCgroups
  if('CATgroups' %in% names(modelList))attr(typeNames,'CATgroups') <- CATgroups
  
  if(burnin >= ng) stop( 'burnin must be < no. MCMC steps, ng' )
  if('censor' %in% names(modelList)){
    for(k in 1:length(censor)){
      if( nrow(censor[[k]]$partition) != 3 )
        stop('censor matrix: 3 rows for value, lo, hi')
      rownames(censor[[k]]$partition) <- c('value','lo','hi')
    }
  }
  
  S <- ncol(ydata)
  if(length(typeNames) == 1)typeNames <- rep(typeNames,S)
  if(length(typeNames) != S) 
    stop('typeNames must be one value or no. columns in y')
  
  ############### factors in y
  
  tmp <- .checkYfactor(ydata, typeNames)
  ydata <- tmp$ydata; yordNames <- tmp$yordNames
  
  tmp <- .buildYdata(ydata, typeNames)
  y   <- tmp$y
  ydataNames <- tmp$ydataNames
  typeNames  <- tmp$typeNames
  CCgroups   <- tmp$CCgroups
  FCgroups   <- tmp$FCgroups
  CATgroups  <- tmp$CATgroups
  
  S <- ncol(y)
  n <- nrow(y)
  
  cat("\nObservations and responses:\n")
  print(c(n, S))
  
  tmp    <- .buildEffort(y, effort, typeNames)
  effort <- tmp
  effMat <- effort$values
  modelList$effort <- effort
  re <- floor( diff( range(log10(effMat),na.rm=T) ) )
  if(re > 2)
    message(paste('sample effort > ', re, ' orders of magnitude--consider units near 1',sep='') )
  
  
  tmp      <- .gjamGetTypes(typeNames)
  typeCols <- tmp$typeCols
  typeFull <- tmp$typeFull
  typeCode <- tmp$TYPES[typeCols]
  allTypes <- sort(unique(typeCols))
  
  tmp <- .gjamXY(formula, xdata, y, typeNames, notStandard)
  x      <- tmp$x; y <- tmp$y; snames <- tmp$snames
  xdata  <- tmp$xdata; xnames <- tmp$xnames
  interBeta   <- tmp$interaction 
  factorBeta  <- tmp$factorAll
  designTable <- tmp$designTable;    xscale <- tmp$xscale
  predXcols   <- tmp$predXcols
  standMat    <- tmp$standMat;      standMu <- tmp$standMu  
  standRows   <- tmp$standRows;    
  xdataNames  <- tmp$xdataNames
  notStandard <- tmp$notStandard[tmp$notStandard %in% xnames]
  
  factorLambda <- interLambda <- NULL
  
  if(!is.null(lambdaPrior)){
    
    lformula <- attr(lambdaPrior$lo,'formula')
    
    tmp <- .gjamXY(lformula, xdata, y, typeNames, notStandard)
    xl   <- tmp$x
    mm <- match(colnames(xl),colnames(xdata))
    wm <- which(is.finite(mm))
    if(length(wm) > 0){
      xdata[,mm[wm]] <- xl[,wm]
    }
    
    xlnames <- tmp$xnames
    interLambda   <- tmp$interaction
    factorLambda <- tmp$factorAll
    
    designTable <- list(beta = designTable, lambda = tmp$designTable)
    
    standMatL    <- tmp$standMat;      standMuL <- tmp$standMu
    standRowsL    <- tmp$standRows;    
    notStandardL <- tmp$notStandard[tmp$notStandard %in% xlnames]
  }
  
  modelList     <- append(modelList, list('formula' = formula,
                                          'notStandard' = notStandard))
  
  Q <- ncol(x)
  
  tmp <- .gjamMissingValues(x, y, factorBeta$factorList, typeNames)
  xmiss  <- tmp$xmiss;   xbound <- tmp$xbound; 
  ymiss  <- tmp$ymiss;   missY <- tmp$missY
  xprior <- tmp$xprior;  yprior <- tmp$yprior
  nmiss  <- nrow(xmiss); mmiss  <- nrow(ymiss)
  x  <- tmp$x; y <- tmp$y
  
  reductList <- .setupReduct(modelList, S, Q, n) ##########
  N <- reductList$N; r <- reductList$r
  if(!is.null(reductList$N))REDUCT <- T
  
  
  tmp <- .gjamHoldoutSetup(holdoutIndex, holdoutN, n)
  holdoutIndex <- tmp$holdoutIndex; holdoutN <- tmp$holdoutN
  inSamples    <- tmp$inSamples;         nIn <- tmp$nIn
  
  tmp <- .gjamSetup(typeNames, x, y, breakList, holdoutN, holdoutIndex,
                    censor=censor, effort=effort) 
  w <- tmp$w; z <- tmp$z; y <- tmp$y; other <- tmp$other; cuts <- tmp$cuts
  cutLo       <- tmp$cutLo; cutHi <- tmp$cutHi; plo <- tmp$plo; phi <- tmp$phi
  ordCols     <- tmp$ordCols; disCols <- tmp$disCols; compCols <- tmp$compCols 
  conCols     <- which(typeNames == 'CON')
  classBySpec <- tmp$classBySpec; breakMat <- tmp$breakMat
  minOrd      <- tmp$minOrd; maxOrd <- tmp$maxOrd; censorCA <- tmp$censorCA
  censorDA    <- tmp$censorDA; censorCON <- tmp$censorCON; 
  ncut <- ncol(cuts);  corCols <- tmp$corCols
  catCols     <- which(attr(typeNames,'CATgroups') > 0)
  sampleW     <- tmp$sampleW
  ordShift    <- tmp$ordShift
  
  sampleW[censorCA] <- 1
  sampleW[censorDA] <- 1
  sampleW[censorCON] <- 1
  sampleWhold <- tgHold <- NULL
  wHold <- NULL
  wmax  <- apply(y/effMat,2,max,na.rm=T)
  pmin  <- -2*abs(wmax)
  
  ploHold <- phiHold <- NULL
  
  if(holdoutN > 0){
    sampleWhold <- sampleW[holdoutIndex,]  #to predict X
    sampleW[holdoutIndex,] <- 1
    tgHold  <- cuts
    wHold   <- w[drop=F,holdoutIndex,]
    
    ploHold <- plo[drop=F,holdoutIndex,]   # if LOHI: updated to current yp
    phiHold <- phi[drop=F,holdoutIndex,]
  }
  
  byCol <- byRow <- F
  if(attr(sampleW,'type') == 'cols')byCol <- T
  if(attr(sampleW,'type') == 'rows')byRow <- T
  indexW <- attr(sampleW,'index')
  
  notCorCols <- c(1:S)
  if(length(corCols) > 0)notCorCols <- notCorCols[-corCols]
  
  ############ 'other' columns
  sigmaDf  <- nIn - Q + S - 1
  sg <- diag(.1,S)
  SO <- S
  
  notOther <- c(1:S)
  sgOther  <- NULL
  if(length(other) > 0){                     
    notOther   <- notOther[!notOther %in% other]
    SO         <- length(notOther)
    sg[other,] <- sg[,other] <- 0
    sgOther    <- matrix( cbind(other,other),ncol=2 )
    sg[sgOther] <- .1
  }
  
  ############## prior on beta
  loB <- hiB <- NULL
  beta <- bg <- matrix(0,Q,S)
  rownames(beta) <- colnames(x)
  BPRIOR <- F
  
  if( !is.null(betaPrior) ){
    colnames(betaPrior$lo) <- .cleanNames(colnames(betaPrior$lo))
    colnames(betaPrior$hi) <- .cleanNames(colnames(betaPrior$hi))
    loB <- betaPrior$lo
    hiB <- betaPrior$hi
    
    bg <- (loB + hiB)/2
    bg[is.nan(bg)] <- 0
    
    wB <- which(!is.na(t(loB[,notOther])), arr.ind=T)[,c(2,1)]
    wB <- rbind(wB, which(!is.na(t(hiB[,notOther])), arr.ind=T)[,c(2,1)])
    colnames(wB) <- c('row','col')
    
    tmp <- .betaPrior(bg, notOther, loB, hiB)
    bg <- tmp$beta; loB <- tmp$loB; hiB <- tmp$hiB
    wB <- tmp$wB; BPRIOR <- tmp$BPRIOR
    bg[is.nan(bg)] <- 0
    
    tmp <- .getPattern(bg[,notOther], wB)
    Brows <- tmp$rows
    Bpattern <- tmp$pattern
    BPRIOR <- T
    bg[!is.finite(bg)] <- 0
  }
  
  zeroBeta <- .factorCoeffs2Zero(factorBeta, snames, betaPrior)  # max zero is missing factor level
  zeroLambda <- NULL
  
  if(byCol){
    inw <- intersect( colnames(y)[indexW], colnames(y)[notOther] )
    indexW <- match(inw,colnames(y)[notOther])
  }
  
  IXX <- NULL
  
  XX    <- crossprod(x)
  IXX <- chol2inv(chol( XX ) )
  
  
  updateBeta <- .betaWrapper(REDUCT, TIME, BPRIOR, notOther, IXX, 
                             betaLim=max(wmax)/2)
  
  ############ dimension reduction
  
  inSamp <- inSamples
  
  CLUST <- T   # dirichlet 
  
  .param.fn <- .paramWrapper(REDUCT, inSamp, SS=length(notOther))
  
  # this function calls .getPars, which in turn calls: 
  # - fnZRcpp containing the full conditional of Z
  # - getPmatKRcpp containing the full conditional of the matrix pl
  # - .sampleP for sampling p|k
  # sigmaerror Inv-Gamma
  # D Inv-Wishart
  # Steps 1 - 3 - 4 - 5- 6
  
  sigmaerror <- .1
  otherpar   <- list(S = S, Q = Q, sigmaerror = sigmaerror, 
                     Z = NA, K =rep(1,S), sigmaDf = sigmaDf)
  sigErrGibbs <- rndEff <- NULL
  
  yp <- y
  wmax <- ymax <- apply(y,2,max)
  wmax <- wmax/effMat
  
  if(REDUCT){
    cat( paste('\nDimension reduced from',S,'X',S,'->',N,'X',r,'responses\n') )
    otherpar$N <- N; otherpar$r <- r; otherpar$sigmaerror <- 0.1
    otherpar$Z <- rmvnormRcpp(N,rep(0,r),1/S*diag(r))
    otherpar$D <- .riwish(df = (2 + r + N), 
                          S = (crossprod(otherpar$Z) +
                                 2*2*diag(rgamma(r,shape=1,rate=0.001))))
    otherpar$K <- sample(1:N,length(notOther),replace=T)
    
    otherpar$alpha.DP <- alpha.DP
    otherpar$pvec     <- .sampleP(N=N, avec=rep(alpha.DP/N,(N-1)),
                                  bvec=((N-1):1)*alpha.DP/N, K=otherpar$K)
    kgibbs <- matrix(1,ng,S)
    sgibbs <- matrix(0,ng, N*r)
    nnames <- paste('N',1:N,sep='-')
    rnames <- paste('r',1:r,sep='-')
    colnames(sgibbs) <- .multivarChainNames(nnames,rnames)
    sigErrGibbs <- rep(0,ng)   
    
    rndEff <- w*0
    
  } else {
    Kindex <- which(as.vector(lower.tri(diag(S),diag=T)))
    nK     <- length(Kindex)
    sgibbs <- matrix(0,ng,nK)
    colnames(sgibbs) <- .multivarChainNames(snames,snames)[Kindex] # half matrix
  }
  
  out <- .param.fn(CLUST=T, x, beta = bg[,notOther], Y = w[,notOther], otherpar)  
  sg[notOther,notOther]    <- out$sg
  otherpar      <- out$otherpar
  
  muw <- w
  
  Y <- w[inSamp,notOther]
  sig <- sg[notOther,notOther]
  
  if(REDUCT){
    Y <- Y - rndEff[inSamp,notOther]
    sig <- sigmaerror
  }
  
  bg[,notOther] <- updateBeta(X = x[inSamp,], Y, sig, beta = bg[,notOther],
                              loB, hiB)
  muw <- x%*%bg
  
  sg[other,] <- sg[,other] <- 0
  diag(sg)[other]          <- 1
  rownames(bg)  <- xnames
  rownames(sg)  <- colnames(sg) <- colnames(bg) <- snames
  colnames(x)   <- xnames
  
  
  
  ############ ordinal data
  
  cutg <- tg <- numeric(0)
  
  ############ setup w
  tmp <- .gjamGetTypes(typeNames)
  typeFull <- tmp$typeFull
  typeCols <- tmp$typeCols
  allTypes <- unique(typeCols)
  Y <- w
  
  LOHI <- F
  if(!LOHI & holdoutN > 0){
    minlo <- apply(plo,2,min)
    minlo[minlo > 0] <- 0
    maxhi <- apply(phi,2,max)
  }
  
  .updateW <- .wWrapper(REDUCT, RANDOM, S, effMat, corCols, notCorCols, typeNames, 
                          typeFull, typeCols, 
                          allTypes, holdoutN, holdoutIndex, censor, 
                          censorCA, censorDA, censorCON, notOther, sampleW, 
                          byRow, byCol,
                          indexW, ploHold, phiHold, sampleWhold, inSamp)
  
  ycount <- rowSums(y)
  
  
  ############ X prediction
  
  tmp <- .xpredSetup(Y, x, bg, interBeta$isNonLinX, factorBeta, 
                     factorBeta$intMat, 
                     standMat, standMu, notOther, notStandard) 
  factorBeta$linFactor <- tmp$linFactor; xpred <- tmp$xpred; px <- tmp$px
  lox <- tmp$lox; hix <- tmp$hix
  
  priorXIV  <- diag(1e-5,ncol(x))
  priorX    <- colMeans(x)
  priorX[abs(priorX) < 1e-10] <- 0
  
  linFactor <- NULL
  
  ################## random groups
  
  if('random' %in% names(modelList)){
    
    RANDOM <- T
    rname  <- modelList$random
    randGroupTab <- table( as.character(xdata[,rname]) )
    
    wss <- names(randGroupTab[randGroupTab <= 2])
    if(length(wss) > 0){
      xdata[,rname] <- .combineFacLevels(xdata[,rname], fname=wss, 
                                         aname = 'rareGroups', vminF=1)
      randGroupTab <- table( as.character(xdata[,rname]) )
    }
    
    randGroups <- names( randGroupTab )
    G <- length(randGroups)
    
    groupIndex  <- match(as.character(xdata[,rname]),randGroups)
    rmm <- matrix(groupIndex,length(groupIndex), S)
    smm <- matrix(1:S, length(groupIndex), S, byrow=T)
    
    randGroupIndex <- cbind( as.vector(smm), as.vector(rmm) )
    colnames(randGroupIndex) <- c('species','group')
    xdata[,rname] <- as.factor(xdata[,rname])
    alphaRandGroup <- matrix(0, S, G)
    rownames(alphaRandGroup) <- snames
    colnames(alphaRandGroup) <- randGroups
    Cmat <- var(w[,notOther]/2)
    Cmat <- Cmat + diag(.1*diag(Cmat))
    Cprior <- Cmat
    CImat <- solve(Cprior)
    Ckeep <- diag(S)
    
    alphaRanSums <- alphaRandGroup*0
    groupRandEff <- w*0
    
    Kindex <- which(as.vector(lower.tri(diag(S),diag=T)))
    nK     <- length(Kindex)
    alphaVarGibbs <- matrix(0,ng,nK)
    colnames(alphaVarGibbs) <- .multivarChainNames(snames,snames)[Kindex] # half matrix
  }
  
  
  ################################## XL prediction: variables in both
  
  
  ############  contrasts, predict F matrix
  
  tmp <- .setupFactors(xdata, xnames, factorBeta)
  ff  <- factorBeta[names(factorBeta) != 'factorList']
  factorBeta <- append(ff,tmp)
  
  ############ E matrix
  emat <- matrix(0,S,S)
  colnames(emat) <- rownames(emat) <- snames
  lo <- hi <- lm <- hm <- ess <- emat
  
  fmat <- factorBeta$fmat
  fnames <- rownames( factorBeta$lCont )
  q2 <- nrow(fmat)
  
  ############ sp richness
  richness <- richFull <- NULL
  RICHNESS <- F
  
  inRichness <- which(!typeNames %in% c('CON','CAT','OC'))
  inRichness <- inRichness[!inRichness %in% other]
  if(length(inRichness) > 2)RICHNESS  <- T
  
  wrich <- y*0 
  wrich[,inRichness] <- 1
  
  presence <- w*0
  
  covx <- cov(x)
  
  ############ sums
  predx  <- predx2 <- xpred*0
  yerror <- ypred  <- ypred2 <- wpred  <- wpred2 <- ymissPred <- ymissPred2 <- y*0
  sumDev <- 0   #for DIC
  sMean  <- sg*0
  ntot   <- 0
  
  ############ gibbs chains
  
  q2 <- length(fnames)
  fSensGibbs <- matrix(0,ng,q2)
  colnames(fSensGibbs) <- fnames
  
  bFacGibbs <- matrix(0,ng,q2*SO)
  colnames(bFacGibbs) <- .multivarChainNames(fnames,snames[notOther])
  
  bgibbs <- matrix(0,ng,S*Q)
  colnames(bgibbs) <- .multivarChainNames(xnames,snames)
  bgibbsUn <- bgibbs                   # unstandardized
  
  covE <- cov( x%*%factorBeta$dCont )  # note that x is standardized
  
  pbar <- txtProgressBar(min=1,max=ng,style=1)
  
  # unstandardize
  
  tmp <- .getUnstandX(x, standRows, standMu[,1],standMat[,1],
                      interBeta$intMat)
  S2U      <- tmp$S2U
  xUnstand <- tmp$xu
  
  if(REDUCT){
    rndTot <- w*0 
  }
  notPA <- which(!typeNames == 'PA' & !typeNames == 'CON')
  
  
  if(length(y) < 10000 | FULL) FULL <- T
  
  if(FULL){
    ygibbs <- matrix(0,ng,length(y))
  }
  if(RICHNESS){
    ypredPres <- ypredPres2 <- ypredPresN <- y*0
    shannon   <- rep(0,n)
  }
  
  for(g in 1:ng){ ########################################################
    
    if(REDUCT){
      
      #   if(g > burnin)CLUST <- F
      
      Y <- w[,notOther]
      if(RANDOM)Y <- Y - groupRandEff[,notOther]
      
      tmp <- .param.fn(CLUST=T, x, beta = bg[,notOther], Y = Y, otherpar)
      sg[notOther,notOther] <- tmp$sg
      otherpar            <- tmp$otherpar
      rndEff[,notOther]   <- tmp$rndEff
      sigmaerror          <- otherpar$sigmaerror
      kgibbs[g,notOther]  <- otherpar$K
      sgibbs[g,]          <- as.vector(otherpar$Z)
      sigErrGibbs[g]      <- sigmaerror
      
      if(length(corCols) > 0){
        if(max(diag(sg)[corCols]) > 5){  #overfitting covariance
          stop(
            paste('\noverfitted covariance, reductList$N = ',N, 
                  'reductList$r = ',r, '\nreduce N, r\n')
          )
        }
      }
      
      sg[sgOther]         <- .1*sigmaerror
      
      sinv <- .invertSigma(sg[notOther,notOther],sigmaerror,otherpar,REDUCT)
      sdg  <- sqrt(sigmaerror)
      
        
      Y <- w[inSamp,notOther] - rndEff[inSamp,notOther]
      if(RANDOM)Y <- Y - groupRandEff[inSamp,notOther]
      bg[,notOther] <- updateBeta(X = x[inSamp,], Y, 
                                  sig = sigmaerror, beta = bg[,notOther],
                                  lo=loB[,notOther], hi=hiB[,notOther])
      muw[inSamp,] <- x[inSamp,]%*%bg
      
    } else {
      Y <- w[inSamp,notOther]
      if(RANDOM)Y <- Y - groupRandEff[inSamp,notOther]
      bg[,notOther] <- updateBeta(X = x[inSamp,], Y, 
                                  sig = sg[notOther,notOther], 
                                  beta = bg[,notOther], 
                                  lo=loB, hi=hiB)
      
      muw[inSamp,] <- x[inSamp,]%*%bg
      
      SS   <- crossprod(w[inSamp,] - muw[inSamp,])
      SI   <- solveRcpp(SS[notOther,notOther])
      sinv <- .rwish(sigmaDf,SI)
      
      sg[notOther,notOther] <- solveRcpp(sinv)
      sgibbs[g,] <- sg[Kindex]
    }
    
    # muw does not include rndEff or groupRandEff
    
    alphaB <- .sqrtRootMatrix(bg,sg,DIVIDE=T)
    
    if(RANDOM){
      
      cw <- w - muw
      if(REDUCT){
        cw <- cw - rndEff
        v  <- 1/sigmaerror*.byGJAM(as.vector(cw), randGroupIndex[,1], 
                                   randGroupIndex[,2], alphaRandGroup*0, 
                                   fun='sum')[notOther,]
        sinv <- diag(1/sigmaerror, SO)
      }else{
        v <- .byGJAM(as.vector(cw), randGroupIndex[,1], 
                     randGroupIndex[,2], alphaRandGroup*0, fun='sum')[notOther,]
        v <- sinv%*%v
      }
      
      alphaRandGroup[notOther,] <- randEffRcpp(v, randGroupTab, 
                                               sinv, CImat)
      if(length(other) > 0)alphaRandGroup[other,] <- 0
      if(g < 100){
        alphaRandGroup[notOther,] <- 
          sweep( alphaRandGroup[notOther,], 2, 
                 colMeans(alphaRandGroup[notOther,]), '-')
      }
      SS  <- crossprod(t(alphaRandGroup[notOther,]))
      SS  <- S*SS + Cmat
      
      testv <- try( chol(SS) ,T)
      if( inherits(testv,'try-error') ){
        tiny  <- .1*diag(SS)
        SS  <- SS + diag(diag(SS + tiny))
      }
      
      Ckeep[notOther,notOther] <- .riwish( df = S*G + 1, SS )
      CImat <- solveRcpp(Ckeep[notOther,notOther])
      
      alphaVarGibbs[g,] <- Ckeep[Kindex]
      groupRandEff <- t(alphaRandGroup)[groupIndex,]
    }
    
    #############not TIME
      
      tmp   <- .updateW( rows=1:n, x, w, y, bg, sg, alpha=alphaB, 
                         cutg, plo, phi, rndEff, groupRandEff, 
                         sigmaerror, wHold )
      w     <- tmp$w
      yp    <- tmp$yp
      plo   <- tmp$plo
      phi   <- tmp$phi
      wHold <- tmp$wHold    #values for w if not held out
      
      
      Y <- w[,notOther]
      if(holdoutN > 0) Y[holdoutIndex,] <- wHold[,notOther]  # if w not held out
      if(RANDOM)Y <- Y - groupRandEff[,notOther]
      
      if( PREDICTX & length(predXcols) > 0){
        
        if( length(interBeta$isNonLinX) > 0 ){
          
          xpred <- .predictY2X_nonLinear(xpred, yy=Y,bb=bg[,notOther],
                                         ss=sg[notOther,notOther],
                                         priorIV = priorXIV,priorX=priorX,
                                         factorObject = factorBeta, interObject = interBeta,
                                         lox, hix)$x
        }
        
        if( length(px) > 0 ){
          wn <- which(!is.finite(xpred),arr.ind=T)
          if(length(wn) > 0){
            tmp <- matrix(priorX,Q,nrow(wn))
            xpred[wn[,1],] <- t(tmp)
          }
          xpred[,px] <- .predictY2X_linear(xpred, yy=Y, bb=bg[,notOther],
                                           ss=sg[notOther,notOther], sinv = sinv,
                                           priorIV = priorXIV, 
                                           priorX=priorX,predCols=px, 
                                           REDUCT=REDUCT, lox, hix)[,px]
          wn <- which(!is.finite(xpred),arr.ind=T)
          if(length(wn) > 0){
            tmp <- matrix(priorX,Q,nrow(wn))
            xpred[wn[,1],] <- t(tmp)
          }
        }
        
        if( length(factorBeta$linFactor) > 0 ){
          
          # predict all factors
          xtmp <- xpred
          xtmp[,factorBeta$findex] <- 
            .predictY2X_linear(xpred, yy=Y, 
                               bb=bg[,notOther],
                               ss=sg[notOther,notOther], sinv = sinv,
                               priorIV = priorXIV, 
                               priorX=priorX,predCols=factorBeta$findex, 
                               REDUCT=REDUCT, lox, hix)[,factorBeta$findex]
          for(k in 1:length(factorBeta$linFactor)){
            
            mm  <- factorBeta$linFactor[[k]]
            tmp <- xtmp[,mm]
            
            tmp[,1] <- 0
            ix  <- apply(tmp,1,which.max)   
            
            tmp <- tmp*0
            tmp[ cbind(1:n,ix) ] <- 1
            tmp <- tmp[,-1,drop=F]
            xpred[,mm[-1]] <- tmp
          }
        }
        xpred[,1] <- 1
      }
    
    setTxtProgressBar(pbar,g)
    
    bgu <- bg                    # unstandardize beta
    if(length(standRows) > 0){
      bgu <- S2U%*%x%*%bg
    }
    
    bgibbsUn[g,] <- bgu          # unstandardized
    bgibbs[g,]   <- bg           # standardized
    
    # Fmatrix centered for factors, 
    # bg is standardized by x, bgu is unstandardized
    
    tmp <- .contrastCoeff(beta=bg[,notOther], 
                          notStand = notStandard[notStandard %in% xnames], 
                          sigma = sg[notOther,notOther], sinv = sinv,
                          stand = standMat, factorObject=factorBeta )
    agg   <- tmp$ag
    beg   <- tmp$eg
    fsens <- tmp$sens
    
    fSensGibbs[g,]  <- sqrt(diag(fsens))
    bFacGibbs[g,] <- agg       # stand for X and W, centered for factors
    
    if(FULL)ygibbs[g,] <- as.vector(yp)
    
    if(g > burnin){
      
      ntot   <- ntot + 1
      ypred  <- ypred + yp
      ypred2 <- ypred2 + yp^2
      
      tmp <- .dMVN(w[,notOther], muw[,notOther], sg[notOther,notOther], log=T)
      
      sumDev <- sumDev - 2*sum(tmp) 
      yerror <- yerror + (yp - y)^2
      
      fmat <- fmat + fsens
      
      sMean  <- sMean + sg
      
      wpred  <- wpred + w
      wpred2 <- wpred2 + w^2
      
      if(RICHNESS){
        
        yy <- yp
        
        if('PA' %in% typeNames){
          wpa <- which(typeNames[inRichness] == 'PA')
          yy[,inRichness[wpa]] <- round(yp[,inRichness[wpa]]) #######
        }
        
        if(length(notPA) > 0){
          w0 <- which(yy[,notPA] <= 0)
          w1 <- which(yy[,notPA] > 0)
          yy[,notPA][w0] <- 0
          yy[,notPA][w1] <- 1
        }
        
        shan <- sweep(yy[,inRichness], 1, rowSums(yy[,inRichness]), '/')
        shan[shan == 0] <- NA
        shan <- -rowSums(shan*log(shan),na.rm=T)
        shannon <- shannon + shan
        
        wpp <- which(yy > 0)
        ypredPres[wpp]  <- ypredPres[wpp] + yp[wpp]
        ypredPres2[wpp] <- ypredPres2[wpp] + yp[wpp]^2
        ypredPresN[wpp] <- ypredPresN[wpp] + 1
        
        presence[,inRichness] <- presence[,inRichness] + yy[,inRichness]
        ones <- round(rowSums(yy[,inRichness]))
        more <- round(rowSums(yy[,inRichness]*wrich[,inRichness,drop=F]))
        richFull <- .add2matrix(ones,richFull)
        richness <- .add2matrix(more,richness)  # only for non-missing
      }
      
      if(RANDOM){
        alphaRanSums <- alphaRanSums + alphaRandGroup
      }
      
      if(PREDICTX & length(predXcols) > 0){
        predx  <- predx + xpred
        predx2 <- predx2 + xpred^2
      }
      
      wa0 <- which(colSums(agg) != 0)
      ess[notOther[wa0],notOther[wa0]]  <- 
        t(agg[,wa0,drop=F])%*%covE%*%agg[,wa0,drop=F] 
      
      emat[notOther[wa0],notOther[wa0]] <- 
        emat[notOther[wa0],notOther[wa0]] + 
        .cov2Cor( ess[notOther[wa0],notOther[wa0]] )
      
      lo[ ess < 0 ] <- lo[ ess < 0 ] + 1
      hi[ ess > 0 ] <- hi[ ess > 0 ] + 1
      
      ess[notOther,notOther] <- ginv(ess[notOther,notOther])
      
      lm[ ess < 0 ] <- lm[ ess < 0 ] + 1  # neg values
      hm[ ess > 0 ] <- hm[ ess > 0 ] + 1  # pos values
      
      if(REDUCT){
        rndTot <- rndTot + rndEff
      }
      
    }
  }     
  
  ################# end gibbs loop ####################
  
  
  otherpar$S <- S 
  otherpar$Q <- Q
  otherpar$snames <- snames
  otherpar$xnames <- xnames
  
  presence <- presence/ntot
  
  if(RICHNESS){
    missRows <- c()
    richNonMiss <- richness/ntot            #only non-missing plots
    yr  <- as.matrix(ydata[,inRichness]) 
    yr[yr > 0] <- 1
    yr <- rowSums(yr,na.rm=T)
    vv  <- matrix(as.numeric(colnames(richNonMiss)),n,
                  ncol(richNonMiss),byrow=T)
    rmu <- rowSums( vv * richNonMiss )/rowSums(richNonMiss)
    
    rsd <- sqrt( rowSums( vv^2 * richNonMiss )/rowSums(richNonMiss) - rmu^2)
    
    vv  <- matrix(as.numeric(colnames(richFull)),n,ncol(richFull),byrow=T)
    rfull <- rowSums( vv * richFull )/rowSums(richFull)
    rfull[missRows] <- NA
    rmu <- rowSums(presence)
    
    shan <- sweep(y[,inRichness], 1, rowSums(y[,inRichness]), '/')
    shan[shan == 0] <- NA
    shanObs <- -rowSums(shan*log(shan),na.rm=T)
    
    richness <- cbind(yr, rmu, rsd, rfull, shanObs, shannon/ntot )
    colnames(richness) <- c('obs','predMu','predSd','predNotMissing',
                            'H_obs', 'H_pred')
    
    ypredPresMu  <- ypredPres/ypredPresN   #predictive mean and se given presence
    ypredPresMu[ypredPresN == 0] <- 0
    yvv <- ypredPres2/ypredPresN - ypredPresMu^2
    yvv[!is.finite(yvv)] <- 0
    ypredPresSe <- sqrt(yvv)
  }
  
  xunstand    <- .getUnstandX(x, standRows, standMu[,1],
                              standMat[,1], interBeta$intMat)$xu
  
  rmspeBySpec <- sqrt( colSums(yerror)/ntot/n )
  rmspeAll    <- sqrt( sum(yerror)/ntot/n/S )
  
  sMean <- sMean/ntot
  
  tmp <- .chain2tab(bgibbs[burnin:ng,], snames, xnames)
  betaStandXmu <- tmp$mu
  betaStandXTable <- tmp$tab
  
  tmp <- .chain2tab(bgibbsUn[burnin:ng,], snames, xnames)
  betaMu <- tmp$mu
  betaTable <- tmp$tab
  
  tmp <- .chain2tab(bFacGibbs[burnin:ng,], snames[notOther], rownames(agg))
  betaStandXWmu <- tmp$mu
  betaStandXWTable <- tmp$tab
  
  tmp <- .chain2tab(fSensGibbs[burnin:ng,,drop=F])
  sensTable <- tmp$tab[,1:4]
  
  yMu <- ypred/ntot
  y22 <- ypred2/ntot - yMu^2
  y22[y22 < 0] <- 0
  ySd <- sqrt(y22)
  
  cMu <- cuts
  cSe <- numeric(0)
  
  wMu <- wpred/ntot
  wpp <- pmax(0,wpred2/ntot - wMu^2)
  wSd <- sqrt(wpp)
  
  meanDev <- sumDev/ntot
  
  tmp <- .dMVN(wMu[,notOther],x%*%betaMu[,notOther],
               sMean[notOther,notOther], log=T)
  pd  <- meanDev - 2*sum(tmp )
  DIC <- pd + meanDev
  
  yscore <- colSums( .getScoreNorm(y[,notOther],yMu[,notOther],
                                   ySd[,notOther]^2),na.rm=T )  # gaussian w
  xscore <- xpredMu <- xpredSd <- NULL
  standX <- xmissMu <- xmissSe <- NULL
  
  if(RANDOM){
    ns <- 500
    simIndex <- sample(burnin:ng,ns,replace=T)
    tmp <- .expandSigmaChains(snames, alphaVarGibbs, otherpar, simIndex=simIndex,
                              sigErrGibbs, kgibbs, REDUCT=F)
    alphaRandGroupVarMu <- tmp$sMu
    alphaRandGroupVarSe <- tmp$sSe
    alphaRandByGroup <- alphaRanSums/ntot
    
  }
  
  if(PREDICTX){
    xpredMu <- predx/ntot
    xpredSd <- predx2/ntot - xpredMu^2
    xpredSd[xpredSd < 0] <- 0
    xpredSd <- sqrt(xpredSd)
    
    xrow <- standRows
    xmu  <- standMu[,1]
    xsd  <- standMat[,1]
    
    xpredMu <- .getUnstandX(xpredMu, xrow, xmu, xsd, intMat)$xu
    xpredSd[,xrow] <- xpredSd[,xrow]*matrix( xsd[xrow], n, length(xrow),
                                             byrow=T ) 
    
    if(Q == 2)xscore <- mean( .getScoreNorm(x[,2],
                                              xpredMu[,2],xpredSd[,2]^2) )
    if(Q > 2)xscore <- colMeans(.getScoreNorm(x[,-1],
                                                xpredMu[,-1],xpredSd[,-1]^2) )
  }
  
  if(length(standRows) > 0){                #unstandardize
    standX <- cbind(standMu[,1],standMat[,1])
    colnames(standX) <- c('xmean','xsd')
    rownames(standX) <- rownames(standMat)
  }
  
  # betaSens, sigma and R
  
  ns <- 500
  simIndex <- sample(burnin:ng,ns,replace=T)
  
  tmp <- .expandSigmaChains(snames, sgibbs, otherpar, simIndex=simIndex,
                            sigErrGibbs, kgibbs, REDUCT)
  corMu <- tmp$rMu; corSe <- tmp$rSe
  sigMu  <- tmp$sMu; sigSe  <- tmp$sSe
  
  whichZero <- which(lo/ntot < ematAlpha & 
                       hi/ntot < ematAlpha,arr.ind=T) #not different from zero
  whConZero <- which(lm/ntot < ematAlpha & 
                       hm/ntot < ematAlpha,arr.ind=T)
  
  ematrix  <- emat/ntot
  fmatrix  <- fmat/ntot
  
  tMu <- tSd <- tMuOrd <- btMu <- btSe <- stMu <- stSe <- numeric(0)
  
  if('PA' %in% typeNames){
    zMu <- yMu
    zSd <- ySd
  }
  
  
  # outputs
  if(length(reductList) == 0)reductList <- list(N = 0, r = 0)
  reductList$otherpar <- otherpar
  
  modelList$effort    <- effort;      modelList$formula <- formula
  modelList$typeNames <- typeNames;    modelList$censor <- censor
  modelList$effort    <- effort; modelList$holdoutIndex <- holdoutIndex
  modelList$REDUCT    <- REDUCT;       modelList$TRAITS <- TRAITS
  modelList$ematAlpha <- ematAlpha; modelList$traitList <- traitList
  modelList$reductList <- reductList; modelList$ng <- ng
  modelList$burnin <- burnin
  
  inputs <- list(xdata = xdata, x = xunstand, standX = standX,
                 standMat = standMat, standRows = standRows, y = y, 
                 notOther = notOther, other = other, breakMat = breakMat, 
                 designTable = designTable, classBySpec = classBySpec, 
                 factorBeta = factorBeta, interBeta = interBeta,
                 linFactor = linFactor, intMat = intMat, RANDOM = RANDOM)
  missing <- list(xmiss = xmiss, xmissMu = xmissMu, xmissSe = xmissSe, 
                  ymiss = ymiss, ymissMu = ymissPred, ymissSe = ymissPred2)
  parameters <- list(betaMu = betaMu, betaTable = betaTable, 
                     betaStandXmu = betaStandXmu, 
                     betaStandXTable = betaStandXTable,
                     betaStandXWmu =  betaStandXWmu,
                     betaStandXWTable = betaStandXWTable,
                     corMu = corMu, corSe = corSe, 
                     sigMu = sigMu, sigSe = sigSe, 
                     ematrix = ematrix, fmatrix = fmatrix,
                     whichZero = whichZero, whConZero = whConZero,
                     wMu = wMu, wSd = wSd, sensTable = sensTable)
  prediction <- list(presence = presence, xpredMu = xpredMu, xpredSd = xpredSd,
                     ypredMu = yMu, ypredSd = ySd, richness = richness)
  chains <- list(sgibbs = sgibbs, bgibbs = bgibbs, bgibbsUn = bgibbsUn,
                 fSensGibbs = fSensGibbs, bFacGibbs = bFacGibbs) 
  fit <- list(DIC = DIC, yscore = yscore, 
              xscore = xscore, rmspeAll = rmspeAll,
              rmspeBySpec = rmspeBySpec)
  if(FULL)chains <- append(chains, list(ygibbs = ygibbs))
  if(RANDOM){
    parameters <- append(parameters,
                         list( randGroupVarMu = alphaRandGroupVarMu,
                               randGroupVarSe = alphaRandGroupVarSe,
                               randByGroup = alphaRandByGroup) )
  }
  if(RICHNESS){
    prediction <- append(prediction, 
                         list(yPresentMu = ypredPresMu, yPresentSe = ypredPresSe))
  }
  if(REDUCT) {
    parameters <- append(parameters, list(rndEff = rndTot/ntot))#, specRand = specRand))
    chains <- append(chains,list(kgibbs = kgibbs, sigErrGibbs = sigErrGibbs))
  }
  
  chains     <- chains[ sort( names(chains) )]
  fit        <- fit[ sort( names(fit) )]
  inputs     <- inputs[ sort( names(inputs) )]
  missing    <- missing[ sort( names(missing) )]
  modelList  <- modelList[ sort( names(modelList) )]
  parameters <- parameters[ sort( names(parameters) )]
  prediction <- prediction[ sort( names(prediction) )]
  
  all <- list(chains = chains, fit = fit, inputs = inputs, missing = missing,
              modelList = modelList, parameters = parameters,
              prediction = prediction)
  all$call <- match.call()
  all <- all[ sort(names(all)) ]
  class(all) <- "gjam"
  
  all
}
