for(g in 1:ng){ ########################################################
  
  if(REDUCT){
    
    #   if(g > burnin)CLUST <- F
    
    Y <- w[,notOther]
    if(RANDOM)Y <- Y - groupRandEff[,notOther] 
    if(TIME)  Y <- Y - mua[,notOther] - mug[,notOther] 
    
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
    
    if(!TIME){
      Y <- w[inSamp,notOther] - rndEff[inSamp,notOther]
      if(RANDOM)Y <- Y - groupRandEff[inSamp,notOther]
      bg[,notOther] <- updateBeta(X = x[inSamp,], Y, 
                                  sig = sigmaerror, beta = bg[,notOther],
                                  lo=loB[,notOther], hi=hiB[,notOther])
      muw[inSamp,] <- x[inSamp,]%*%bg
      
    } else {
      
      mua  <- Umat%*%Amat 
      mug  <- Vmat%*%Lmat
      
      Y <- w[,notOther] - mua[,notOther] - mug[,notOther] - rndEff[,notOther]
      if(RANDOM)Y <- Y - groupRandEff[,notOther]
      bg[,notOther] <- updateBeta(X = x[tindex[,2],], Y = Y[tindex[,2],], 
                                  sig = sigmaerror, beta = bg[,notOther],
                                  rows = Brows, pattern = Bpattern,
                                  lo=loB[,notOther], hi=hiB[,notOther])
      mub <- x%*%bg
      Y   <- w - mub - mua - rndEff
      if(RANDOM)Y <- Y - groupRandEff
      
      Lmat[,notOther] <- updateBeta(X = Vmat[tindex[,2],], 
                                    Y = Y[tindex[,2],notOther], sig=sigmaerror, 
                                    beta = Lmat[,notOther],
                                    rows = Lrows, pattern = Lpattern, 
                                    lo=loLmat, hi=hiLmat, ixx=F)
      
      #     Lmat[,notOther] <- .updateBetaMet(X = Vmat[tindex[,2],], 
      #                                       Y[tindex[,2],notOther], 
      #                                       B = Lmat[,notOther],
      #                           lo=loLmat, hi=hiLmat, loc = wL, REDUCT, 
      #                           sig=sigmaerror,sp=spL)
      mug  <- Vmat%*%Lmat
      Y    <- w - mub - mug - rndEff
      if(RANDOM)Y <- Y - groupRandEff
      Amat <- updateBeta(X = Umat[tindex[,2],], Y[tindex[,2],], sig=sigmaerror, 
                         rows = Arows, pattern = Apattern, 
                         beta = Amat,
                         lo=loAmat, hi=hiAmat, ixx=F)
      #     Amat <- .updateBetaMet(X = Umat[tindex[,2],], Y[tindex[,2],notOther], 
      #                                       B = Amat,
      #                                       lo=loAmat, hi=hiAmat, loc = wA, REDUCT, 
      #                                      sig=sigmaerror,sp=rexp(nA,1/spA))
      mua <- Umat%*%Amat
      
      #     if(g %in% gcheck){
      #       g2   <- g - 1
      #       spA <- apply(alphaGibbs[g1:g2,],2,sd)/2 + tinyg
      #       spL <- apply(ggibbs[g1:g2,],2,sd)/2 + tinyg
      #       if(g < 200)g1 <- g
      #     }
      
      muw <- mub + mug + mua + rndEff
    }
    
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
  
  if( 'OC' %in% typeCode ){
    tg <- .updateTheta(w,tg,cutLo,cutHi,ordCols,
                       holdoutN,holdoutIndex,minOrd,maxOrd) # var scale
    cutg <- .gjamCuts2theta(tg,ss = sg[ordCols,ordCols]) # corr scale
    breakMat[ordCols,1:lastOrd] <- cutg
    cgibbs[g,] <- as.vector( cutg[,-c(1,2,ncut)] )
    
    plo[,ordCols] <- cutg[cutLo]
    phi[,ordCols] <- cutg[cutHi]
  }
  
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
  
  if(TIME){
    
    #    muw does not include groupRandEff
    tmp <- .updateW(w,plo,phi,wpropTime,xl,yp,Lmat,Amat,mub,rndEff, groupRandEff,
                    sdg,muw,Umat,Vmat,sinv)
    w <- tmp$w; muw <- tmp$muw; yp <- tmp$yp; Umat <- tmp$Umat; Vmat <- tmp$Vmat
    
    groups <- NULL
    
    for(k in allTypes){
      
      wk <- which(typeCols == k)
      nk <- length(wk)
      wo <- which(wk %in% notOther)
      wu <- which(typeCols[notOther] == k)
      wp <- w[, wk, drop=F]
      yp <- yp[, wk, drop=F]
      
      if(typeFull[wk[1]] == 'countComp')groups <- CCgroups
      if(typeFull[wk[1]] == 'fracComp')groups  <- FCgroups
      if(typeFull[wk[1]] == 'categorical')groups <- CATgroups
      
      glist <- list(wo = wo, type = typeFull[wk[1]], yy = y[,wk,drop=F], 
                    wq = wp, yq = yp, cutg = cutg, 
                    censor = censor, censorCA = censorCA, 
                    censorDA = censorDA, censorCON  = censorCON, 
                    eff = effMat[,wk,drop=F], groups = groups, 
                    k = k, typeCols = typeCols, notOther = notOther, 
                    wk = wk, sampW = sampleW[,wk])
      tmp <- .gjamWLoopTypes( glist )
      w[,wk]  <- tmp[[1]]
      yp[,wk] <- tmp[[2]]
    }
    
    #predict X
    
    ww <- w
    ww[ww < 0] <- 0
    mua  <- Umat%*%Amat
    mug  <- Vmat%*%Lmat
    muw <- mua + mub + mug + rndEff
    
    xtmp <- xpred
    xtmp[,-1] <- .tnorm(n*Qall,-3,3,xpred[,-1],.1)
    
    # factors
    if( length(linFactor) > 0 ){
      
      for(k in 1:length(linFactor)){
        
        mm  <- linFactor[[k]]
        wcol <- sample(mm,n,replace=T)
        xtmp[,mm[-1]] <- 0
        xtmp[ cbind(1:n, wcol) ] <- 1
        
      }
    }
    
    if(length(intMat) > 0){     #  interactions
      xtmp[,intMat[,1]] <- xtmp[,intMat[,2]]*xtmp[,intMat[,3]]
    }
    
    ae     <- mua + rndEff
    Vnow   <- Vmat
    mubNow <- xpred[,xnames]%*%bg
    mubNew <- xtmp[,xnames]%*%bg
    
    Vnow[tindex[,2],] <- ww[tindex[,1],gindex[,'colW']]*
      xpred[tindex[,2],xlnames][,gindex[,'rowG']]
    Vnow[timeZero+1,] <- ww[timeZero,gindex[,'colW']]*
      xpred[timeZero+1,xlnames][,gindex[,'rowG']]
    mugNow <- Vnow%*%Lmat
    muNow  <- mubNow + mugNow + ae
    
    Vnew[tindex[,2],] <- ww[tindex[,1],gindex[,'colW']]*
      xtmp[tindex[,2],xlnames][,gindex[,'rowG']]
    Vnew[timeZero+1,] <- ww[timeZero,gindex[,'colW']]*
      xtmp[timeZero+1,xlnames][,gindex[,'rowG']]
    mugNew <- Vnew%*%Lmat
    muNew  <- mubNew + mugNew + ae
    
    if(REDUCT){
      pnow <- dnorm(w[,notOther],muNow[,notOther],sdg,log=T)
      pnew <- dnorm(w[,notOther],muNew[,notOther],sdg,log=T)
      a1   <- exp( rowSums(pnew - pnow) )
    }else{
      pnow <- .dMVN(w[tindex[,2],notOther],muNow,sinv=sinv,log=T) 
      pnew <- .dMVN(w[tindex[,2],notOther],muNew,sinv=sinv,log=T) 
      a1   <- exp(pnew - pnow)
    }
    z    <- runif(length(a1),0,1)
    za   <- which(z < a1)
    if(length(za) > 0){
      xpred[za,] <- xtmp[za,]
      Vmat[za,] <- Vnew[za,]
      muw[za,]  <- muNew[za,]
      mub[za,]  <- mubNew[za,]
      mug[za,]  <- mugNew[za,]
    }
    
    if(nlmiss > 0)xl[xlmiss] <- xpred[xmiss]
    
    if(nmiss > 0){
      
      x[xmiss] <- xpred[xmiss]
      
      tmp    <- .getUnstandX(x, standRows, standMu[,1],
                             standMat[,1], intMat)            
      S2U    <- tmp$S2U
      XX     <- crossprod(x)
      IXX    <- solveRcpp(XX)
    }
    
    ggibbs[g,]     <- Lmat[wL]
    alphaGibbs[g,] <- Amat[wA]
    
  } else{ #############not TIME
    
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
    
    if(nmiss > 0){
      
      x[xmiss] <- .imputX_MVN(x,Y,bg[,notOther],xmiss,sinv,xprior=xprior,
                              xbound=xbound)[xmiss]
      tmp      <- .getUnstandX(x, standRows, standMu[,1],
                               standMat[,1], intMat)            
      S2U    <- tmp$S2U
      XX     <- crossprod(x)
      IXX    <- solveRcpp(XX)
    }
    
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
  }
  
  setTxtProgressBar(pbar,g)
  
  bgu <- bg                    # unstandardize beta
  if(length(standRows) > 0){
    if(TIME){
      bgu <- S2U%*%mub
      lambda[ gindex[,c('rowG','colW')]] <- Lmat[wL]
      lambdas <- S2UL%*%mug      # unstandardized lambda
      lgibbs[g,] <- lambdas[,notOther]
    }else{
      bgu <- S2U%*%x%*%bg
    }
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
  
  if(TRAITS){
    Atrait <- bg%*%t(specTrait[,colnames(yp)])  # standardized
    Strait <- specTrait[,colnames(yp)]%*%sg%*%t(specTrait[,colnames(yp)])
    bTraitGibbs[g,] <- Atrait
    mgibbs[g,] <- Strait
    
    minv <- ginv(Strait)
    
    tmp <- .contrastCoeff(beta=Atrait, 
                          notStand = notStandard[notStandard %in% xnames], 
                          sigma = Strait, sinv = minv,
                          stand = standMat, factorObject=factorBeta )
    tagg   <- tmp$ag
    bTraitFacGibbs[g,] <- tagg # stand for X and W, centered for factors
  }
  
  if(TIME){
    
    tmp <- .contrastCoeff(beta=lambda[,notOther], 
                          notStand = notStandardL[notStandardL %in% xlnames], 
                          sigma = sg[notOther,notOther],sinv = sinv,
                          stand=standMatL, factorObject=factorLambda)
    lgg   <- tmp$ag
    leg   <- tmp$eg
    lsens <- tmp$sens
    
    lss <- sqrt(diag(lsens))
    
    if(g == 1){
      if( !all(names(lss) %in% colnames(gsensGibbs)) )
        colnames(gsensGibbs) <- names(lss)
    }
    
    gsensGibbs[g,names(lss)] <- lss
    
    alpha[ aindex[,c('toW','fromW')] ] <- Amat[wA]
    asens <- Amat[,notOther]%*%sinv%*%t(Amat[,notOther])
    asens <- sqrt(diag(asens))
    asensGibbs[g,] <- asens
  }
  
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
    
    if(mmiss > 0){
      ymissPred[ymiss]  <- ymissPred[ymiss] + y[ymiss]
      ymissPred2[ymiss] <- ymissPred2[ymiss] + y[ymiss]^2
    }
    if(nmiss > 0){
      xmissSum  <- xmissSum + x[xmiss]
      xmissSum2 <- xmissSum2 + x[xmiss]^2
    }
    
    if(PREDICTX & length(predXcols) > 0){
      predx  <- predx + xpred
      predx2 <- predx2 + xpred^2
    }
    
    wa0 <- which(colSums(agg) != 0)
    ess[notOther[wa0],notOther[wa0]]  <- 
      t(agg[,wa0,drop=F])%*%covE%*%agg[,wa0,drop=F] 
    if(TIME){
      wa0 <- which(colSums(lgg) != 0)
      ess[notOther[wa0],notOther[wa0]]  <- 
        ess[notOther[wa0],notOther[wa0]] +
        t(lgg[,wa0,drop=F])%*%covL%*%lgg[,wa0,drop=F] 
    }
    
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
    
    if(TRAITS){
      yw     <- sweep(yp,1,rowSums(yp),'/')
      yw[yw <= 0]   <- 0
      yw[is.na(yw)] <- 0
      Ttrait <- .gjamPredictTraits(yw,specTrait[,colnames(yp)], traitTypes)
      tpred  <- tpred + Ttrait
      tpred2 <- tpred2 + Ttrait^2
    }
  }
}     

################# end gibbs loop ####################