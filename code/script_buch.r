rm ( list = ls() )
library("mvtnorm")

#### INPUT OF THE DATA ####

data.input = read.csv ( "dataset_bubi.csv" )

# Transform the data so that they are supported on the real line
y = as.matrix ( log ( data.input + 1e-8 ) )
NAs = which ( is.na(y) )

# Here we compute the coordinates for each of the cells
lonRange = c ( 9.05, 9.34444 )
latRange = c ( 45.36587, 45.57 )
lat = (seq(0,nrow(data.input)-1) %% 48) * (latRange[2] - latRange[1]) / 47 + latRange[1]
lon = floor(seq(0,nrow(data.input)-1) / 48) * (lonRange[2] - lonRange[1]) / 53 + lonRange[1]

#### CONSTRUCTION OF THE COVARIATES ####

harmonics = 5 # Number of harmonics we want to consider
times = seq ( 1, ncol(data.input) ) / ncol(data.input) # Vector of the times

x = matrix ( nrow = length(times), ncol = harmonics * 2 + 1 )
x[,1] = rep ( 1, times = length(times) )

for ( h in 1:harmonics ) {
  x[,2*h] = cos ( 2*pi*h*times )
  x[,2*h+1] = sin ( 2*pi*h*times )
}

# Weekday-weekend flag; 0 = weekday, 1 = weekend
phi = c ( rep(0, times=71), rep(1, times=47), rep(0, times=119), rep(1, times=47), rep(0, times=51) )

# Overall covariate matrix
h = matrix ( nrow = length(times), ncol = 2 * ncol(x) )
for ( t in 1:length(times) ) {
  if ( phi[t] == 0 ) h[t,] = c ( x[t,], rep(0, times = ncol(x) ) )
  else h[t,] = c ( x[t,], x[t,] )
}

#### HYPERPARAMETERS ####

sigma1 = 1 # Variance of yi
muBeta = rep ( 0, times=2*ncol(x) ) # Prior mean of the beta coefficients
sigmaBeta = diag ( 2*ncol(x) ) # Prior variance-covariance matrix of the beta coefficients
alpha0 = 1 # Mass parameter of the Dirichlet process

#### SOME MATRICES TO BE PRECOMPUTED ####

sigma2 = sigma1*sigma1
sigmaBetaInv = solve(sigmaBeta)
detSigmaBeta = det(sigmaBeta)
sigmaInv = t(h) %*% h / sigma2 + sigmaBetaInv
sigma = solve ( sigmaInv )
detSigma = det(sigma)
muBetaSigmaBetaMuBeta = t(muBeta) %*% sigmaBetaInv %*% muBeta / 2

#### FIT OF THE MODEL ####

cells = nrow ( y )

# State of the Markov chain
#ALL'INIZIO HO UJ SOLO LABELS, OSSIA 1
labels = rep ( 1, times = cells )
beta = matrix ( rep(0, times=ncol(h)), nrow = 1, ncol = ncol(h) )
#ALL'INIZIO clusterSize è UN VETTORE DI UN ELEMENTO INIZIALLIZZATO A 2592
clusterSize = c ( cells )
clusterMatrix = matrix ( rep(1, times=cells*cells), nrow = cells, ncol = cells )
y[NAs] = rep ( 0, times = length(NAs) )

# Number of iterations
posteriorDraws = 1
burnInIterations = 0

for ( niter in 1:(posteriorDraws + burnInIterations) ) { # MCMC loop
  # Step 1 : resample the cluster assignments
  #for each cell i = 1,...,n do
  for ( i in 1:cells ) { # Loop on cells
    cat ( "\014" )
    cat ( "iter: ", niter, " cell: ", i, sep = "" )
    #EQUATION 12 IN BUCH REPORT
    betaBar = solve ( sigmaInv, t( y[i,] %*% h / sigma2 + t(muBeta) %*% sigmaBetaInv ) )
    
    log.prob.new = log(alpha0) - length(times) / 2 * log ( 2*pi*sigma2 ) - log ( detSigmaBeta ) / 2 + log ( detSigma ) / 2
                 - t(y[i,]) %*% y[i,] / (2*sigma2) - muBetaSigmaBetaMuBeta + t(betaBar) %*% sigmaInv %*% betaBar / 2
    #EQUATION 10 IN BUCH REPORT
    prob.new = exp ( log.prob.new )
    #HERE I TAKE FROM LABELS EX: (1, 2, 2, 1, 3, 4) ONLY (1, 2, 3, 4)
    uniqueLabels = unique(labels)
    #ORA MI CALCOLO LA PROB DI ASSEGANRE LA MIA NUOVA OSSERVAZIONE AD UN CLUSTER GIà ESISTENTE,
    #PER FARLO MI CREO UN VETTORE VUOTO DELLA LUNGHEZZA DEL NUMERO DI CLUSTER CHE HO GIà
    #E MI CALCOLO IN UN CICLO LA PROB DI FINIRE NEL SINGOLO CLUSTER (FACCIO STO RAGIONAMENTO PER OGNI CLUSTER)
    prob.old = rep ( 0, times = length(uniqueLabels) ) # Probability of assigning to an existing cluster
    twopisigma = (2 * pi * sigma2) ^ (-length(times) / 2)
    
    #ALL'INIZIO AVRò SOLO UN'ITERAZIONE POICHé labels è STATO INIZIALLIZZATO A TUTTI 1
    #E PERTANTO HO UN SOLO uniquelabels
    for ( c in 1:length(uniqueLabels) ) {
      if ( uniqueLabels[c] != labels[i] ) {
        #EQUATION 9 IN BUCH REPORT (LAST TERM IN THE EXPONENTIAL)
        dx = y[i,] - h %*% beta[c,]
        #EQUATION 9 IN BUCH REPORT -> WE UPDATE THE PROB OF A NEW OBSERVATION TO BELONG TO AN
        #                             EXISTING CLUSTER
        #ALL'INIZIO AVRò SOLO UNA PROB OLD PERCHé HO UN SOLO CLUSTER E AVRò SOLO UN ELEMENTO 
        #IN clusterSize, PERTANTO clusterSize[c] HA SENSO PERCHé ACCEDO ALL'UNICO ELEMENTO
        prob.old[c] = clusterSize[c] * twopisigma * exp ( -t(dx) %*% dx / (2*sigma2) )
      }
    }
    
    # MI CREO nci,−i
    # The i-th cell is removed from its cluster
    clusterMatrix[,i] = rep ( 0, times = cells )
    clusterMatrix[i,] = rep ( 0, times = cells )
    clusterSize[labels[i]] = clusterSize[labels[i]] - 1
    
    #IF nci,−i=0 RESAMPLE THE COEFFICIENT beta_c
    if ( clusterSize[labels[i]] == 0 )
      beta[labels[i],] = rep ( NA, times = ncol(h) )
    
    # Then, pick the new cluster for the i-th cell
    a = runif ( 1, 0, sum(prob.old) + prob.new )
    if ( a < prob.new ) {
      labels[i] = max(labels) + 1
      clusterMatrix[i,i] = 1
      clusterSize[labels[i]] = 1
    }
    else {
      a = a - prob.new
      for ( j in (1:cells)[-i] ) {
        if ( a < prob.old[j] ) {
          labels[i] = labels[j];
          clusterMatrix[j,i] = 1;
          clusterMatrix[i,] = clusterMatrix[j,]
          clusterMatrix[,i] = clusterMatrix[i,]
          clusterSize[labels[i]] = clusterSize[labels[i]] + 1
          break
        }
        else a = a - prob.old[j]
      }
    }
    
    # If the cluster of i contains only i, we need to resample its parameters
    if ( clusterSize[labels[i]] == 1 ) {
      while ( nrow(beta) < labels[i] )
        beta = rbind ( beta, rep ( NA, times = ncol(h) ) )
      
      beta[labels[i],] = rmvnorm ( 1, mean = betaBar, sigma = sigma )
    }
  } # End loop on cells
  
  # Step 2 : resample the cluster parameters
  #for each unique cluster label c do
  for ( l in unique(labels) ) {
    sumy = rep ( 0, times = length(times) )
    for ( i in which(labels == l) )
      sumy = sumy + y[i,]
    #STESSA FORMULA DI PRIMA SOLO CHE AL POSTO DI y[i,] %*% h AVREMO sumy %*% h PERCHé 
    #STO PRENDENDO TUTTE LE RIGHE CHE APPARTENGONO AD UN DATO CLUSTER E LE STO SOMMANDO TRA LORO
    betaBar = solve ( sigmaInv, t( sumy %*% h / sigma2 + t(muBeta) %*% sigmaBetaInv ) )
    beta[labels[i],] = rmvnorm ( 1, mean = betaBar, sigma = sigma )
  } 
} # End MCMC loop