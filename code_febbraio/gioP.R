X = # design matrix
  Y = # plots x species matrix
  #usa colnames(x)
  formula <- as.formula( ~ x[,1]+x[,2]) # dopo la tilde devi mettere le tue variabili che vuoi (come fai sempre con lm tipo)

rl <- list(r =r , N = N_stick) #r sono latent factors, N il livello di troncamento

ml   <- list(ng = 2000, burnin = 50, typeNames = 'PA', reductList = rl) #ng iterazioni, typeNames il tipo di dati


fit<-gjamReduct(formula, xdata = X, ydata = Y, modelList = ml)
