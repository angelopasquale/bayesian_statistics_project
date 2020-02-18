#X = # design matrix
#  Y = # plots x species matrix
  #usa colnames(x)
  source(file = "gjamHfunctions.R")
source(file="gjam.R")
  colnames(x)<-c("x_0","x_1")
x<-as.data.frame(x)
attach(x)
  formula <- as.formula( ~x_0+x_1) # dopo la tilde devi mettere le tue variabili che vuoi (come fai sempre con lm tipo)

rl <- list(r =r , N = N_stick) #r sono latent factors, N il livello di troncamento

ml   <- list(ng = 2000, burnin = 50, typeNames = 'PA', reductList = rl) #ng iterazioni, typeNames il tipo di dati


fit<-gjamReduct(formula, xdata = x, ydata = Y, modelList = ml)
