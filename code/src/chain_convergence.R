chain_convergence<-function(chain){
  # # #GRAFICI DELLA CATENA
  # #Analysis of output
   chain<-result$chain
   chain<-as.vector(chain)
   chain<-as.numeric(chain)
   chain<-as.data.frame(chain)
   chain<-as.mcmc(chain)
   chain<-ggs(chain)
  # # #Traceplot
  
  x11()
   ggs_traceplot(chain) + xlab('Iterations') + ylab('A[1,2] ') + theme(text = element_text(size=15, family="LM Roman 10")) + ggtitle('Traceplot')
  # #Running Mean
   x11()
   ggs_running(chain)
  #
  #
  # #Autocorrelation
  x11()
  ggs_autocorrelation(chain)
  #
  x11()
  ggs_histogram(chain)

  effectiveSize(chain$value)
}