rm ( list = ls() )

library("Rcpp")
library("mcclust")
library("ggplot2")
library("raster")
library("leaflet")
library("gtable")
library("grid") 

#### INPUT OF THE DATA ####

data.input = read.csv ( "dataset/dataset.csv" )
ncells = nrow(data.input)
data.input = data.input[1:ncells,1:327]

# Here we compute the coordinates for each of the cells
latRange = c ( 45.36587, 45.57 )
lat = (seq(0,nrow(data.input)-1) %% 48) * (latRange[2] - latRange[1]) / 47 + latRange[1]

lonRange = c ( 9.05, 9.34444 )
lon = floor(seq(0,nrow(data.input)-1) / 48) * (lonRange[2] - lonRange[1]) / 53 + lonRange[1]

# Some cells (actually, only the cell 1569) have only NAs; we ignore them because
# we have absolutely no information about them
count.na = function ( x ) { sum(is.na(x)) }
#na.rows = which ( apply ( data.input, 1, count.na ) == ncol ( data.input ) )
na.rows = c()
if ( length(na.rows) > 0 ) data.input = data.input[-na.rows,]

# Transform the data so that they are supported on the real line
y = as.matrix ( log ( data.input + 1e-8 ) )

#### CONSTRUCTION OF THE COVARIATES ####

harmonics = c(1:14) # Number of harmonics we want to consider
times = seq ( 1, ncol(data.input) ) / ncol(data.input) # Vector of the times

x = matrix ( nrow = length(times), ncol = length(harmonics) * 2 + 1 )
x[,1] = rep ( 1, times = length(times) )

if ( length(harmonics) > 0 )
  for ( h in 1:length(harmonics) ) {
     x[,2*h] = cos ( 2*pi*harmonics[h]*times )
     x[,2*h+1] = sin ( 2*pi*harmonics[h]*times )
  }

# Weekday-weekend flag; 0 = weekday, 1 = weekend
phi = c ( rep(0, times=71), rep(1, times=48), rep(0, times=120), rep(1, times=48), rep(0, times=40) )

# Overall covariate matrix
h = matrix ( nrow = length(times), ncol = 2 * ncol(x) )
for ( t in 1:length(times) ) {
  if ( phi[t] == 0 ) h[t,] = c ( x[t,], rep(0, times = ncol(x) ) )
  else h[t,] = c ( rep(0, times = ncol(x)), x[t,] )
}

#### FIT OF THE MODEL ####

posteriorDraws = 1000
burnInIterations = 100
trim = 1

sourceCpp ( "fitAnovaDDP.cpp", cacheDir = "." )
set.seed ( 15 )
fitOutput = fitAnovaDDP ( y, h, posteriorDraws, burnInIterations, trim, a0 = 0, b0 = 0, kappa0 = 0, nu0 = 0 )

#### GRAPHICAL GOODNESS OF FIT CHECK ####
convergencePlot = function ( x, title = "" ) {
  layout ( matrix(c(1,1,2,3),nrow=2,ncol=2,byrow=T) )
  hist ( x, main = title, breaks = 20 )
  plot ( x, type = "l", main = "Traceplot", xlab="", ylab="" )
  acf ( x, main = "Autocorrelation" )
}

pdf ( "report/figures/tauBeta.pdf", width=7, height=6 )
convergencePlot ( 1 / fitOutput$tauBeta[,1], "sigmaBeta" )
dev.off()
pdf ( "report/figures/sigma.pdf", width=7, height=6 )
convergencePlot ( fitOutput$sigma, "sigma" )
dev.off()

cell = 500
ypred = rep ( 0, times = ncol(y) )
yhigh = rep ( 0, times = ncol(y) )
ylow = rep ( 0, times = ncol(y) )
betaMean = matrix ( rep(0, times=ncol(h)), ncol = 1, nrow = ncol(h) )

for ( b in 1:length(fitOutput$beta) ) {
  betaMean = betaMean + fitOutput$beta[[b]][fitOutput$labels[[b]][cell] + 1, ]
}

betaMean = betaMean / length(fitOutput$beta)
ypred = h %*% betaMean

for ( t in 1:ncol(y) ) {
  residual.high = function ( c ) {
    sum = 0.0
    for ( b in 1:length(fitOutput$beta) ) {
      sum = sum + pnorm ( c, mean = h[t,] %*% fitOutput$beta[[b]][fitOutput$labels[[b]][cell]+ 1, ], sd = sqrt(fitOutput$sigma[b]) )
    }
    sum = sum / length(fitOutput$beta) - 0.95
    sum
  }

  residual.low = function ( c ) {
    sum = 0.0
    for ( b in 1:length(fitOutput$beta) ) {
      sum = sum + pnorm ( c, mean = h[t,] %*% fitOutput$beta[[b]][fitOutput$labels[[b]][cell]+ 1, ], sd = sqrt(fitOutput$sigma[b]) )
    }
    sum = sum / length(fitOutput$beta) - 0.05
    sum
  }

  yhigh[t] = uniroot ( residual.high, interval = c(ypred[t] - 2, ypred[t] + 2), extendInt = "yes" )$root
  ylow[t] = uniroot ( residual.low, interval = c(ypred[t] - 2, ypred[t] + 2), extendInt = "yes" )$root
}

data.toplot = data.frame ( t = times * 327, y = y[cell,], what = rep ( 0, times = length(times)) )
data.toplot = rbind ( data.toplot, data.frame ( t = times * 327, y = ypred, what = rep ( 1, times = length(times)) ) )
data.toplot = rbind ( data.toplot, data.frame ( t = times * 327, y = ylow, what = rep ( 2, times = length(times)) ) )
data.toplot = rbind ( data.toplot, data.frame ( t = times * 327, y = yhigh, what = rep ( 3, times = length(times)) ) )

pdf ( paste('./report/figures/predictiveCell', cell, '.pdf', sep = ''), width = 7, height = 3.5 )
plot1 = ggplot ( data = data.toplot, aes (x = t, y = y, group = what, color = as.factor(what)) ) +
 geom_line() +
 geom_vline ( xintercept = seq(0.25,327.25,24), linetype = 3 ) + 
 theme_bw() + theme ( panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) +
 scale_x_continuous ( name = "t [h]", limits=c(0,327), breaks=seq(12.25,327.25,24), minor_breaks=NULL,
                      labels = rep(c("Wed", "Thu", "Fri", "Sat", "Sun", "Mon", "Tue"),times=2) ) +
 scale_y_continuous ( name = "log(Erlang)", limits=c(-15,10) ) +
 scale_color_manual ( values = c("black", "red", "blue", "purple" ),
                      labels = c("observed", "pred. mean", "pred. 0.05", "pred. 0.95"), name = paste("Cell n.", cell, "") )
plot1
dev.off()

locationMap = leaflet() %>%
  fitBounds ( lat1=min(lat), lat2=max(lat), lng1=min(lon), lng2=max(lon) ) %>%
  addProviderTiles(providers$Stamen.TonerLines) %>%
  addRectangles(lat1 = latRange[1], lat2 = latRange[2], lng1 = lonRange[1], lng2 = lonRange[2], color="black", fill = F, weight = 4, opacity = 1 ) %>%
addCircleMarkers(lat = lat[10], lng = lon[10], stroke = F, color = "red", fillOpacity = 1, radius = 5 )
locationMap

#### CLUSTERING ####
 
clustering = minbinder ( fitOutput$clusterMatrixPost )
clustering.labels = clustering$cl
for ( i in na.rows )
  clustering.labels = c(clustering.labels[1:i-1], NA, clustering.labels[i:length(clustering.labels)] )

cellplot.data = matrix ( clustering.labels, nrow = 48, ncol = 54, byrow=F );
cellplot = raster ( nrows=48, ncols=54, xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), vals=cellplot.data )
pal = colorNumeric(c("#2962FF", "#FFD600", "#DD2C00", "#00C853", "#AA00FF"), values(cellplot), na.color = "white")
map = leaflet() %>%
  fitBounds ( lat1=min(lat), lat2=max(lat), lng1=min(lon), lng2=max(lon) ) %>%
  addRasterImage ( cellplot, pal, opacity=1, method="ngb" ) %>%
  addProviderTiles(providers$Stamen.TonerLines)
map

data.toplot = data.frame ()
for ( i in 1:nrow(y) ) {
  data.toplot = rbind ( data.toplot,
                        data.frame ( t = times*327, y = as.numeric(t(data.input[i,])), cellIdx = i, cluster = clustering.labels[i] ) )
}
pdf ( "report/figures/erlangClusterColored.pdf", width=7, height=4 )
plot2 = ggplot ( data = data.toplot, aes (x = t, y = y, group = cellIdx, color = as.factor(cluster)) ) +
  geom_line() +
  geom_vline ( xintercept = seq(0.25,327.25,24), linetype = 3 ) + 
  theme_bw() + theme ( panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom" ) +
  scale_x_continuous ( name = "t [h]", limits=c(0,327), breaks=seq(12.25,327.25,24), minor_breaks=NULL,
                       labels = rep(c("Wed", "Thu", "Fri", "Sat", "Sun", "Mon", "Tue"),times=2) ) +
  scale_y_continuous ( name = "Erlang" ) + 
  scale_color_manual ( values = c("#2962FF", "#FFD600", "#DD2C00", "#00C853", "#AA00FF"), name = "cluster" )
plot2
dev.off()
