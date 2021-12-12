
# function to define the distribution of points
# arguments:
# n: integer, number of knots in the interval
# lenHB: integer, length of the signal
# rpeak: integer, index of the peak R
# punti_medi: vettore di 10  interi, 10 indici dei punti P Q R S T
# plot: boolean, if true plot the signal with 
# frazioni: vettore di 3 float, frazione entro cui dividere i knots nella fascia R, P ,T
#
# returns il vettore dei knots entro cui è stato diviso l'intervallo (di lunghezza a volte maggiore di 1 o 2 di n)
find_knots <- function(n,lenHB,rpeak,punti_medi, plot = FALSE, frazioni = c('R' = 2/5, 'P' = 1/6, 'T' = 1/6)) {
  
  # knots in curva R
  nodiR = seq(from = rpeak + punti_medi[5] - 1,to = rpeak + punti_medi[7] + 2,length.out = n*frazioni[1])
  lenR = length(nodiR)
  
  # knots in curva P
  nodiP = seq(from = rpeak + punti_medi[2] - 3,to = rpeak + punti_medi[3] + 3,length.out = n*frazioni[2])
  lenP = length(nodiP)
  
  # knots in curva T
  nodiT = seq(from = rpeak + punti_medi[9] - 3,to = rpeak + punti_medi[10] + 3,length.out = n*frazioni[3])
  lenT = length(nodiT)
  
  # tratti mancanti
  n - lenR - lenP - lenT    # nodi mancanti che devo assegnare
  
  tratto1 = 1:(nodiP[1]-1)  # nodi tra 1 e curva P
  len1 = length(tratto1)
  
  tratto2 = (nodiP[lenP]+1):(nodiR[1]-1)  # nodi tra P e R
  len2 = length(tratto2)
  
  tratto3 = (nodiR[lenR]+1):(nodiT[1]-1)  # nodi tra R e T
  len3 = length(tratto3)
  
  tratto4 = (nodiT[lenT]+1):(length(HB_mean))  # nodi tra T e end
  len4 = length(tratto4)
  
  other = c(tratto1,tratto2,tratto3,tratto4)  
  length(other)
  
  interval = (length(other) / (n*(1-frazioni[1]-frazioni[2]-frazioni[3])))  # splitto uniformemente i restanti valori
  
  nodi_other = other[seq(1, length(other), interval)]
  nodi_other=c(nodi_other,end)
  length(nodi_other)
  
  knots = sort(c(nodiR,nodiP,nodiT,nodi_other))
  
  # rimuovo knots troppo vicini
  indeces = rep(0, length(knots)-1)
  
  eps = nodiR[2]-nodiR[1]
  
  for (i in 1:(length(knots)-1)){
    if (knots[i+1] - knots[i] < eps) indeces[i] = i
  }
  
  if (length(which(indeces > 0)) > 0) {
    knots = knots[-indeces]
  }
  
  
  # aggiungo un knot anche nel picco
  knots = sort(c(knots,rpeak))
  print(paste('Number of knots: ',length(knots)))
  
  return(knots)
}


smoothing <- function(signal, rpeak, punti_medi) {
  
  time = 1:length(signal)
  
  # generalized cross-validation
  nknots <- 10:25
  gcv <- numeric(length(nknots))
  for (i in 1:length(nknots)){
    nodi <- find_knots(n = nknots[i], lenHB = length(signal), rpeak = rpeak, punti_medi = punti_medi)
    basis <- create.bspline.basis(c(0,length(signal)), breaks = nodi, norder = 3)
    gcv[i] <- smooth.basis(time, signal, basis)$gcv
  }
  par(mfrow=c(1,1))
  plot(nknots,gcv)
  opt_n_knots = nknots[which.min(gcv)]
  print('Optimal:')
  
   
  nodi <- find_knots(n = opt_n_knots, lenHB = length(signal), rpeak = rpeak, punti_medi = punti_medi)
  
  # plot dei nodi
  # plot del segnale con punti di intersse e knots
  x11()
  plot(signal, lwd = 2, col = 'black', type = 'l', ylim = c(-0.2,0.5))
  points(nodi,rep(-0.2,length(nodi)), pch = 4)
  abline(h = -0.2, col = 'red')
  abline(v = c(punti_medi[c(2,3,5,7,9,10)]+rpeak), lty = 2, col = 'red')
  
  # smoothin
  basis <- create.bspline.basis(c(0,length(signal)), breaks = nodi, norder = 3)
  x11()
  plot(basis)
  
  basismat <- smooth.basis(argvals=time, y=signal, fdParobj=basis)
  smooth <- eval.fd(time, basismat$fd)
  
  #first derivative
  smooth_1 <- eval.fd(time, basismat$fd, Lfd=1)
  
  x11()
  plot(time,signal,xlab="t")
  points(time,smooth ,type="l",col="green",lwd=2)
  title('Smoothed data')
  legend('topright',c('Observed data','Smoothed data'), col = c('black','green'),lty=c(NA,1),pch = c(1,NA))
  
  # plot(time,smooth_1,type = 'l', col = 'blue')
  
  return(list(Smoothed_data = smooth, First_derivative = smooth_1))
}



