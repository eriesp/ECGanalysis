smoothing_modified <- function(signal, rpeak, punti_medi) {
  
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
  plot(signal, lwd = 2, col = 'black', type = 'l', ylim = c(-0.2,1))
  points(nodi,rep(-0.2,length(nodi)), pch = 4)
  abline(h = -0.2, col = 'red')
  abline(v = c(punti_medi[c(2,3,5,6,9,10)]+rpeak), lty = 2, col = 'red')
  
  # smoothin
  basis <- create.bspline.basis(c(0,length(signal)), breaks = nodi, norder = 3)
  
  basismat <- smooth.basis(argvals=time, y=signal, fdParobj=basis)
  smooth <- eval.fd(time, basismat$fd)
  
  
  x11()
  plot(time,signal,xlab="t")
  points(time,smooth ,type="l",col="green",lwd=2)
  title('Smoothed data')
  legend('topright',c('Observed data','Smoothed data'), col = c('black','green'),lty=c(NA,1),pch = c(1,NA))
  
  # plot(time,smooth_1,type = 'l', col = 'blue')
  
  return(fdaObj = basismat$fd)
}

save(smoothing_modified, file = "R_functions/smoothin_modified.RData")