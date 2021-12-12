
setwd("~/UNI/NECSTCamp/StPetersbirg_git/ECGanalysis")

library(fda)
library(dplyr)


df <- read.csv('HB_I03/Signals.csv', header = TRUE)
head(df)

# HB_mean = df['mean']

df %>% select(-mean)

HB_mean = rowMeans(df)



plot(HB_mean, lwd = 2, col = 'black', type = 'l')
for ( c in df ) lines( c, type="l", col = 'grey' )
lines(HB_mean, lwd = 2, col = 'black', type = 'l')


punti = read.csv('HB_I01/Peaks.csv', header = TRUE)
head(punti)

punti_medi = round(colMeans(punti))

rpeak = 90 +1

plot(HB_mean, lwd = 2, col = 'black', type = 'l')
points(rpeak+punti_medi,HB_mean[rpeak + punti_medi])
abline(v = rpeak, col = 'red')

# provo a spostare punti messi very male
abline(v = 45) # questo deve essere il nuovo P_peak
rpeak + punti_medi[1:3] - 45
punti_medi[1:3] = punti_medi[1:3] - 15
abline(v = punti_medi[1:3] + rpeak, col = 'red')

# T_peak
abline(v = 158)
158 - rpeak - punti_medi[8]
punti_medi[8] = 158 - rpeak

#T_offset
abline(v = 170)
170 - rpeak - punti_medi[10]
punti_medi[10] = 170 - rpeak


# ok punti sistemati 

#------- Smoothing -----
end = length(HB_mean)
time = 1:end

plot(HB_mean, lwd = 2, col = 'black', type = 'l')
abline(v = rpeak + punti_medi[5] - 2)
abline(v = rpeak + punti_medi[7] + 3)

curvaP = HB_mean[1:(rpeak + punti_medi[5] - 2)]
plot(curvaP) 

curvaT = HB_mean[(rpeak + punti_medi[7] + 3):end]
plot(curvaT)

curvaR = HB_mean[(rpeak + punti_medi[5] - 1):(rpeak + punti_medi[7] + 2)]  
plot(curvaR)

length(curvaP) + length(curvaR) + length(curvaT)

##### curva P ####-------------------
#### trovo numero knots ideale per curvaP
time = 1:length(curvaP)
basis <- create.bspline.basis(rangeval=c(0,length(curvaP)), nbasis=18, norder=4)
plot(basis)

basismat <- eval.basis(time, basis)
smooth <- basismat %*% lsfit(basismat, curvaP, intercept=FALSE)$coef

plot(time,curvaP,xlab="t",ylab="observed data")
points(time,smooth ,type="l",col="green",lwd=2)

# gcv
nbasis <- 5:20
gcv <- numeric(length(nbasis))
for (i in 1:length(nbasis)){
  basis <- create.bspline.basis(c(0,length(curvaP)), nbasis[i], 4)
  gcv[i] <- smooth.basis(time, curvaP, basis)$gcv
}
par(mfrow=c(1,1))
plot(nbasis,gcv)
nbasis[which.min(gcv)]

# dato che il numero di basi rimane alto provo a fare smoothing con penalizzazione
basis <- create.bspline.basis(rangeval=c(0,length(curvaP)), nbasis=18, norder=4)

# lambda <- c(1e-2,1e-1,1,1e1,1e2)
lambda <- c(50,10,1,1e-1,1e-2,1e-3,1e-4,1e-5)
gcv <- numeric(length(lambda))
for (i in 1:length(lambda)){
  functionalPar <- fdPar(fdobj=basis, Lfdobj=1, lambda=lambda[i])  
  gcv[i] <- smooth.basis(time, curvaP, functionalPar)$gcv
}
par(mfrow=c(1,1))
plot(log10(lambda),gcv)
lam_gcv = lambda[which.min(gcv)]
lam_gcv

functionalPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=lam_gcv)  # I want to penalize the first derivative
# functional parameter, having arguments: 
# basis, order of the derivative to be penalized, smoothing parameter.

Xss <- smooth.basis(time, curvaP, functionalPar)
Xss$df

Xss0 <- eval.fd(time, Xss$fd, Lfd=0)

plot(time,curvaP,xlab="t",ylab="observed data")
points(time, Xss0 ,type="l",col="green",lwd=2)



# famo che tenemo 12 basi e basta in questo intervallo

##### curva T ####-------------------
time = 1:length(curvaT)
basis <- create.bspline.basis(rangeval=c(0,length(curvaT)), nbasis=9, norder=4)
plot(basis)

Xsp <- smooth.basis(argvals=time, y=curvaT, fdParobj=basis)
Xsp0 <- eval.fd(time, Xsp$fd)

plot(time,curvaT,xlab="t",ylab="observed data")
points(time,Xsp0 ,type="l",col="green",lwd=2)

nbasis <- 5:20
gcv <- numeric(length(nbasis))
for (i in 1:length(nbasis)){
  basis <- create.bspline.basis(c(0,length(curvaT)), nbasis[i], 4)
  gcv[i] <- smooth.basis(time, curvaT, basis)$gcv
}
par(mfrow=c(1,1))
plot(nbasis,gcv)
nbasis[which.min(gcv)]



