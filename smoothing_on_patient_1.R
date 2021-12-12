
setwd("~/UNI/NECSTCamp/StPetersbirg_git/ECGanalysis")

library(fda)

df1 <- read.csv('HB_I01/Signals.csv', header = TRUE)
head(df1)

# changing names of the columns
vec_of_names = numeric(dim(df1)[2])

for (i in 1:dim(df1)[2]) {
  vec_of_names[i] <- paste(names(df1)[i],"I01", sep = "_")
}
vec_of_names
names(df1) = vec_of_names


df2 <- read.csv('HB_I02/Signals.csv', header = TRUE)
head(df2)

# same here
vec_of_names = numeric(dim(df2)[2])

for (i in 1:dim(df2)[2]) {
  vec_of_names[i] <- paste(names(df2)[i],"I02", sep = "_")
}
vec_of_names
names(df2) = vec_of_names

# joining the two datasets because coming from the same patitent

df = cbind(df1, df2)
head(df)
dim(df)

# computing the mean of all the heartbeats
HB_mean = rowMeans(df)

plot(HB_mean, lwd = 2, col = 'black', type = 'l')
for ( c in df ) lines( c, type="l", col = 'grey' )
lines(HB_mean, lwd = 2, col = 'black', type = 'l')

# what if i compute the median of it?
library(matrixStats)

HB_median = rowMedians(as.matrix(df))
HB_median
lines(HB_median, lwd = 2, col = 'red', type = 'l')

# considernado che la mediana di funzioni costanti non è per forza costante forse è meglio lavorare
# con la media 


# provo a trovare media degli P_offset, P_onset, R_onset, R_peaks, 

punti1 = read.csv('HB_I01/Peaks.csv', header = TRUE)
head(punti1)

punti2 = read.csv('HB_I02/Peaks.csv', header = TRUE)
head(punti2)

punti = rbind(punti1,punti2)
head(punti)
dim(punti)

punti_medi = round(colMeans(punti))

plot(HB_mean, lwd = 2, col = 'black', type = 'l')
points(85+punti_medi,HB_mean[85 + punti_medi])


#======= SMOOTHING ==============================================
# inizio a provare un po' di smoothing
end = dim(df)[1]
time = 1:end

# inizio con le spline con una distribuzione uniforme di knots
basis <- create.bspline.basis(rangeval=c(0,end), nbasis=50, norder=4)
plot(basis)

basismat <- eval.basis(time, basis)
smooth1 <- basismat %*% lsfit(basismat, HB_mean, intercept=FALSE)$coef

plot(time,HB_mean,xlab="t",ylab="observed data")
points(time,smooth1 ,type="l",col="green",lwd=2)

# forse troppe basi, non cattura il picco R

# provo con meno basi

basis <- create.bspline.basis(rangeval=c(0,end), nbasis=40, norder=4)
plot(basis)

basismat <- eval.basis(time, basis)
smooth1 <- basismat %*% lsfit(basismat, HB_mean, intercept=FALSE)$coef

plot(time,HB_mean,xlab="t",ylab="observed data")
points(time,smooth1 ,type="l",col="green",lwd=2)


# generalized cross validation to understand the best num of basis

nbasis <- 45:60
gcv <- numeric(length(nbasis))
for (i in 1:length(nbasis)){
  basis <- create.bspline.basis(c(0,end), nbasis[i], 4)
  gcv[i] <- smooth.basis(time, HB_mean, basis)$gcv
}
par(mfrow=c(1,1))
plot(nbasis,gcv)
nbasis[which.min(gcv)]

# inserirsco il valore trovato

basis <- create.bspline.basis(rangeval=c(0,end), nbasis=59, norder=4)
plot(basis)

basismat <- eval.basis(time, basis)
smooth1 <- basismat %*% lsfit(basismat, HB_mean, intercept=FALSE)$coef

plot(time,HB_mean,xlab="t",ylab="observed data")
points(time,smooth1 ,type="l",col="green",lwd=2)

# ok picco molto buono, ma numero di basi esagerato
# ---------------------------------------------
# ora cerco altra maniera per distribuire di knots:

# divido il mio segnale in tre intervalli:
plot(HB_mean, lwd = 2, col = 'black', type = 'l')
abline(v = 85 + punti_medi[5] - 2)
abline(v = 85 + punti_medi[7] + 3)

curvaP = HB_mean[1:(85 + punti_medi[5] - 2)]
plot(curvaP) 

curvaT = HB_mean[(85 + punti_medi[7] + 3):169]
plot(curvaT)

curvaR = HB_mean[(85 + punti_medi[5] - 1):(85 + punti_medi[7] + 2)]  
plot(curvaR)

length(curvaP) + length(curvaR) + length(curvaT)
# ---------------------------------
#### trovo numero knots ideale per curvaP
time = 1:length(curvaP)
basis <- create.bspline.basis(rangeval=c(0,length(curvaP)), nbasis=13, norder=4)
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

# 18

basis <- create.bspline.basis(rangeval=c(0,length(curvaP)), nbasis=18, norder=4)
plot(basis)

basismat <- eval.basis(time, basis)
smooth <- basismat %*% lsfit(basismat, curvaP, intercept=FALSE)$coef

plot(time,curvaP,xlab="t",ylab="observed data")
points(time,smooth ,type="l",col="green",lwd=2)

# mmh

# -------curva T----------------------
#### vado di curvaT

time = 1:length(curvaT)
basis <- create.bspline.basis(rangeval=c(0,length(curvaT)), nbasis=13, norder=4)
plot(basis)

basismat <- eval.basis(time, basis)
smooth <- basismat %*% lsfit(basismat, curvaT, intercept=FALSE)$coef

plot(time,curvaT,xlab="t",ylab="observed data")
points(time,smooth ,type="l",col="green",lwd=2)

# gcv
nbasis <- 5:20
gcv <- numeric(length(nbasis))
for (i in 1:length(nbasis)){
  basis <- create.bspline.basis(c(0,length(curvaT)), nbasis[i], 4)
  gcv[i] <- smooth.basis(time, curvaT, basis)$gcv
}
par(mfrow=c(1,1))
plot(nbasis,gcv)
nbasis[which.min(gcv)]

# 13

basis <- create.bspline.basis(rangeval=c(0,length(curvaT)), nbasis=13, norder=4)
plot(basis)

basismat <- eval.basis(time, basis)
smooth <- basismat %*% lsfit(basismat, curvaT, intercept=FALSE)$coef

plot(time,curvaT,xlab="t",ylab="observed data")
points(time,smooth ,type="l",col="green",lwd=2)

#------ cruvaR -------------------
# ma in realt� sta cosa non mi serve

time = 1:length(curvaR)
basis <- create.bspline.basis(rangeval=c(0,length(curvaR)), nbasis=20, norder=4)
plot(basis)

basismat <- eval.basis(time, basis)
smooth <- basismat %*% lsfit(basismat, curvaR, intercept=FALSE)$coef

plot(time,curvaR,xlab="t",ylab="observed data")
points(time,smooth ,type="l",col="green",lwd=2)

# gcv
nbasis <- 15:25
gcv <- numeric(length(nbasis))
for (i in 1:length(nbasis)){
  basis <- create.bspline.basis(c(0,length(curvaR)), nbasis[i], 4)
  gcv[i] <- smooth.basis(time, curvaR, basis)$gcv
}
par(mfrow=c(1,1))
plot(nbasis,gcv)
nbasis[which.min(gcv)]

# 23

basis <- create.bspline.basis(rangeval=c(0,length(curvaR)), nbasis=23, norder=4)
plot(basis)

basismat <- eval.basis(time, basis)
smooth <- basismat %*% lsfit(basismat, curvaR, intercept=FALSE)$coef

plot(time,curvaR,xlab="t",ylab="observed data")
points(time,smooth ,type="l",col="green",lwd=2)

#------------ ok quindi ora distribuisco i knots in questa maniera ---------

# curvaP -> ordine 18
# curvaT -> ordine 13

n = 15
nodi = c(seq(from = 1, to = 85 + punti_medi[5] - 2, length.out = 19),
         seq(from = 85 + punti_medi[5] - 1,to = 85 + punti_medi[7] + 2,length.out = n),
         seq(from = 85 + punti_medi[7] + 3,to = length(HB_mean), length.out = 14))

length(nodi)

end = length(HB_mean)
time = 1:end

# inizio con le spline con una distribuzione uniforme di knots
basis <- create.bspline.basis(rangeval=c(1,end), breaks = nodi, norder=4)
plot(basis)

basismat <- eval.basis(time, basis)
smooth1 <- basismat %*% lsfit(basismat, HB_mean, intercept=FALSE)$coef

plot(time,HB_mean,xlab="t",ylab="observed data")
points(time,smooth1 ,type="l",col="green",lwd=2)



# provo con il gcv con questo nuovo set di nodi
nknots <- 10:30
gcv <- numeric(length(nknots))
for (i in 1:length(nknots)){
  
  #creazione del vettore di nodi
  nodi = c(seq(from = 1, to = 85 + punti_medi[5] - 2, length.out = 19),
           seq(from = 85 + punti_medi[5] - 1,to = 85 + punti_medi[7] + 2,length.out = nknots[i]),
           seq(from = 85 + punti_medi[7] + 3,to = length(HB_mean), length.out = 14))
  
  # creazione della base
  basis <- create.bspline.basis(c(1,end), breaks = nodi, norder = 4)
  gcv[i] <- smooth.basis(time, HB_mean, basis)$gcv
}
par(mfrow=c(1,1))
plot(nknots,gcv)
nknots[which.min(gcv)]

# n = 20

n = 15
nodi = c(seq(from = 1, to = 85 + punti_medi[5] - 2, length.out = 19),
         seq(from = 85 + punti_medi[5] - 1,to = 85 + punti_medi[7] + 2,length.out = n),
         seq(from = 85 + punti_medi[7] + 3,to = length(HB_mean), length.out = 14))

length(nodi)

end = length(HB_mean)
time = 1:end

# inizio con le spline con una distribuzione uniforme di knots
basis <- create.bspline.basis(rangeval=c(1,end), breaks = nodi, norder=4)
plot(basis)

basismat <- eval.basis(time, basis)
smooth1 <- basismat %*% lsfit(basismat, HB_mean, intercept=FALSE)$coef

plot(time,HB_mean,xlab="t",ylab="observed data")
points(time,smooth1 ,type="l",col="green",lwd=2)

# io lascio n = 15

# calcolo la derivata

# numerical derivative
der_obs <- (HB_mean[3:end]-HB_mean[1:(end-2)])/(time[3:end]-time[1:(end-2)])

# derivative from spline approximation
basismat_der<- eval.basis(time, basis, Lfdobj=1)
der_smooth <- basismat_der %*% lsfit(basismat_der, HB_mean, intercept=FALSE)$coef

# plot
plot(time,der_smooth ,type="l",col="blue",lwd=2)
points(time[2:(end-1)],der_obs,xlab="t",ylab="first differences x",type="l")


