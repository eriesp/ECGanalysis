
setwd("~/UNI/NECSTCamp/StPetersbirg_git/ECGanalysis")

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

# ora cerco altra maniera per distribuire di knots:
# se n è il numero di knots
# ne metto almeno 3 per la curva P
# ne metto almeno 3 per la curva T
# ne metto almeno n/3 per QRT

# no aspe
# se li metto esponenzialmente distanti dal picco r?

plot(time, HB_mean, type = 'l')
abline(v = c(85-1,85+1), col = 'red')
abline(v = c(85-2,85+2), col = 'red')
abline(v = c(85-5,85+5), col = 'orange')
abline(v = c(85-10,85+10), col = 'yellow')



abline(v = c(85-25, 85-30,85), col = 'pink')
abline(v = c(85-25, 85+25), col = 'pink')
abline(v = c(85-25, 85+25), col = 'pink')


# ora devo trovare come distribuirli in base al numero tot di knots che voglio
