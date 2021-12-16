
setwd("~/UNI/NECSTCamp/StPetersbirg_git/ECGanalysis")
load("R_functions/smoothin_and_knots_func.RData")

library(pracma)
library(fda)
library(fdakma)    

start = 1

names = vector(mode="character", length=(75-start + 1 ))

for (i in start:75) {
  if (i < 10) names[i-start+1] = paste("HB_I0",i,sep = "")
  else names[i-start+1] = paste("HB_I",i, sep="")
}

i = 1

df <- read.csv(paste(names[i],'/Smooth.csv',sep = ''), header = TRUE)
head(df)

y0 = df$Smoothed_data
y1 = df$First_derivative

time = 1:length(y0)

plot(time, y0, type = 'l')


# provo a trovare i picchi di questa func perché gli altri che avevo fatto manualmente sono andati persi
picchi = findpeaks(y0)

plot(y0, type = 'l')
points(picchi[,2],picchi[,1], pch = 3, col = 'red')

picchi = findpeaks(-y0)
plot(y0, type = 'l')
points(picchi[,3:4],y0[picchi[,3:4]], pch = 3, col = 'red')

peaks = read.csv(paste(names[i],'/Peaks.csv',sep = ''), header = TRUE)

punti_medi = round(apply(peaks, 2, mean, na.rm = TRUE))

# r_peak 
rpeak = ceiling(length(y0) / 2)

plot(y0, lwd = 2, col = 'black', type = 'l')
points(rpeak+punti_medi[c(2,3,5,6,10)],y0[rpeak + punti_medi[c(2,3,5,6,10)]], pch = 3, col = 'red')
abline(v = rpeak, col = 'red')


x11()
par(mfrow = c(2,5))
count = 1
i = 39

while (count <= 10){
  file = paste(names[i],'/Smooth.csv',sep = '')
  if (file.exists(file)){
    df = read.csv(file, header = TRUE)
    
    y0 = df$Smoothed_data
    
    rpeak = ceiling(length(y0) / 2)
    
    peaks = read.csv(paste(names[i],'/Peaks.csv',sep = ''), header = TRUE)
    
    punti_medi = round(apply(peaks, 2, mean, na.rm = TRUE))
    
    plot(y0, lwd = 2, col = 'black', type = 'l')
    points(rpeak+punti_medi[c(1,2,3,5,6,10)],y0[rpeak + punti_medi[c(1,2,3,5,6,10)]], pch = 3, col = 'red')
    abline(v = rpeak, col = 'red')
    
    i = i+1
    count = count + 1 
  }
  else {
    i = i +1
  }
}


# prendo di fissi P_onset, R onset, R offset e fine e faccio il landpmark 

signal = y0
time = 1:length(y0)
rpeak = ceiling(length(y0) / 2)

nodi <- find_knots(n = 25, lenHB = length(signal), rpeak = rpeak, punti_medi = punti_medi)

basis <- create.bspline.basis(c(0,length(signal)), breaks = nodi, norder = 3)
x11()
plot(basis)

basismat <- smooth.basis(argvals=time, y=signal, fdParobj=basis)
smooth <- eval.fd(time, basismat$fd)

#first derivative
smooth_1 <- eval.fd(time, basismat$fd, Lfd=1)

plot(y0)
lines(smooth, type = 'l', col = 'green', lwd = 2)

#### LANDMARK REGISTRATION  #####
i = 1

peaks = read.csv(paste(names[i],'/Peaks.csv',sep = ''), header = TRUE)
punti_medi = round(apply(peaks, 2, mean, na.rm = TRUE))

landmark = punti_medi[c(2,5,6)]

i = i +1
while (i <= 75){
  file = paste(names[i],'/Smooth.csv',sep = '')
  if (file.exists(file)){
    
    peaks = read.csv(paste(names[i],'/Peaks.csv',sep = ''), header = TRUE)
    punti_medi = round(apply(peaks, 2, mean, na.rm = TRUE))
    
    landmark = rbind(landmark,punti_medi[c(2,5,6)])
    
    i = i+1
  }
  else {
    i = i +1
  }
}

dim(landmark)

# smoothing e registration

signal = y0
time = 1:length(y0)
rpeak = ceiling(length(y0) / 2)

nodi <- find_knots(n = 25, lenHB = length(signal), rpeak = rpeak, punti_medi = punti_medi)

wbasis <- create.bspline.basis(c(0,length(signal)), breaks = nodi, norder = 3)

Wfd0   <- fd(matrix(0,wbasis$nbasis,1),wbasis)
WfdPar <- fdPar(Wfd0, 1, 1e-4)

fdobj   <- smooth.basis(argvals, densY, wbasis,
                        fdnames = c("x", "samples", "density"))$fd

regDens   <- landmarkreg(fdobj, landmark[1:k,], WfdPar=WfdPar, monwrd=TRUE)

