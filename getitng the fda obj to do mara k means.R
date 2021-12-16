# k means della Mara

setwd("~/UNI/NECSTCamp/StPetersbirg_git/ECGanalysis")
load("R_functions/smoothin_and_knots_func.RData")

library(fda)
library(fdakma)


# no ok non posso perché hanno tutti lunghezze differenti

# rifaccio lo smoothing di tutte le func ottenendo il fda object

start = 1

names = vector(mode="character", length=(75-start + 1 ))

for (i in start:75) {
  if (i < 10) names[i-start+1] = paste("HB_I0",i,sep = "")
  else names[i-start+1] = paste("HB_I",i, sep="")
}

indeces = c(1,2,3,4,6,8,9,13,15,19,21,23,26,27,29,33,34,35,37,38,39,40,
            41,42,45,47,49,50,52,53,54,55,56,57,58,59,63,64,66,68,70,72)

length(indeces)



i = 42



df <- read.csv(paste(names[indeces[i]],'/Smooth.csv',sep = ''), header = TRUE)
dim(df)

y0 = df$Smoothed_data

plot(y0, lwd = 2, col = 'black', type = 'l')

punti = read.csv(paste(names[indeces[i]],'/Peaks.csv',sep = ''), header = TRUE)
# head(punti)

# punti_medi = round(colMeans(punti))

punti_medi = round(apply(punti, 2, mean, na.rm = TRUE))

# r_peak 
rpeak = ceiling(length(y0) / 2)

# plot dei punti importanti
plot(y0, lwd = 2, col = 'black', type = 'l')
points(rpeak+punti_medi[c(1,2,3,5,6,8,9,10)],y0[rpeak + punti_medi[c(1,2,3,5,6,8,9,10)]], col = 'red')
abline(v = rpeak, col = 'red')

# questo problema dei punti è ingestibile anyway 
# li shifto
k = 65
abline(v = k)
punti_medi[3] = k - rpeak

k = 100
abline(v = k)
punti_medi[1] = k - rpeak

k = 40
abline(v = k)
punti_medi[2] = k - rpeak

k = 165
abline(v = k)
punti_medi[9] = k - rpeak

k = 180
abline(v = k)
punti_medi[8] = k - rpeak

k = 200
abline(v = k)
punti_medi[10] = k - rpeak

# R onset
k = 115
abline(v = k)
punti_medi[5] = k - rpeak

# R offset
k = 120
abline(v = k)
punti_medi[6] = k - rpeak


plot(y0, lwd = 2, col = 'black', type = 'l')
points(rpeak+punti_medi[c(1,2,3,5,6,8,9,10)],y0[rpeak + punti_medi[c(1,2,3,5,6,8,9,10)]], col = 'red')
abline(v = rpeak, col = 'red')

write.csv(data.frame(punti_medi),paste(names[indeces[i]],'/Punti_medi.csv',sep = ''), row.names = FALSE)

load("R_functions/smoothin_modified.RData")
# smoothing
smoothing_modified = smoothing_modified(y0,rpeak = rpeak, punti_medi = punti_medi)
# smoothing_modified

# list.of.fda = list(i = smoothing_modified)
list.of.fda = c(list.of.fda, list(smoothing_modified))
length(list.of.fda)


### ora che ho gli fobj di tutti li valuto su intervallo 1:200
# e li aggiungo tutti alla stessa matrice

time = seq(1,169,length.out = 200)

basismat.fd = list.of.fda[1]$i

smooth <- eval.fd(time, basismat.fd)
smooth_1 <- eval.fd(time, basismat.fd, Lfd=1)

plot(time,smooth,type = 'l')
plot(time,smooth_1,type = 'l')

# metto in una matrice tutte le funzioni che ho smoothate
Y0 = t(smooth)
Y1 = t(smooth_1)

for (i in 2:length(list.of.fda)){
  file = paste(names[indeces[i]],'/Smooth.csv',sep = '')
  df = read.csv(file, header = TRUE)
  
  time = seq(1,dim(df)[1],length.out = 200)
  
  basismat.fd = list.of.fda[[i]]
  
  smooth <- eval.fd(time, basismat.fd)
  smooth_1 <- eval.fd(time, basismat.fd, Lfd=1)
  
  Y0 = rbind(Y0,t(smooth))
  Y1 = rbind(Y1,t(smooth_1))
  
}

matplot(t(Y0), type = 'l', xlab = 'index')


# k - mean effettivo

time = 1:200

fdakma_example <- kma(
  x=time, y0=Y0, n.clust = 3, 
  warping.method = 'affine', 
  similarity.method = 'd0.pearson',  # similarity computed as the cosine
  # between the functions
  # (correlation)
  center.method = 'k-means'
  #seeds = c(1,21) # you can give a little help to the algorithm...
)

kma.show.results(fdakma_example)

# k = 5
fdakma_example <- kma(
  x=time, y0=Y0, n.clust = 5, 
  warping.method = 'affine', 
  similarity.method = 'd0.pearson',  # similarity computed as the cosine
  # between the functions
  # (correlation)
  center.method = 'k-means'
  #seeds = c(1,21) # you can give a little help to the algorithm...
)

kma.show.results(fdakma_example)


# k = 8
fdakma_example <- kma(
  x=time, y0=Y0, n.clust = 8, 
  warping.method = 'affine', 
  similarity.method = 'd0.pearson',  # similarity computed as the cosine
  # between the functions
  # (correlation)
  center.method = 'k-means'
  #seeds = c(1,21) # you can give a little help to the algorithm...
)

kma.show.results(fdakma_example)

# con k = 8 la varibilità non migliora


matplot(t(Y1), type = 'l', xlab = 'index')

# k = 5
fdakma_example <- kma(
  x=time, y0=Y0, y1 = Y1, n.clust = 5, 
  warping.method = 'affine', 
  similarity.method = 'd1.pearson',  # similarity computed as the cosine
  # between the derivative functions
  # (correlation)
  center.method = 'k-means'
  #seeds = rep(1,4) # you can give a little help to the algorithm...
)

kma.show.results(fdakma_example)
