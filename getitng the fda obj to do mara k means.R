# k means della Mara

setwd("~/UNI/NECSTCamp/StPetersbirg_git/ECGanalysis")
load("R_functions/smoothin_and_knots_func.RData")

library(fda)
library(fdakma)


# hanno tutti lunghezze differenti

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

load('lista_dei_fobj.RData')

time = seq(1,169,length.out = 500)

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
  
  time = seq(1,dim(df)[1],length.out = 500)
  
  basismat.fd = list.of.fda[[i]]
  
  smooth <- eval.fd(time, basismat.fd)
  smooth_1 <- eval.fd(time, basismat.fd, Lfd=1)
  
  Y0 = rbind(Y0,t(smooth))
  Y1 = rbind(Y1,t(smooth_1))
  
}

matplot(t(Y0), type = 'l', xlab = 'index')
matplot(t(Y1), type = 'l', xlab = 'index')

# saving these points to export them
write.csv(data.frame(Y0),file = 'Y0.csv', row.names = FALSE)
write.csv(data.frame(Y1),file = 'Y1.csv', row.names = FALSE)


#### rescaling the data to i punti noti su intervallo 1:200 ####
# punti noti: P_onset P_offset, R_onset, R_offset, T_offset
# dei punti medi sono [2,3,5,6,10]


for (i in 1:length(indeces)){
  smooth_vec_filename = paste(names[indeces[i]],'/Smooth.csv',sep = '')
  smooth_vec = read.csv(smooth_vec_filename, header = TRUE)
  rpeak = ceiling(dim(smooth_vec)[1] / 2)
  
  file = paste(names[indeces[i]],'/Punti_medi.csv',sep = '')
  punti = read.csv(file, header = TRUE) + rpeak
  
  punti = punti*200/(dim(smooth_vec)[1])
  
  write.csv(punti, file, row.names = FALSE)
}

# concateno tutti questi punti in un'unica matrice
points.matrix = t(read.csv('HB_I01/Punti_medi.csv', header = TRUE))

for (i in 2:length(indeces)){

  file = paste(names[indeces[i]],'/Punti_medi.csv',sep = '')
  punti = t(read.csv(file, header = TRUE))
  
  points.matrix = rbind(points.matrix,punti)
}

write.csv(points.matrix, file = "punti_per_smoothed_data.csv")

#### K-means ####

Y0 = read.csv('Y0.csv', header = TRUE)
Y1 = read.csv('Y1.csv', header = TRUE)


# k - mean effettivo

time = 1:500

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
fdakma_example$labels

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





kma.compare_example <- kma.compare (
  x=time, y0=Y0, y1=Y1, n.clust = 3:8, 
  warping.method = c('affine'), 
  similarity.method = 'd1.pearson',
  center.method = 'k-means',
  plot.graph=1)


# k = 5
fdakma_example <- kma(
  x=time, y0=Y0, y1 = Y1, n.clust = 6, 
  warping.method = 'affine', 
  similarity.method = 'd1.pearson',  # similarity computed as the cosine
  # between the derivative functions
  # (correlation)
  center.method = 'k-means'
  #seeds = rep(1,4) # you can give a little help to the algorithm...
)

kma.show.results(fdakma_example)


kma.compare_example <- kma.compare (
  x=time, y0=Y0, y1=Y1, n.clust = 4:7, 
  warping.method = c('affine'), 
  similarity.method = 'd0.pearson',
  center.method = 'k-means',
  plot.graph=1)

table(fdakma_example$labels,indeces)
