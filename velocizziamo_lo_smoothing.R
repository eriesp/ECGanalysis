
setwd("~/UNI/NECSTCamp/StPetersbirg_git/ECGanalysis")
load("R_functions/smoothin_and_knots_func.RData")

library(fda)

start = 1

names = vector(mode="character", length=(75-start + 1 ))

for (i in start:75) {
  if (i < 10) names[i-start+1] = paste("HB_I0",i,sep = "")
  else names[i-start+1] = paste("HB_I",i, sep="")
}

i = 75

df <- read.csv(paste(names[i],'/Signals.csv',sep = ''), header = TRUE)
head(df)

dim(df)

if (length(which(df[,dim(df)[2]] == 0)) == dim(df)[1]) {
  df = df[-dim(df)[2]]
}
dim(df)

HB_mean = rowMeans(df)

plot(HB_mean, lwd = 2, col = 'black', type = 'l')
for ( c in df ) lines( c, type="l", col = 'grey' )
lines(HB_mean, lwd = 2, col = 'black', type = 'l')


punti = read.csv(paste(names[i],'/Peaks.csv',sep = ''), header = TRUE)
head(punti)

# punti_medi = round(colMeans(punti))

punti_medi = round(apply(punti, 2, mean, na.rm = TRUE))

# r_peak 
rpeak = ceiling(length(HB_mean) / 2)

# plot dei punti importanti
plot(HB_mean, lwd = 2, col = 'black', type = 'l')
points(rpeak+punti_medi,HB_mean[rpeak + punti_medi])
abline(v = rpeak, col = 'red')

# questo problema dei punti è ingestibile anyway 
# li shifto
abline(v = 105)
punti_medi[3] = 105 - rpeak

abline(v = 100)
punti_medi[1] = 100 - rpeak

abline(v = 45)
punti_medi[2] = 45 - rpeak

abline(v = 170)
punti_medi[9] = 170 - rpeak

abline(v = 180)
punti_medi[8] = 180 - rpeak

abline(v = 190)
punti_medi[10] = 190 - rpeak

# S peak
abline(v = 120)
punti_medi[7] = 120 - rpeak

# R offset
abline(v = 140)
punti_medi[6] = 140 - rpeak


# plot dei punti importanti
plot(HB_mean, lwd = 2, col = 'black', type = 'l')
points(rpeak+punti_medi,HB_mean[rpeak + punti_medi])
abline(v = rpeak, col = 'red')

# smoothing
my_smoothin = smoothing(HB_mean,rpeak = rpeak, punti_medi = punti_medi)
df_smooth = data.frame(my_smoothin)
write.csv(df_smooth,paste(names[i],'/Smooth.csv',sep = ''), row.names = FALSE)
