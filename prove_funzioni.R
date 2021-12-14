
setwd("~/UNI/NECSTCamp/StPetersbirg_git/ECGanalysis")
load("R_functions/smoothin_and_knots_func.RData")

library(fda)


df <- read.csv('HB_I03/Signals.csv', header = TRUE)
head(df)

HB_mean = rowMeans(df)

plot(HB_mean, lwd = 2, col = 'black', type = 'l')
for ( c in df ) lines( c, type="l", col = 'grey' )
lines(HB_mean, lwd = 2, col = 'black', type = 'l')

# ottengo i punti d'interesse
punti = read.csv('HB_I01/Peaks.csv', header = TRUE)
head(punti)

punti_medi = round(colMeans(punti))

# r_peak lo conosco dal file jupyter
rpeak = 90 +1

# midifico P
# rpeak + punti_medi[1:3] - 45
punti_medi[1:3] = punti_medi[1:3] - 15
# abline(v = punti_medi[1:3] + rpeak, col = 'red')

# T_peak
# abline(v = 158)
# 158 - rpeak - punti_medi[8]
punti_medi[8] = 158 - rpeak

#T_offset
# abline(v = 170)
# 170 - rpeak - punti_medi[10]
punti_medi[10] = 170 - rpeak


# --------- smoothing -------------

my_smoothin = smoothing(HB_mean,rpeak = rpeak, punti_medi = punti_medi)

