setwd("~/UNI/NECSTCamp/StPetersburgINCART")

# dataset con i battiti
df1 <- read.csv('HB_I01/Signals.csv', header = TRUE)
head(df1)

# dataset con i punti di intresse
peaks = read.csv('HB_I01/Peaks.csv', header = TRUE)
head(peaks)

matplot(df1['X0'], type = 'l')
abline(v = 85, col = "red")
abline(v = peaks[1,1] + 85, col = 'blue')
abline(v = peaks[1,4] + 85, col = 'green')
abline(v = peaks[1,7] + 85, col = 'yellow')
abline(v = peaks[1,8] + 85, col = 'pink')

time = 1:169

library(fda) 

firstHB = df1$X0


basis <- create.fourier.basis(rangeval=c(0,169),nbasis=31)
data_HB.fb <- Data2fd(y = firstHB,argvals = time,basisobj = basis)
plot.fd(data_HB.fb)
points(time,firstHB ,type="l",col="blue",lwd=1)

