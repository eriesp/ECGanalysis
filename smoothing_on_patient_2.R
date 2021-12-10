

setwd("~/UNI/NECSTCamp/StPetersbirg_git/ECGanalysis")

library(fda)

df1 <- read.csv('HB_I03/Signals.csv', header = TRUE)
head(df1)

# changing names of the columns
vec_of_names = numeric(dim(df1)[2])

for (i in 1:dim(df1)[2]) {
  vec_of_names[i] <- paste(names(df1)[i],"I03", sep = "_")
}
vec_of_names
names(df1) = vec_of_names


df2 <- read.csv('HB_I04/Signals.csv', header = TRUE)
head(df2)

# same here
vec_of_names = numeric(dim(df2)[2])

for (i in 1:dim(df2)[2]) {
  vec_of_names[i] <- paste(names(df2)[i],"I04", sep = "_")
}
vec_of_names
names(df2) = vec_of_names

df3 <- read.csv('HB_I05/Signals.csv', header = TRUE)
head(df2)

# same here
vec_of_names = numeric(dim(df2)[2])

for (i in 1:dim(df2)[2]) {
  vec_of_names[i] <- paste(names(df2)[i],"I05", sep = "_")
}
vec_of_names
names(df2) = vec_of_names

# joining the two datasets because coming from the same patitent

df = cbind(df1, df2, df3)
head(df)
dim(df)

# computing the mean of all the heartbeats
HB_mean = rowMeans(df)

plot(HB_mean, lwd = 2, col = 'black', type = 'l')
for ( c in df ) lines( c, type="l", col = 'grey' )
lines(HB_mean, lwd = 2, col = 'black', type = 'l')