# provo a distribuire io i nodi per tutti

# devono esserce un po' attorno a P
# un po' attorno a T
# molti attorno a R

end = length(HB_mean)
time = 1:end

n = 17 # a caso 

frazioni = c('R' = 2/5, 'P' = 1/6, 'T' = 1/6)

# knots in curva R
nodiR = seq(from = rpeak + punti_medi[5] - 1,to = rpeak + punti_medi[7] + 2,length.out = n*frazioni[1])
lenR = length(nodiR)

# knots in curva P
nodiP = seq(from = rpeak + punti_medi[2] - 3,to = rpeak + punti_medi[3] + 3,length.out = n*frazioni[2])
lenP = length(nodiP)

# knots in curva T
nodiT = seq(from = rpeak + punti_medi[9] - 3,to = rpeak + punti_medi[10] + 3,length.out = n*frazioni[3])
lenT = length(nodiT)

n - lenR - lenP - lenT    # nodi mancanti che devo assegnare

tratto1 = 1:(nodiP[1]-1)  # nodi tra 1 e curva P
len1 = length(tratto1)

tratto2 = (nodiP[lenP]+1):(nodiR[1]-1)  # nodi tra P e R
len2 = length(tratto2)

tratto3 = (nodiR[lenR]+1):(nodiT[1]-1)  # nodi tra R e T
len3 = length(tratto3)

tratto4 = (nodiT[lenT]+1):(length(HB_mean))  # nodi tra T e end
len4 = length(tratto4)

other = c(tratto1,tratto2,tratto3,tratto4)  
length(other)


interval = (length(other) / (n*(1-frazioni[1]-frazioni[2]-frazioni[3])))  # splitto uniformemente i restanti valori

nodi_other = other[seq(1, length(other), interval)]
nodi_other=c(nodi_other,end)
length(nodi_other)

knots = sort(c(nodiR,nodiP,nodiT,nodi_other))
length(knots)


plot(HB_mean, lwd = 2, col = 'black', type = 'l', ylim = c(-0.2,0.5))
points(knots,rep(-0.2,length(knots)), pch = 4)
abline(h = -0.2, col = 'red')
abline(v = c(punti_medi[c(2,3,5,7,9,10)]+rpeak), lty = 2, col = 'red')

indeces = rep(0, length(knots)-1)

eps = nodiR[2]-nodiR[1]

for (i in 1:(length(knots)-1)){
  if (knots[i+1] - knots[i] < eps) indeces[i] = i
}

if (length(which(indeces > 0)) > 0) {
  knots = knots[-indeces]
}

knots = sort(c(knots,rpeak))
length(knots)


# preR = knots[knots < punti_medi[5] + rpeak]
# postR = knots[knots > punti_medi[7] + rpeak]
# 
# for (i in 1:(length(preR)-1)){
#   if (preR[i+1] - preR[i] < eps){ preR = preR[preR != preR[i]]}
#     
# }
# 
# for (i in 1:(length(postR)-1)){
#   if (postR[i+1] - postR[i] < eps) {postR = postR[postR != postR[i]]}
#     
# }
# 
# knots = sort(c(preR,postR,nodiR))
# length(knots)

plot(HB_mean, lwd = 2, col = 'black', type = 'l', ylim = c(-0.2,0.5))
points(knots,rep(-0.2,length(knots)), pch = 4)
abline(h = -0.2, col = 'red')
abline(v = c(punti_medi[c(2,3,5,7,9,10)]+rpeak), lty = 2, col = 'red')


# proviamo lo smoothing con questi nodi
basis <- create.bspline.basis(rangeval=c(1,end), breaks = knots, norder=3)
plot(basis)

Xsp <- smooth.basis(argvals=time, y=HB_mean, fdParobj=basis)
Xsp0 <- eval.fd(time, Xsp$fd)

plot(time,HB_mean,xlab="t",ylab="observed data")
points(time,Xsp0 ,type="l",col="green",lwd=2)


# prova funzionamento funzione
nodi = find_knots(25,end,rpeak,punti_medi)
nodi
