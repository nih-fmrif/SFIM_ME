source('/usr/local/abin/AFNIio.R')
#myBetas<-read.AFNI('/data/SFIMJGC/TALK_fMRIClassME/PrcsData/SBJ02/D03_MEICA/SBJ02_S02Run10.chComp.EXTRA.Beta021+orig')
#TE<-c(15.4,29.7,44)
#xS0<-c(1,1,1)
myBetas<-read.AFNI('/data/SFIMJGC/TALK_fMRIClassME/PrcsData/SBJ03T01/D03_MEICA/SBJ03T01_RestE5A3.chComp.EXTRA.Beta030+orig')
TE<-c(11.1,23.5,35.9,48.3,60.7)
xS0<-c(1,1,1,1,1)

names(myBetas)
#I=32
#J=55
#K=8
I=32
J=52
K=9
#I=20
#J=19
#K=8
Betas<-myBetas$brk[I+1,J+1,K+1,]

xR2<-TE/mean(TE)

summary(lm(Betas~0+xR2))
summary(lm(Betas~0+xS0))

FM_results<-lm(Betas~xR2)
FM_slope<-FM_results$coefficients[2]
FM_inter<-FM_results$coefficients[1]
Beta_estimate = (FM_slope*xR2)+FM_inter
plot(TE,Betas,col='green')
lines(TE,Betas,col='green')
points(TE,Beta_estimate,col='red')
lines(TE,Beta_estimate,col='red')
summary(lm(Betas~xR2))