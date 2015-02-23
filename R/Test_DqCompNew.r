# Example for 
# http://stackoverflow.com/questions/25085537/permute-groups-in-a-data-frame-r#
#

source("R/Test_Dq_fun.r")


# Test with other data
#
setwd("Simul")
bName <- "neuFish65_256"
Dq1 <- readNeutral_calcDq(paste0(bName,"T100mfSAD.txt"))

plotDq_ReplaceR(Dq1,0.2,0.04,0.0001)


# Compare replacementRate = 1 and 0.1 
#
require(plyr)
require(dplyr)

Dqq <- filter(Dq1,MortalityRate==0.2,DispersalDistance==0.04,ColonizationRate==0.0001,
              ReplacementRate==1 | ReplacementRate==0.1, q<11 & q>-11)
names(Dqq)

# Add repetitions

simbyrep <- nrow(Dqq)/3

Dqq$rep <- rep( 1:3,each=simbyrep)

plotDq(Dqq,"as.factor(paste(ReplacementRate,rep))")


adt <- pairwiseAD_Dif(Dqq,"Dq",Dqq[,c(4,10)])
adt


# # Compare using permutations
# #
# require(kSamples)
# ks <- ad.test(list(Dqq$Dq[1:35],Dqq$Dq[36:70]),method="simulated",Nsim=10000)
# ks

require(reshape2)
Dq2<- dcast(Dqq,ReplacementRate+rep~q,value.var="Dq")
require(statmod)
compareTwoGrowthCurves(Dq2$ReplacementRate,Dq2[,3:ncol(Dq2)],nsim=1000)

# Compare ReplacementRate 1 and 0.01

Dqq <- filter(Dq1,MortalityRate==0.2,DispersalDistance==0.04,ColonizationRate==0.0001,
              ReplacementRate==1 | ReplacementRate==0.01, q<11 & q>-11)

simbyrep <- nrow(Dqq)/3

Dqq$rep <- rep( 1:3,each=simbyrep)

plotDq(Dqq,"as.factor(paste(ReplacementRate,rep))")


adt <- pairwiseAD_Dif(Dqq,"Dq",Dqq[,c(4,10)])
adt

Dq2<- dcast(Dqq,ReplacementRate+rep~q,value.var="Dq")
compareTwoGrowthCurves(Dq2$ReplacementRate,Dq2[,3:ncol(Dq2)],nsim=1000)

# Visualizar ecdf en las que se basas AD test

ggplot(Dqq,aes(x = Dq,colour=paste(ReplacementRate,rep))) + stat_ecdf(geom="smooth") 
