# Example for 
# http://stackoverflow.com/questions/25085537/permute-groups-in-a-data-frame-r#
#
#DqStr <- "Group   q        Dq       SD.Dq
#    1 -3.0 0.7351 0.0067
#    1 -2.5 0.6995 0.0078
#    1 -2.0 0.6538 0.0093
#    2 -3.0 0.7203 0.0081
#    2 -2.5 0.6829 0.0094
#    2 -2.0 0.6350 0.0112"
#    
#    Dq1 <- read.table(textConnection(DqStr), header=TRUE)
#    
#    g <-unique(Dq1$q)
#    Dq2<- data.frame()
#    for(n in g)
#    {
#      Dqq <- Dq1[Dq1$q==n,]
#      Dqq$Group <-sample(Dqq$Group)
#      Dq2 <- rbind(Dq2,Dqq)    
#    }
#
#    n <- -3
#    print(n)
#    Dqq <- Dq1[Dq1$q==g,]
#    Dqq$Group <-sample(Dqq$Group)
#    Dq2 <- rbind(Dq2,Dqq)    
#    print(Dq2[,1:2])
#    library(plyr)
#    ddply(Dq1,.(q), function(x) { x$Group <- sample(x$Group)
#                                  data.frame(x)})
#
#Dq1$Group <- with(Dq1, ave(Group, q, FUN = sample))
#Dq1

# Test for comparing two curves with sd
#

source("R/Test_Dq_fun.r")


# Test with other data
#
setwd("Simul")
Dq1 <- readNeutral_calcDq(paste0(bName,"T100mfSAD.txt"))
Dq1$Time <- rep( 1:(nrow(Dq1)/35),each=35)
Dq1$n <- 6

plotDq_ReplaceR(Dq1,0.2,0.04,0.0001)

Dqq <- sel_ReplaceR(Dq1,0.2,0.04,0.001,1)
names(Dqq)


Dqq<- Dqq[1:70,c(5:8,10)]
plotDq(Dqq,"as.factor(Time)")

#undebug(compareTwoGCurves)

# Compare using permutations
#
compareTwoGCurves(Dqq$Time,Dqq[,3:5],Dqq$q,1000)


Dq2 <- sel_ReplaceR(Dq1,0.2,0.04,0.001,0.01)
Dq2<- Dq2[1:70,c(5:8,10)]
plotDq(Dq2,"as.factor(Time)")

Dq2 <- rbind(Dq2[1:35,],Dqq[1:35,])
plotDq(Dq2,"as.factor(Time)")
names(Dq2)

compareTwoGCurves(Dq2$Time,Dq2[,3:5],Dq2$q,1000)


# Compare using differences an Bartlett test
#

g <- unique(Dq2$Time)
Diq1 <- with(Dq2, Dq2[Time==g[1],])
Diq2 <- with(Dq2, Dq2[Time==g[2],]) 
diffNormT(Diq1[,3:5],Diq2[,3:5])


g <- unique(Dqq$Time)
Diq1 <- with(Dqq, Dqq[Time==g[1],])[,3:5]
Diq2 <- with(Dqq, Dqq[Time==g[2],])[,3:5]
diffNormT(Diq1,Diq2)


# Compare using monte carlo
#
undebug(compareTwoCurvesMC)
compareTwoCurvesMC(Dqq$Time,Dqq$Dq,Dqq$SD.Dq)
compareTwoCurvesMC(Dq2$Time,Dq2$Dq,Dq2$SD.Dq,10000)
Diq2 <- with(Dq2, Dq2[q<0,])
compareTwoCurvesMC(Diq2$Time,Diq2$Dq,Diq2$SD.Dq,10000)

# Another dataset

Dq1 <- readNeutral_calcDq(paste0(bName,"T100mfOrd.txt"))
Dq1$Time <- rep( 1:(nrow(Dq1)/35),each=35)
Dq1$n <- 6

plotDq_ReplaceR(Dq1,0.2,0.04,0.0001)

Dqq <- sel_ReplaceR(Dq1,0.2,0.04,0.001,1)
names(Dqq)


Dqq<- Dqq[1:70,c(5:8,10)]
plotDq(Dqq,"as.factor(Time)")

#undebug(compareTwoGCurves)

compareTwoGCurves(Dqq$Time,Dqq[,3:5],Dqq$q,1000)


Dq2 <- sel_ReplaceR(Dq1,0.2,0.04,0.001,0.01)
Dq2<- Dq2[1:70,c(5:8,10)]
plotDq(Dq2,"as.factor(Time)")

Dq2 <- rbind(Dq2[1:35,],Dqq[1:35,])
plotDq(Dq2,"as.factor(Time)")


compareTwoGCurves(Dq2$Time,Dq2[,3:5],Dq2$q,1000)

Dq2$SD.Dq <- Dq2$SD.Dq/10 
plotDq(Dq2,"as.factor(Time)")
compareTwoGCurves(Dq2$Time,Dq2[,3:5],Dq2$q,1000)


g <- unique(Dqq$Time)

Diq1 <- with(Dqq, Dqq[Time==g[1],])
Diq2 <- with(Dqq, Dqq[Time==g[2],]) 
diffNormT(Diq1[,3:5],Diq2[,3:5])


g <- unique(Dq2$Time)

Diq1 <- with(Dq2, Dq2[Time==g[1],])
Diq2 <- with(Dq2, Dq2[Time==g[2],]) 


names(Dq2)
diffNormT(Diq1[,3:5],Diq2[,3:5])



Dq2$SD.Dq[is.nan(Dq2$SD.Dq)] <- 0
c <- compareTwoCurvesMC(Dq2$Time,Dq2$Dq,Dq2$SD.Dq,1000)
hist(c$dist)
c$stat
c
debug(compareTwoCurvesMC)

## COMPARE USING MULTIPLE T TEST


debug(compareTwoCurvesT)
c <- compareTwoCurvesT(Dq2$Time,Dq2$Dq,Dq2$SD.Dq,6)
Dq2$p.value <- c
Dq2[Dq2$p.value<0.05,]

Fisher.test(na.omit(c))  # p-value = 0.017

c <- compareTwoCurvesT(Dqq$Time,Dqq$Dq,Dqq$SD.Dq,6)
c[c<0.05]

Fisher.test(na.omit(c))  # p-value = 0.017

c[c<0.05]




### TEST KS for Dq comparison  
#
# THIS Discard the SD but works the way it should!
#
Dq1 <- readNeutral_calcDq(paste0(bName,"T100mfSAD.txt"))


names(Dq1)
plotDq_ReplaceR(Dq1,0.2,0.04,0.0001)

Dqq <- sel_ReplaceR(Dq1,0.2,0.04,0.001,1)
Dqq$rep <- rep( 1:(nrow(Dqq)/35),each=35)

mk1 <- pairKS_DqByRep(Dqq$Dq,Dqq$rep)
propNotDiffSAD(mk1) # =.36

Dqq <- sel_ReplaceR(Dq1,0.2,0.04,0.001,0.1)
Dqq$rep <- rep( 1:(nrow(Dqq)/35),each=35)
mk1 <- pairKS_DqByRep(Dqq$Dq,Dqq$rep)
propNotDiffSAD(mk1) # =1

Dqq <- sel_ReplaceR(Dq1,0.2,0.04,0.001,0.01)
Dqq$rep <- rep( 1:(nrow(Dqq)/35),each=35)
mk1 <- pairKS_DqByRep(Dqq$Dq,Dqq$rep)
propNotDiffSAD(mk1) # =1


Dqq <- sel_ReplaceR(Dq1,0.2,0.04,0.001,0.001)
Dqq$rep <- rep( 1:(nrow(Dqq)/35),each=35)
mk1 <- pairKS_DqByRep(Dqq$Dq,Dqq$rep)
propNotDiffSAD(mk1) # =1

Dqq <- sel_ReplaceR(Dq1,0.2,0.04,0.001,0)
Dqq$rep <- rep( 1:(nrow(Dqq)/35),each=35)
mk1 <- pairKS_DqByRep(Dqq$Dq,Dqq$rep)
propNotDiffSAD(mk1) # =1

# Para DqSAD solamente replacementRate=1 da que los que provienen de simulaciones iguales 
# son distintos

# Dq1 <- readNeutral_calcDq(paste0(bName,"T100mfOrd.txt"))

# RELACIONAR VARAIBILIDADES TEMPORALES - ESPACIALES con PARAMETROS
# Como varia SAD DqSAD DqSRS en tiempo en relacion de con parametros/procesos 
# (nicho vs no nicho)
#


require(plyr)
Dq2 <- ddply(Dq1,.(MortalityRate,DispersalDistance,ColonizationRate,ReplacementRate,q),summarise,Dq=mean(Dq))

plotDq_ReplaceR(Dq2,0.2,0.04,0.0001)


# Test pairwise diferences in SAD
#
mk1 <- pairKS_Dq(Dq2)
propNotDiffSAD(mk1) # 0.84 is not so good 

# Select 1 from each combination
#
Dqq <- sel_ReplaceR(Dq1,0.2,0.04,0.001,0)
Dqq$rep <- rep( 1:(nrow(Dqq)/35),each=35)
Diq1 <- Dqq[Dqq$rep==sample(10,1),]

Dqq <- sel_ReplaceR(Dq1,0.2,0.04,0.001,0.1)
Dqq$rep <- rep( 1:(nrow(Dqq)/35),each=35)
Diq1 <- rbind(Diq1,Dqq[Dqq$rep==sample(10,1),])

Dqq <- sel_ReplaceR(Dq1,0.2,0.04,0.001,0.01)
Dqq$rep <- rep( 1:(nrow(Dqq)/35),each=35)
Diq1 <- rbind(Diq1,Dqq[Dqq$rep==sample(10,1),])

Dqq <- sel_ReplaceR(Dq1,0.2,0.04,0.001,0.001)
Dqq$rep <- rep( 1:(nrow(Dqq)/35),each=35)
Diq1 <- rbind(Diq1,Dqq[Dqq$rep==sample(10,1),])

Dqq <- sel_ReplaceR(Dq1,0.2,0.04,0.001,1)
Dqq$rep <- rep( 1:(nrow(Dqq)/35),each=35)
Diq1 <- rbind(Diq1,Dqq[Dqq$rep==sample(10,1),])

mk1 <- pairKS_Dq(Diq1)
propNotDiffSAD(mk1) # 0.3 is very good 

# TO DO a function to compare one of each parameter combination.

hh <- function(x) {
n <- nrow(x)/35
rp <- rep( 1:n,each=35)
return(x[rp==sample(n,1),])
}
debug(hh)
Dq2 <- ddply(Dq1,.(MortalityRate,DispersalDistance,ColonizationRate,ReplacementRate),hh)

debug(pairKS_Dq)
mk1 <- pairKS_Dq(Dq2,35)
mk1 <- pairKS_Dq(Dq1,35)
propNotDiffSAD(mk1) # 0.3 is very good 
debug(compMethod_DqKS_Time)
m <- compMethod_DqKS_Time(bName,100,spMeta) 

