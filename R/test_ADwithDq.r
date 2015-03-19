# Multifractal analysis of multispecies spatial distributions - Anderson-Darling test

#I use Dq from the simple spatial patterns with uniform and logseries SAD (DqT)


load(".RData")

simul  <- F # variable to perform or not the simulations

oldcd <-getwd()
source("R/Neutral_fun.r")

# Set the location of the binary 
#
mfBin <- "~/Dropbox/cpp/SpatialAnalysis/mfsba" # /mfSBA /multiSpeciesSBA

require(pander)
panderOptions('table.split.table',Inf)
panderOptions("table.style", "multiline") 
options("scipen"=100, "digits"=4)
knitr::opts_chunk$set(comment = NA)

# Define the function


# test the permutation with Anderson-Darling statistic for Dq
# using the global variable DqT
#
test_ADPerm <-function(type="DqSRS",nsp=8,nsp1=64,maxq=24){

  require(ggplot2)
  require(kSamples)
  require(plyr)
  require(dplyr)

  # select 1 set for comparison
  # 
  Dq1<- filter(DqT,Side==256,NumSp==nsp,SAD=="Logseries",Type==type,abs(q)<=maxq) %>%
        group_by(SAD,Type,q) %>% summarize(SD.Dq=sd(Dq),Dq=mean(Dq),count=n(),Curve=paste(type,nsp))%>%
        select(Curve,Dq,SD.Dq,q)
        
  # Change location add a constant

  Dq2 <- mutate(Dq1,Dq=Dq+mean(Dq)*0.01,Curve=paste(type,nsp,"+ 0.1"))


  g <- ggplot(rbind(Dq1,Dq2),aes(x=q,y=Dq,colour=Curve)) + theme_bw()+ geom_point() +  ggtitle("Add a constant")
  print(g)  
  g <- ggplot(rbind(Dq1,Dq2),aes(x=Dq,colour=Curve))+theme_bw()+stat_ecdf()
  print(g)


  ad <- ad.test(Dq1$Dq,Dq2$Dq,method="simul")
  print(ad)


  # Delete rows with uneven q

  Dq2 <- filter(Dq1,(q%%2)==0) %>% mutate(Curve=paste(type,nsp," even"))

  g <- ggplot(rbind(Dq1,Dq2),aes(x=q,y=Dq,colour=Curve)) + theme_bw()+ geom_point() + ggtitle("Only even q")
  print(g)
  g <- ggplot(rbind(Dq1,Dq2),aes(x=Dq,colour=Curve))+theme_bw()+stat_ecdf()
  print(g)

  ad <- ad.test(Dq1$Dq,Dq2$Dq,method="simul")
  print(ad)



  ## small differences in shape but same endings
  #
  # first stretch

  p1 <- 3
  p2 <- nrow(Dq1)-2
  q1 <- Dq1$q[p1]
  q2 <- Dq1$q[p2]
  Dq2 <- filter(Dq1,(q<q1 | q>q2) | q==0)

  sp <-splinefun(Dq2$q,Dq2$Dq)
  Dq2 <- mutate(Dq1,Dq=sp(q),Curve=paste(type,nsp," spline")) #%>% filter(abs
  g <- ggplot(rbind(Dq1,Dq2),aes(x=q,y=Dq,colour=Curve)) + theme_bw()+ geom_point() + ggtitle("Spline with similar endings")
  print(g)
  g <- ggplot(rbind(Dq1,Dq2),aes(x=Dq,colour=Curve))+theme_bw()+stat_ecdf()
  print(g)

  ad <- ad.test(Dq1$Dq,Dq2$Dq,method="simul")
  print(ad)

  ## linear shape with same ending
  #
  p1 <- 4
  p2 <- nrow(Dq1)-3
  q1 <- Dq1$q[p1]
  q2 <- Dq1$q[p2]
  D1 <- Dq1$Dq[p1]
  D2 <- Dq1$Dq[p2]

  y <- c(Dq1$Dq[1:(p1-1)],D1 + (D2-D1)/(q2-q1)*(Dq1$q[p1:p2]-q1),Dq1$Dq[(p2+1):nrow(Dq1)])

  Dq2 <- mutate(Dq1,Dq=y,Curve=paste(type,nsp,"linear"))

  g <- ggplot(rbind(Dq1,Dq2),aes(x=q,y=Dq,colour=Curve)) + theme_bw()+ geom_point() + ggtitle("Linear with similar endings")
  print(g)
  g <- ggplot(rbind(Dq1,Dq2),aes(x=Dq,colour=Curve))+theme_bw()+stat_ecdf()
  print(g)

  ad <- ad.test(Dq1$Dq,Dq2$Dq,method="simul")
  print(ad)


  ## Diffferent no. species
  #

  Dq2<- filter(DqT,Side==256,NumSp==nsp1,SAD=="Logseries",Type==type,abs(q)<=maxq) %>%
        group_by(SAD,Type,q) %>% summarize(SD.Dq=sd(Dq),Dq=mean(Dq),count=n(),Curve=paste(type,nsp1))%>%
        select(Curve,Dq,SD.Dq,q)
        
  g <- ggplot(rbind(Dq1,Dq2),aes(x=q,y=Dq,colour=Curve)) + theme_bw()+ geom_point() + ggtitle("Diffferent no. species")
  print(g)
  g <- ggplot(rbind(Dq1,Dq2),aes(x=Dq,colour=Curve))+theme_bw()+stat_ecdf()
  print(g)

  ad <- ad.test(Dq1$Dq,Dq2$Dq,method="simul")
  print(ad)

  #
  #
  type1 <- paste0("rnz",type)
  Dq2<- filter(DqT,Side==256,NumSp==nsp,SAD=="Logseries",Type==type1,abs(q)<=maxq) %>%
        group_by(SAD,Type,q) %>% summarize(SD.Dq=sd(Dq),Dq=mean(Dq),count=n(),Curve=paste(type,nsp," rnz")) %>%
        select(Curve,Dq,SD.Dq,q)
  g <- ggplot(rbind(Dq1,Dq2),aes(x=q,y=Dq,colour=Curve)) + theme_bw()+ geom_point() + ggtitle("Randomized spatial pattern")
  print(g)
  g <- ggplot(rbind(Dq1,Dq2),aes(x=Dq,colour=Curve))+theme_bw()+stat_ecdf()
  print(g)

  ad <- ad.test(Dq1$Dq,Dq2$Dq,method="simul")
  print(ad)


  #
  #

  Dq2<- filter(DqT,Side==256,NumSp==nsp,SAD=="Logseries",Type==type1,abs(q)<=maxq) %>%
        group_by(SAD,Type,q) %>% summarize(SD.Dq=sd(Dq),Dq=mean(Dq),count=n(),Curve=paste(unique(Type),unique(NumSp)," sample"))%>%
        select(Curve,Dq,SD.Dq,q)

  Dq2 <- mutate(Dq2,Dq=rnorm(nrow(Dq2),mean=Dq,sd=SD.Dq))

  g <- ggplot(rbind(Dq1,Dq2),aes(x=q,y=Dq,colour=Curve)) + theme_bw()+ geom_point() + ggtitle("Sample from different Dq")
  print(g)
  g <- ggplot(rbind(Dq1,Dq2),aes(x=Dq,colour=Curve))+theme_bw()+stat_ecdf()
  print(g)

  ad <- ad.test(Dq1$Dq,Dq2$Dq,method="simul")
  print(ad)

  #
  #

  Dq2 <- mutate(Dq1,Dq=rnorm(nrow(Dq1),mean=Dq,sd=SD.Dq),Curve=paste(type,nsp," sample"))

  g <- ggplot(rbind(Dq1,Dq2),aes(x=q,y=Dq,colour=Curve)) + theme_bw()+ geom_point() + ggtitle("Sample from the same Dq")
  print(g)
  g <- ggplot(rbind(Dq1,Dq2),aes(x=Dq,colour=Curve))+theme_bw()+stat_ecdf()
  print(g)

  ad <- ad.test(Dq1$Dq,Dq2$Dq,method="simul")
  print(ad)

}



# Test Anderson with 8 sp

test_ADPerm("DqSRS",8,256)



# Test Anderson with 64 sp

test_ADPerm("DqSRS",64,256)


# Test Anderson with 64 sp 

test_ADPerm("DqSAD",64,256)



