
require(dplyr)

Dq1 <- group_by(DqT,Type,Side,NumSp,SAD) %>% summarize(count=n())

Dq1 <- filter(DqT, SAD=="Uniform",Side==256) %>% group_by(Type,Side,NumSp,SAD) %>% summarize(count=n())

Dq1 <- filter(DqT, SAD=="Uniform",Side==256,NumSp==64) %>% group_by(Type,Side,NumSp,SAD,rep) %>% summarize(count=n())

Dq1 <- filter(DqT, SAD=="Uniform",Side==256,NumSp==64) 
Dq3 <- distinct(Dq1,Type,Side,NumSp,SAD,rep,q)

Dq2 <- filter(DqT, !(SAD=="Uniform" & Side==256 & NumSp==64)) 
DqT1 <-rbind(Dq2,Dq3)  

Dq1 <- group_by(DqT1,Type,Side,NumSp,SAD) %>% summarize(count=n())

DqT_bak <- DqT
DqT <- DqT1
rm(Dq1,Dq2,Dq3,DqT1)
save.image()
