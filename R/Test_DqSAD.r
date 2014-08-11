# Two problems 
# 
# Why DqSAD is not strictly decreasing
rm(list = ls())

load(".RData")

setwd("Simul")

Dq1 <- readNeutral_calcDq(paste0(bName,"T100mfSAD.txt"))
plotDq_ReplaceR(Dq1,0.2,0.04,0.001,"DqSAD")


Dqq <- sel_ReplaceR(Dq1,0.2,0.04,0.001,0)
Dqq <- sel_ReplaceR(Dq1,0.2,0.04,0.001,0.1)

Dqq[12:24, c("q","Dq","SD.Dq")]

# At q=6 it began to increase!!!

# Make simulations without pomac and use plotDqFit

par <- read.table("sim.par",quote="",stringsAsFactors=F)
# Change base name

# Number of time steps 
par[par$V1=="nEvals",]$V2 <- 100
# Change interval to measure Density and Diversity
par[par$V1=="inter",]$V2 <- 100
par[par$V1=="modType",]$V2 <- 4 # Hierarchical saturated
par[par$V1=="sa",]$V2 <- "S" # Save a snapshot of the model
par[par$V1=="baseName",]$V2 <- paste0("Exp",spMeta,"")
par[par$V1=="pomac",]$V2 <- 0 # 0:one set of parms 
# 1:several simulations with pomac.lin parameters 


write.table(par, "sim.par",sep="\t",row.names=F,col.names=F,quote=F)


# Run the model with the parameter ReplacementRate turned on = Hierarchical 
# and exponential dispersal 
#
# Set the location of the binary xp"

# make simulations

s <- paste(neuBin,"sim.par","fishE.inp")
s

if(simul) system(s,wait=F)


