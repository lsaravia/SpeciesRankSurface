
# Generate inp parameters file for simulation
#
genNeutralParms <- function(fname,side,mprob,birth,mortality,dispersal,colonization,replace=0){
  S <- length(mprob)
  parm <- data.frame( n=c(rep(0,3),1:S) ,prob=c(rep(0,3),prob=mprob))
  parm$text <- "" 
  parm$z <- 0
  parm$text[1] <- paste(side,side,sep="\t")
  parm$text[2] <- S
  parm$text[3] <- paste(0,birth,mortality,dispersal,colonization,replace,sep="\t")
  nff <- S+3
  parm$text[4:nff] <-   with(parm[4:nff,], paste(n,z,z,z,format(sort(prob),scientific=F),sep="\t"))                    
  if(!grepl(".inp",fname)) fname <-paste0(fname,".inp")

  write.table(parm$text,fname,sep="\t",row.names=F,col.names=F,quote=F)
}


# Generates pomac.lin paramenter file for multiple simulations
#
genPomacParms <- function(fname,GrowthR,MortR,DispD,ColonR,ReplaceR,numRep=1)
  {
  pom <- expand.grid(GrowthRate=GrowthR,MortalityRate=MortR, DispersalDistance=DispD, 
    ColonizationRate=ColonR,ReplacementRate=ReplaceR)
  pom <- pom[rep(seq_len(nrow(pom)), numRep), ]
  if(!grepl(".lin",fname)) fname <-paste0(fname,".lin")
  write.table(pom,fname,sep="\t",row.names=F,col.names=T,quote=F)
}

# Merge data from 2 simulations 
# 
# eks: 1 row data.frame with 2 sets of parms for which ks.test was done (fixed TIME)
# den1: full data in long format
#
mergePairSAD <- function(eks, den1){
  # eks <- mks[61,]
  
  den2 <-merge(den1,eks,by.x=c(1:4),by.y=c(2:5))[,1:6]
  den2 <- den2[den2$value>0,]
  den2$parms <-paste(eks[,2:5],collapse="_")
  den2$Rank <- nrow(den2) - rank(den2$value) +1
  
  den3 <-merge(den1,eks,by.x=c(1:4),by.y=c(6:9))[,1:6]
  den3 <- den3[den3$value>0,]
  den3$parms <-paste(eks[,6:9],collapse="_")
  den3$Rank <- nrow(den3) - rank(den3$value) +1
  den2 <- rbind(den2,den3)
  
}

# pairwise KS test for all combinations of parameters in denl= density in long format
# parameters are in columns 1:4
# 
#
pairKS_SAD <- function(denl){
  
  parms <- unique(denl[,1:4])
  combo <- combn(nrow(parms),2)
  
  require(plyr)
  mks <-adply(combo,2, function(x) {
    p1 <- parms[x[1],]
    p2 <- parms[x[2],]
    d1 <- merge(denl,p1)
    d2 <- merge(denl,p2)
    ks <- ks.test(d1$value,d2$value)
    out <-data.frame(p1,p2,t(paste(p1,p2,sep="_")),p.value=ks$p.value,stringsAsFactors=F)
    ln <-length(names(p1))*3
    names(out)[1:ln]<-c(paste0(abbreviate(names(p1)),1),paste0(abbreviate(names(p1)),2),names(p1))
    return(out)      
  })
  mks$p.adjust <- p.adjust(mks$p.value, method="hommel")
  return(mks)
}




# Read simulation output and change from wide to long format NO TIME
#
meltDensityOut_NT <- function(fname,num_sp){
  den <- read.delim(fname)
  names(den)[1:5]<-c("GrowthRate","MortalityRate","DispersalDistance","ColonizationRate","ReplacementRate")
  
  # from 7 to 473 there are species densities
  # Put simpler names to variables to identify species
  names(den)[7:(6+num_sp)]<-as.character(1:num_sp)
  
  require(reshape2)
  
  den1 <- melt(den,id.vars=c("MortalityRate","DispersalDistance","ColonizationRate","ReplacementRate"),measure.vars=c(7:473),variable.name="Species")
  
}
