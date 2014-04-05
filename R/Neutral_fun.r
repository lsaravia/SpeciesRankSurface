
# Generate inp parameters file for simulation
#
genNeutralParms <- function(fname,side,mprob,birth,mortality,dispersal,colonization){
  S <- length(mprob)
  parm <- data.frame( n=c(rep(0,3),1:S) ,prob=c(rep(0,3),prob=mprob))
  parm$text <- "" 
  parm$z <- 0
  parm$text[1] <- paste(side,side,sep="\t")
  parm$text[2] <- S
  parm$text[3] <- paste(0,birth,mortality,dispersal,colonization,sep="\t")
  nff <- S+3
  parm$text[4:nff] <-   with(parm[4:nff,], paste(n,z,z,z,format(sort(prob),scientific=F),sep="\t"))                    
  fname <- paste0(fname,S,".inp")
  write.table(parm$text,fname,sep="\t",row.names=F,col.names=F,quote=F)
}


# Generates pomac.lin paramenter file for multiple simulations
#
genPomacParms <- function(fname,GrowthR,MortR,DispD,ColonR,ReplaceR,numRep=1)
  {
  pom <- expand.grid(GrowthRate=GrowthR,MortalityRate=MortR, DispersalDistance=DispD, 
    ColonizationRate=ColonR,ReplacementRate=ReplaceR)
  pom <- pom[rep(seq_len(nrow(pom)), numRep), ]
  if(!grepl("lin",fname)) fname <-paste0(fname,".lin")
  write.table(pom,fname,sep="\t",row.names=F,col.names=T,quote=F)
}