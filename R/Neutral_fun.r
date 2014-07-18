
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

# Calculate Ranks for every parameter combination, for ggplot2  
#
calcRankSAD  <- function(den)
{
  require(plyr)
  hh <- function(x) { 
  x1 <- x[x$value>0,]
  x1$parms <- paste(unique(x[,1:4]),collapse="_")
  x1$Rank <- nrow(x1) - rank(x1$value) +1
  return(x1)
  }
  ddply(den, .(MortalityRate,DispersalDistance,ColonizationRate,ReplacementRate),hh )
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

# Read simulation output and set variable names in wide format 
#
readWideDensityOut <- function(fname,num_sp){
  den <- read.delim(fname)
  names(den)[1:5]<-c("GrowthRate","MortalityRate","DispersalDistance","ColonizationRate","ReplacementRate")
  return(den)
}


# Proportion of not different SAD at 0.05 Hommel adjusted level 
#
propNotDiffSAD <- function(mk) nrow(mk[mk$p.adjust>0.05,c(10:13,14:15)])/nrow(mk)

# Proportion of not different SRS at 0.05 Hommel adjusted level 
#
propNotDiffSRS <- function(mk) nrow(mk[mk$adj.P.Value>0.05,])/nrow(mk)


# Calculates Dq from a data.frame read from the output of neutral model 
# auxiliar function for the next one
#
calcDq_frame <- function(pp)
{
  pp$Dq  <- with(pp,ifelse(q==1,alfa,Tau/(q-1)))
  pp$SD.Dq  <- with(pp,ifelse(q==1,SD.alfa,abs(SD.Tau/(q-1))))
  pp$R.Dq <- with(pp,ifelse(q==1,R.alfa,R.Tau))
  return(pp[,c(1:6,15:17)])
}              
# Reads the output of multifractal spectra of neutral model
# an calculates Dq 
#
readNeutral_calcDq <-function(fname)
{
  md1 <- read.table(fname,header=F,skip=1)
  md1 <- md1[,c(2:16)]
  names(md1)<-c("MortalityRate","DispersalDistance","ColonizationRate","ReplacementRate","Time","q","Tau","alfa","f(alfa)","R.Tau","R.alfa","R.f","SD.Tau","SD.alfa","SD.f")
  
  md1 <-calcDq_frame(md1)
}


# Read Dq output from neutral model and calculate pairwise differences 
#
# Dqf: data frame from readNeutral_calcDq 
# qNumber: number of q used
#
compDq_frame <- function(Dqf,qNumber)
{
  if( !require(statmod) & !require(reshape2))
    stop("required statmod and reshape2")
  
  # Subset for testing 
  #
  #Dqf <- with(Dqf,Dqf[MortalityRate==.2 & DispersalDistance==0.04 & ColonizationRate==0.001, ])
  
  # Set the number of repetitions using nrow and number of q  
  # 
  Dqf$rep <- rep( 1:(nrow(Dqf)/qNumber),each=qNumber)
  
  # Build variable for comparisons
  Dqf$factor <- do.call(paste, c(Dqf[,1:4],sep="_"))
    
  # Prepare data.frame in wide format 
  #
  Dq2 <- melt(Dqf, id.vars=c("q","rep","factor"),measure.var="Dq")
  Dq2 <- dcast(Dq2, factor+rep~ q)
  
  # Compare SRS curves
  #
  c2 <- compareGrowthCurves(Dq2$factor,Dq2[,3:37],nsim=1000)
}

# Plot Dq with fixed parameters except ReplacementRate
#
plotDq_ReplaceR <- function(Dqf,MortR,DispD,ColonR)
{
  require(plyr)
  c3 <- with(Dqf,Dqf[MortalityRate==MortR & DispersalDistance==DispD & ColonizationRate==ColonR, ])
  c3$factor <- do.call(paste, c(c3[,1:4],sep="_"))
  c3 <- ddply(c3, .(factor,q), summarize, SD.Dq=sd(Dq),Dq=mean(Dq))

  require(ggplot2)
  gp <- ggplot(c3, aes(x=q, y=Dq, colour=factor)) +
    geom_errorbar(aes(ymin=Dq-SD.Dq, ymax=Dq+SD.Dq), width=.1) +
    geom_point() + theme_bw()
  print(gp)

}

# Compare 2 communities
#
H_bss <- function (comm, R = 2000, index = "shannon") 
{
  if (length(rownames(comm)) > 2) 
    stop("Please supply community of only 2 samples.\n\n")
  require(vegan)
  RT <- rowSums(comm)
  DF <- sum(rowSums(comm > 0)) - 2
  CT <- colSums(comm)
  Hx <- as.numeric(diversity(comm, index = index))[1]
  Hy <- as.numeric(diversity(comm, index = index))[2]
  Ho <- abs(Hx - Hy)
  Hb <- vector(mode = "numeric", length = R)
  for (i in 1:R) {
    Mb <- r2dtable(1, RT, CT)
    Mb <- as.matrix(Mb[[1]])
    Hd <- as.numeric(abs(diversity(Mb, index = index)[1] - 
                           diversity(Mb, index = index)[2]))
    Hb[i] <- Hd
  }
  Q <- quantile(Hb, c(0.025, 0.975))
  Z <- abs((Ho - mean(Hb))/sd(Hb))
  P <- mean(abs(Hb) > Ho)
  Pz <- pt(Z, df = DF, lower.tail = FALSE) * 2
  result <- list(H.orig = c(Hx, Hy), statistic = Ho, z.score = Z, 
                 CI = Q, p.sim = P, p.z = Pz, df = DF, Index = index)
  class(result) = "Hbss"
  return(result)
}
