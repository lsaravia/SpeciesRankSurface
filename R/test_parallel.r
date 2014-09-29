require(doParallel)
cl <- makeCluster(2)
registerDoParallel(cl)
require(foreach)

cc <-foreach(i=1:nrow(simul)) %dopar% 
{
  p <- data.frame(a1=i,a=runif(5),b=runif(5))
  c <- data.frame(c1=i,c=runif(5),d=runif(5))
  cc <- list("p"=p, "c"=c)
}
po <- data.frame()
co <- data.frame()
for(i in 1:nrow(simul)) {
  po <- rbind(po, cc[[i]]$p)
  co<- rbind(co,cc[[i]]$c)
}

library(rbenchmark)
benchmark( foreach(i=1:nrow(simul)) %dopar% 
{
  p <- data.frame(a1=i,a=runif(5000),b=runif(5000))
  c <- data.frame(c1=i,c=runif(5000),d=runif(5000))
  cc <- list("p"=p, "c"=c)
}, 
foreach(i=1:nrow(simul)) %do% 
{
  p <- data.frame(a1=i,a=runif(5000),b=runif(5000))
  c <- data.frame(c1=i,c=runif(5000),d=runif(5000))
  cc <- list("p"=p, "c"=c)
})
