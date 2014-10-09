# Compare Dq without repetition the most interesting are 
# compareTwoCurvesT and KS.test


# Modified functions from package statmod to compare curves
#

# y1 is a matrix with columns mean,sd,n
#
meanAbsT <-function (y1, y2) 
{
  if (is.null(dim(y1)) || is.null(dim(y2))) 
    return(NA)
  y1 <- as.matrix(y1)
  y2 <- as.matrix(y2)
  if (ncol(y1) != ncol(y2)) 
    stop("Number of time points must match")
  m1 <- y1[,1]
  m2 <- y2[,1]
  v1 <- y1[,2]*y1[,2]
  v2 <- y2[,2]*y2[,2]
  n1 <- y1[,3]
  n2 <- y2[,3]
  s <- ((n1 - 1) * v1 + (n2 - 1) * v2)/(n1 + n2 - 2)
  t.stat <- abs(m1 - m2)/sqrt(s * (1/n1 + 1/n2))
  weighted.mean(t.stat, w = (n1 + n2 - 2)/(n1 + n2), na.rm = TRUE)
}

# Comparing two curves with SD for each point
# using group permutations and an average T statistics
#
# NOT working the curves are always equal!!!
#
compareTwoGCurves <- function (group, y, q, nsim = 100, fun = meanAbsT) 
{
    group <- as.vector(group)
    g <- unique(group)
    if (length(g) != 2) 
        stop("Must be exactly 2 groups")
    stat.obs <- fun(y[group == g[1], , drop = FALSE], y[group == 
        g[2], , drop = FALSE])
    asbig <- 0
    for (i in 1:nsim) {
        pgroup <- ave(group, q, FUN = sample)
        stat <- fun(y[pgroup == g[1], , drop = FALSE], y[pgroup == 
            g[2], , drop = FALSE])
        if (abs(stat) >= abs(stat.obs)) 
            asbig <- asbig + 1
    }
    list(stat = stat.obs, p.value = asbig/nsim)
}


# Comparing two curves with SD for each point  
# testing for (y1-y1)/(y1+y2)/2 have a normal distribution
# or that it is equal to 0 so they come from the same population
#
# NOT WORK for Dq curves very similar are always different
#
diffNormT <-function (y1, y2) 
{
  if (is.null(dim(y1)) || is.null(dim(y2))) 
    return(NA)
  y1 <- as.matrix(y1)
  y2 <- as.matrix(y2)
  if (ncol(y1) != ncol(y2)) 
    stop("Number of time points must match")
  m1 <- y1[,1]
  m2 <- y2[,1]
  v1 <- y1[,2]*y1[,2]
  v2 <- y2[,2]*y2[,2]
  n1 <- y1[,3]
  n2 <- y2[,3]
  
  s <- sqrt(v1/n1 + v2/n2)
  d <- (m1 - m2)/(s*(m1+m2)/2)
  qqnorm(d)
  print(shapiro.test(d))
  print(t.test((m1-m2),mu=0))
}

# Comparing two curves with SD for each point using monte-carlo 
# assuming each point have a normal distribution
# (Faster than average T )
#
# NOT WORKING The Dq curves are always equal!!
#
compareTwoCurvesMC <- function (group, y, sd,nsim = 1000) 
{
  vstat <- double(nsim)
  group <- as.vector(group)
  g <- unique(group)
  if (length(g) != 2) 
    stop("Must be exactly 2 groups")
  y1 <- y[group == g[1]]
  y2 <- y[group == g[2]]
  sd1 <-sd[group == g[1]]
  sd2 <-sd[group == g[2]]
  stat.obs <- sum((y1-y2)^2) # mean(abs(y1-y2))
  asbig <- 0
  l <- length(y1)
  for (i in 1:nsim) {
    z1=rnorm(l,y1,sd1)
    z2=rnorm(l,y2,sd2)
    
    stat <- sum((z1-z2)^2)
    vstat[i] <- stat
    if (abs(stat) >= abs(stat.obs)) 
      asbig <- asbig + 1
  }
  list(dist = vstat, stat = stat.obs, p.value = asbig/nsim)
}

# Comparing two curves with SD for each point with multiple T-tests 
# using Hommel procedure to adjust p
# 
# WORKS but difficult to say if the curves are different because it 
# test only the points.
# In general for Dq the range q=0.5-3 is significative
# and the rest is not. 
#
# If I use the fisher test to combine p it gives always significative
#
compareTwoCurvesT <- function (group, y, sd,n) 
{
#  group <- as.vector(group)
  g <- unique(group)
  if (length(g) != 2) 
    stop("Must be exactly 2 groups")
  y1 <- y[group == g[1]]
  y2 <- y[group == g[2]]
  sd1 <-sd[group == g[1]]
  sd2 <-sd[group == g[2]]
  sd1 <- sd1*sd1
  sd2 <- sd2*sd2
  stat <- (y1-y2)/sqrt((sd1/n)+(sd2/n))
  df <- (sd1/n+sd2/n)^2/( (sd1/n)^2 / (n-1) + (sd2/n)^2 / (n-1)) 
  #p.adjust <- p.adjust(pt(stat,df),method="hommel")
  p <- pt(stat,df)
  #print(ks.test(y1,y2))
  #print(t.test(y1-y2,mu=0))
  return(list("stat"=stat,"p"=p))
}

# Fisher procedure to combine p.value from multiple T test in one p.value
#
Fisher.test <- function(p) {
  Xsq <- -2*sum(log(p))
  p.val <- 1-pchisq(Xsq, df = 2*length(p))
  return(c(Xsq = Xsq, p.value = p.val))
}