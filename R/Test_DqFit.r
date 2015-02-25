#
# Testing the fit of Dq applying Liebovitch 1989 oscillations
#

source("R/Neutral_fun.r")

setwd("Simul")

zq <- readDq_fit(256,64,"Logseries")
g <-plotDqFitG(zq)
g

# Calculate small sample corrected  AIC
AICc <- function(mo)
{
  k <-attributes(logLik(mo))$df-1
  n <-attributes(logLik(mo))$nobs
  AIC(mo)+2*k*(k+1)/(n-k-1)
}

# Calculate AICc differences between linear an linear+quadratic models
#
linQuadAic <- function(y,x){
  mod1 <- lm(y~x)
  mod2 <- lm(y~x+I(x^2))
  AICc(mod1)-AICc(mod2)
}


# Calculate AICc differences between linear an linear+periodic models (tiene problemas)
#
linPerAic <- function(y,x){
  mod1 <- lm(y~x)
  mod3 <- lm(y~x+I(sin(2*pi*x)),data=zq1)
  
  c1 <- round(coef(mod3)[1],4)
  c2 <- round(coef(mod3)[2],4)
  c3 <- round(coef(mod3)[3],4)
  mod3 <- nls(y~a + b*x+ c*sin(2*pi*x/d),
              data=zq1,
              start=list(a=c1,b=c2,c=c3,d=1))

  AICc(mod1)-AICc(mod3) # r2_l = mod1$r.squared 
}


# Compare linear and quadratic models using AICc
#
require(dplyr)
aDif <- group_by(zq,Type,q) %>% summarize(AICDif=linQuadAic(logTr,LogBox))

ggplot(aDif,aes(q,AICDif,colour=Type)) + geom_point() + geom_line()


ggplot(aDif,aes(q,AICDif,shape=Type)) + geom_point()


# Select 
zq1 <- filter(zq,q==-1,Type=="DqSAD")

# Test quadratic term
#
mod1 <- lm(logTr~LogBox,data=zq1)
summary(mod1)
coef(mod1)
mod2 <- lm(logTr~LogBox+I(LogBox^2),data=zq1)
summary(mod2)
coef(mod2)[3]
mod3 <- lm(logTr~LogBox+I(sin(2*pi*LogBox)),data=zq1)

c1 <- round(coef(mod3)[1],4)
c2 <- round(coef(mod3)[2],4)
c3 <- round(coef(mod3)[3],4)
mod3 <- nls(logTr~a + b*LogBox+ c*sin(2*pi*LogBox/d),
            data=zq1,
            start=list(a=c1,b=c2,c=c3,d=1))
summary(mod3)
coef(mod3)
AICc(mod1)
AICc(mod2) 
AICc(mod3) 

coef(mod1)
coef(mod3)

x <- seq(min(zq1$LogBox),max(zq1$LogBox),0.01)
plot(zq1$LogBox,zq1$logTr,data=zq1)
lines(x,predict(mod1,list(LogBox=x)),col="brown")
lines(x,predict(mod2,list(LogBox=x)),col="red")
lines(x,predict(mod3,list(LogBox=x)),col="blue")

anova(mod1,mod2)
AIC(mod1,mod2)
AICc(mod1)-AICc(mod2)
AICc(mod1)-AICc(mod3)

###### FUNCTION TO CALCULATE DIF IN AICc for all the models for use with dplyr/plyr


# Select 
zq1 <- filter(zq,q==-1,Type=="DqSRS")

# Test quadratic term
#
mod1 <- lm(logTr~LogBox,data=zq1)
summary(mod1)
mod2 <- lm(logTr~LogBox+I(LogBox^2),data=zq1)
summary(mod2)
anova(mod1,mod2)
AIC(mod1,mod2)
AICc(mod1)-AICc(mod2)


# Select 
zq1 <- filter(zq,q==-1,Type=="rnzSRS")

# Test quadratic term
#
mod1 <- lm(logTr~LogBox,data=zq1)
summary(mod1)
mod2 <- lm(logTr~LogBox+I(LogBox^2),data=zq1)
summary(mod2)
AIC(mod1,mod2)
AICc(mod1)-AICc(mod2)

# Nonlinear estimation has many problems! 

oscill_frac <- function(x){
  (a * x ^b)*sin(2*pi*log(x)/log(c)) 
}

zq1$Tr <- exp(zq1$logTr)
nm1 <-nls(Tr ~ oscill_frac(BoxSize),data=zq1,start = list(a = 1, b = 1, c = 10),
          algorithm = "port")

