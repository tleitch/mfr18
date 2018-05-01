# Realized Variance Example
library(highfrequency)
# Load Realized Variance Data
load("spxRV.RData")
# Clean and make a copy
spxRV2=spxRV[!is.na(spxRV)]
# Perform linear fit
spxHarFit= harModel(data=spxRV2 , periods = c(1,5,22), RVest = c("rCov"), 
             type="HARRV",h=1,transform=NULL)
plot(spxHarFit)

# Linear fit check for HAR model
# Grab dates to avoid missing dates in RV data, use in independent variable selection
spxdates=index(spxRV2["2004"])
# Get fit, adjust index for fitted values to make into xts data type
rv2harFit=lm(spxRV2["2004"]~as.xts(spxHarFit$fitted.values)[as.Date(spxdates)] +0)
summary(rv2harFit)




## Parkinson method example
library(quantmod)
# Load data
sp500=getSymbols('^GSPC',from="1999-01-01",warnings=F,verbose = F)
# Get Open High Low Close from quantmod
ohlc<-OHLC(GSPC)[,1:4]
# Change names to be compatible with TTR function "volatility"
names(ohlc)<-c("Open","High","Low","Close")
library(TTR)
# Fit data
vParkinson <- volatility(ohlc, n=1, calc="parkinson")

#plot(vParkinson)  # optional, run interactively

# fit a linear model of Parkinson predictions to realized variance for 2004
# 1st, get dates from RV data because it has missing dates
spxdates=index(spxRV2["2004"])

# Now fit (remember parky is vol forecast and spxRV2 is variance)
rv2parkyFit=lm(spxRV2["2004"]~vParkinson[spxdates]^2 +0)
summary(rv2parkyFit)

# Garch Fit Example
library(quantmod)
library(rugarch)
sp500=getSymbols('^GSPC',from="1999-01-01",warnings=F,verbose = F)
spRets=dailyReturn(GSPC,type='log',leading=TRUE)  # type=log gives lognormal returns
# First we need a model spec - see ugarchspec, default is standard GARCH
baseSpec=ugarchspec()
# Now call ugarch with spec and data to get GARCH fit
spGFit=ugarchfit(spec=baseSpec,data=spRets)
spGFit

# plot(spGFit) # optional, run interactively

# Build linear model of RV vs GARCH forecast
spxdates=index(spxRV2["2004"])
plot(spxRV2["2004"],main="RV vs GARCH Variance Forecast")
points(as.xts(sigma(spGFit))[as.Date(spxdates)]^2)
rv2garchFit=lm(spxRV2["2004"]~as.xts(sigma(spGFit))[as.Date(spxdates)]^2 +0)
summary(rv2garchFit)
plot(rv2garchFit)


# Multi plot
plot(spxRV2["2004"],main="RV vs GARCH Variance Forecast")
points(as.xts(sigma(spGFit))[as.Date(spxdates)]^2,pch="X")

plot(spxRV2["2004"],main="RV vs HAR Variance Forecast")
points(as.xts(spxHarFit$fitted.values))

plot(spxRV2["2004"],main="RV vs Parky Variance Forecast")
points(as.xts((vParkinson[spxdates]^2)/252))  # Parkinson returns annualized vol +> trun into daily

# Models all together
plot(as.xts((vParkinson[spxdates]^2)/252),main="GARCH vs HAR vs Parky")
lines(as.xts(spxHarFit$fitted.values),col="green")
lines(as.xts(sigma(spGFit))[as.Date(spxdates)]^2,col="red")  # Parkinson returns annualized vol +> trun into daily



## Now look at PNL resulting from limit strategies  ##

## This function calculates the position for $100k of risk 
## at 99th %-ile given a vol forecast 
## and the realized returns

cumPNL  <-function(rets,vols,conf=.01,lag=0){
    # get forecast VaR
    pnl=cumPNL=position=vector(length(rets)-lag,mode="numeric")
    for(i in (1+lag):(length(rets))){
        VaR=as.numeric(-vols[i]*qnorm(conf))
        position[i-lag]=100000/VaR
        pnl[i-lag]=position[i-lag]*rets[i-lag]
        if(i>2){
            cumPNL[i-lag]=cumPNL[i-1-lag]+pnl[i-lag]
        } else {
            cumPNL[i-lag]=pnl[i-lag]
        }
    }
    
    return(xts(data.frame(pnl,cumPNL,position),order.by=index(rets)[(1+lag):length(rets)]))
    
}


## Cacluate position P&L for strategies for time period 2004
parkyPNL=cumPNL(spRets[as.Date(spxdates)],vParkinson[as.Date(spxdates)]/sqrt(252))
garchPNL=cumPNL(spRets["2004"],sigma(spGFit)["2004"])
harPNL=cumPNL(spRets[as.Date(spxdates)],sqrt(as.xts(spxHarFit$fitted.values)[as.Date(spxdates)]))

## plot Cummulative P&L for position strategy
plot(parkyPNL$cumPNL,ylim=c(-500000,1000000),main="Cummulative P&L for Limit Strategies")
lines(garchPNL$cumPNL,col="red")
lines(harPNL$cumPNL,col="green")

## Now look at position forecast 
plot(parkyPNL$position,main="Trade Position for Limit Strategies")
lines(garchPNL$position,col="red")
lines(harPNL$position,col="green")





### HW 2 Examples Start here  #########################################################################

### RiskMetric   Utility Function #######
rmVol <- function(X_t){
    N=length(X_t)
    require(timeSeries)
    sigma_t <- c()
    is.timeSeries(sigma_t)
    w<- 0.94^(0:74) #lambda = 0.94, T=77 as specified in RiskMetrics
    sigma_t[77]<- sum(w*X_t[76:2]^2)/sum(w) #compute the 77th term 
    for (s in 78:N){ #and all the others
        sigma_t[s]<- 0.94*sigma_t[s-1]+0.06*X_t[s-1]^2
    }
    sigma_t<- sqrt(sigma_t) #here it is the conditioned volatility serie
    sigma_t=as.xts(sigma_t,order.by=index(X_t)) 
    sigma_t
    
}


## RM Aanslysis #################################
rmVolatilities=rmVol(spRets)
spxdates=index(spxRV2["2004"])
plot(spxRV2["2004"],main="RV vs RiskMetrics Variance Forecast")
points(rmVolatilities[as.Date(spxdates)]^2)
rv2riskmetricsFit=lm(spxRV2["2004"]~rmVolatilities[as.Date(spxdates)]^2 + 0)
summary(rv2riskmetricsFit)
plot(rv2riskmetricsFit)

garchPNL=cumPNL(spRets["2004"],rmVolatilities["2004"])
# Models all together
## plot Cummulative P&L for position strategy
plot(parkyPNL$cumPNL,ylim=c(-500000,1000000),main="Cummulative P&L for Limit Strategies")
lines(garchPNL$cumPNL,col="red")
lines(harPNL$cumPNL,col="green")
lines(garchPNL$cumPNL,col="blue")

## Now look at position forecast 
plot(parkyPNL$position,main="Trade Position for Limit Strategies")
lines(garchPNL$position,col="red")
lines(harPNL$position,col="green")
lines(garchPNL$position,col="blue")





### Conditional Coverage and Independence  ##################################

## utility functions to calculate LRuc and LRind, the likelihood ratios for unconditional coverage and independence
exceptions = function(VaR, returns){
    #Generates a vector with 0's in the days in which exceptions do not occur
    #and the daily return in the ones where they do
    #VaR = Value-at-Risk estimates for n days
    #returns = daily returns data of the n days
    
    output = vector(mode = "numeric", length = length (VaR))
    
    for (i in seq(length(VaR))){
        if (-returns[i] > VaR[i]){
            output[i] = returns[i]
        }
    }
    return(output)
}

LRuc = function(excepts, alpha){
    #Generates Kupiec Unconditional Coverage Test statistic - equation 2.28
    #excepts = vector with 0's in the days in which exceptions do not occur
    #and the daily return in the ones where they do
    #alpha = alpha used in the VaR estimation
    
    n = length(excepts)
    m = sum(excepts != 0)
    stat = 2 * log  ((1 - m/n)^(n - m) * (m/n)^m) - 2 * log((1 - alpha)^(n - m) * alpha^m)  
    
    return (stat)
}


LRind = function(excepts){
    #Generates Christoffersen Serial Independence Test statistic - equation 2.29
    #excepts = vector with 0's in the days in which exceptions do not occur
    #and the daily return in the ones where they do
    
    n00 = 0
    n01 = 0
    n10 = 0
    n11 = 0
    
    for (i in seq(length(excepts) - 1)){
        
        if (excepts[i] == 0){
            if (excepts[i+1] == 0) n00 = n00 + 1
            
            else n01 = n01 + 1
        }
        
        if (excepts[i] != 0){
            if (excepts[i+1] == 0) n10 = n10 + 1
            
            else n11 = n11 + 1
        }
    }
    
    pi01 = n01 / (n00 +n01)
    pi11 = n11 / (n10 +n11)
    pi = (n01 + n11) / (n00 +n01 + n10 +n11)
    LR = 2 * log  ((1 - pi01)^(n00) * pi01^n01 * (1 - pi11)^(n10) * pi11^n11) - 2 * log((1 - pi)^(n00 + n10) * pi^(n01 + n11))
    
    return (LR)
}




### Unconditional Coverage  (event count) ####################################
# get vector of VaR estimates at 95th %-ile for standard Garch
varEstGch95=-sigma(spGFit)["2004"]*qnorm(.05)
gch95Exceptions=exceptions(varEstGch95,spRets["2004"])
gchUC95=LRuc(gch95Exceptions,.05)   ### reject at LRuc >2.7


## Now 99th
varEstGch99=-sigma(spGFit)["2004"]*qnorm(.01)
gch99Exceptions=exceptions(varEstGch99,spRets["2004"])
gchUC99=LRuc(gch99Exceptions,.01)


# get vector of VaR estimates at 95th %-ile for Parky
varEstPrky95=-vParkinson["2004"]*qnorm(.05)/sqrt(252)
prky95Exceptions=exceptions(varEstPrky95,spRets["2004"])
prkyUC95=LRuc(prky95Exceptions,.05)


varEstPrky99=-vParkinson["2004"]*qnorm(.01)/sqrt(252)
prky99Exceptions=exceptions(varEstPrky99,spRets["2004"])
prkyUC99=LRuc(prky99Exceptions,.01)

# get vector of VaR estimates at 95th %-ile for Parky
spxdates=index(spxRV2["2004"])
varEstHar95=-sqrt(as.xts(spxHarFit$fitted.values)[as.Date(spxdates)])*qnorm(.05)
har95Exceptions=exceptions(as.vector(varEstHar95),as.vector(spRets[as.Date(spxdates)]))
harUC95=LRuc(har95Exceptions,.05)

varEstHar99=-sqrt(as.xts(spxHarFit$fitted.values)[as.Date(spxdates)])*qnorm(.01)
har99Exceptions=exceptions(as.vector(varEstHar99),as.vector(spRets[as.Date(spxdates)]))
harUC99=LRuc(har99Exceptions,.01)

### Independence (events not correlated) ####
# garch
gchIND95=LRind(gch95Exceptions)
pchisq(gchIND95,df=1)
gchIND99=LRind(gch99Exceptions)   ### reject if LRind>2.7

# parky
prkyIND95=LRind(prky95Exceptions)
pchisq(prkyIND95,df=1)
prkyIND99=LRind(prky99Exceptions)

# har
harIND95=LRind(har95Exceptions)
pchisq(harIND95,df=1)
harIND99=LRind(har99Exceptions)

# riskmetrics
rmIND95=LRind(rm95Exceptions)
pchisq(rmIND95,df=1)
rmIND99=LRind(rm99Exceptions)

### Conditional Coverage = LRind + LRuc  and again reject if LRcc >2.7
## Do it yourself but only on chosen model
