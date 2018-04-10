# Realized Variance Example
library(highfrequency)
# Load Realized Variance Data
load("spxRV.RData")
# Clean and make a copy
spxRV2=spxRV[!is.na(spxRV)]
# Perform linear fit
spxHarFit= harModel(data=spxRV2 , periods = c(1,5,22), RVest = c("rCov"), 
             type="HARRV",h=1,transform=NULL)
summary(spxHarFit)
plot(spxHarFit)

# Linear fit check for HAR model
# Grab dates to avoid missing dates in RV data, use in independent variable selection
spxdates=index(spxRV2["2004"])
# Get fit, adjust index for fitted values to make into xts data type
rv2harFit=lm(spxRV2["2004"]~as.xts(spxHarFit$fitted.values)[as.Date(spxdates)] )
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

plot(vParkinson)  # optional, run interactively

# fit a linear model of Parkinson predictions to realized variance for 2004
# 1st, get dates from RV data because it has missing dates
spxdates=index(spxRV2["2004"])

# Now fit (remember parky is vol forecast and spxRV2 is variance)
rv2parkyFit=lm(spxRV2["2004"]~vParkinson[spxdates]^2 + 0)
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


# plot(spGFit) # optional, run interactively

# Build linear model of RV vs GARCH forecast
spxdates=index(spxRV2["2004"])
plot(spxRV2["2004"],main="RV vs GARCH Variance Forecast")
points(as.xts(sigma(spGFit))[as.Date(spxdates)]^2)
rv2garchFit=lm(spxRV2["2004"]~as.xts(sigma(spGFit))[as.Date(spxdates)]^2)
summary(rv2garchFit)
plot(rv2garchFit,which=1)
# plot forecast
plot(ugarchforecast(spGFit))

lvgSpec=ugarchspec(list(model="gjrGARCH"))
# Now call ugarch with levergage spec and data to get GARCH fit
spLvgFit=ugarchfit(spec=lvgSpec,data=spRets)
plot(spxRV2["2004"],main="RV vs GARCH w/Leverage Variance Forecast")
points(as.xts(sigma(spLvgGFit))[as.Date(spxdates)]^2)
rv2garchLvgFit=lm(spxRV2["2004"]~as.xts(sigma(spLvgGFit))[as.Date(spxdates)]^2)
summary(rv2garchLvgFit)
plot(rv2garchLvgFit,which=1)
# plot garch lvg forecast
plot(ugarchforecast(spLvgFit))


# Multi plot
plot(spxRV2["2004"],main="RV vs GARCH Variance Forecast")
points(as.xts(sigma(spGFit))[as.Date(spxdates)]^2,pch="X")

plot(spxRV2["2004"],main="RV vs HAR Variance Forecast")
points(as.xts(spxHarFit$fitted.values))

plot(spxRV2["2004"],main="RV vs Parky Variance Forecast")
points(as.xts((vParkinson[spxdates]^2)/252))  # Parkinson returns annualized vol +> trun into daily

# Models all together
plot(as.xts(sigma(spGFit))[as.Date(spxdates)]^2,main="GARCH vs HAR vs Parky")
lines(as.xts(spxHarFit$fitted.values)[spxdates])
lines(as.xts((vParkinson[spxdates]^2)/252))  # Parkinson returns annualized vol +> trun into daily



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




