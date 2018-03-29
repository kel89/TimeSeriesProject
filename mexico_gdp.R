# mexico_gdp.r
# Kenneth Lipke
# Exploration and Initial modeling of Mexican GDP from 2013

setwd("~/Documents/STSCI 4550 -- Time Series/Project")
rm(list=ls())

# Read the data
df <- read.csv("gdp_2013.csv", header=T)
gdp <- ts(df[names(df)[2]], frequency=4, start=c(1993, 1)) # starts in Q1 of 1993

# plot
plot.ts(gdp, main="Mexican GDP from (in 2013 pesos)", ylab="GDP")

## For report we should get summary statistics
summary(gdp)

# Clearly not stationary (let determine the kind of trend)
library(CADFtest)
q <- round(sqrt(length(gdp))) # max lag
CADFtest(gdp, type="trend", max.lag.y=q, criterion = "BIC") # fail to reject, so stochastic

# Lets try in logs frst
log_gdp <- log(gdp)
CADFtest(log_gdp, type="trend", max.lag.y=q, criterion = "BIC") # deterministic --> work in differences
# we could do a trend and work with residuals, but not sure how we would handle predictions (I think it
# makes sense to work in differences)

# Take in differences
dlgdp <- diff(log_gdp)
plot.ts(dlgdp, main="Log GDP in First Differences", ylab="Not Sure?")

# determine if stationary
CADFtest(dlgdp, type="drift", max.lag.y=q, criterion = "BIC") # reject --> stationary (good)

# Look for seasonal trends
monthplot(dlgdp) # seems like it --> go in seasonal differences

# Seasonal differences
sdlgdp <- diff(dlgdp, 4)
plot.ts(sdlgdp)
CADFtest(sdlgdp, type="drift", max.lag.y=q, criterion = "BIC") # super reject --> very stationary (good)

# Start modeling
par(mfrow=c(2,1))
acf(sdlgdp)
pacf(sdlgdp)
par(mfrow=c(1,1))

# I would go with just seasonal AR term, and see what happens, I realize quite odd
library(forecast)
sar <- Arima(log_gdp, order=c(0,1,0), seasonal=c(1,1,0))
summary(sar)
acf(sar$residuals)  
Box.test(sar$residuals, lag=q, type="Ljung-Box")  # reject, not good enough

# I think we need another SAR term
sar2 <- Arima(log_gdp, order=c(0,1,0), seasonal=c(2,1,0))
summary(sar2)
acf(sar2$residuals)  
Box.test(sar2$residuals, lag=q, type="Ljung-Box") # this test passes

# We should probably try with an MA model so we have something to compare prediction accuracy to
  
  
  