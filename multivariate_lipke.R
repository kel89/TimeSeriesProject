# multivariate_lipke.R

rm(list = ls())
setwd("~/Documents/STSCI 4550 -- Time Series/TimeSeriesProject")

# Questions:
# For co-integration, what do we do with seasonal differences?
# Such a thing as seasonal ADL model?
# The VAR model with GDP and CPI, what is that MQ test doing?


library(readxl)
library(CADFtest)
library(forecast)
library(vars)
library(stargazer)
library(MTS)

##############
#####Data#####
##############
GDP<- read_excel("Data_Roberto.xlsx",sheet="GDP",skip=17)
GDP<-ts(GDP[,2],start=c(1993,1),frequency=4)
GDP<-window(GDP,start=c(1997,1),end=c(2017,4))

CPI<-read_excel("Data_Roberto.xlsx",sheet ="CPI",skip=17)
CPI<-ts(CPI[,2],start=c(1969,1),frequency=12)
CPI<-aggregate(CPI,nfrequency=4)/3
CPI<-window(CPI,start=c(1997,1),end=c(2017,4))

FX<-read_excel("Data_Roberto.xlsx",sheet ="FX",skip=17)
FX<-ts(FX[,2],start=c(1968,1),frequency=12)
FX<-aggregate(FX,nfrequency=4)/3
FX<-window(FX,start=c(1997,1),end=c(2017,4))

MB<-read_excel("Data_Roberto.xlsx",sheet ="Monetary Base",skip=17)
MB<-MB[-c(1),]
MB<-ts(MB[,2],start=c(1986,1),frequency=12)
MB<-aggregate(MB,nfrequency=4)/3
MB<-window(MB,start=c(1997,1),end=c(2017,4))

UE<-read_excel("Data_Roberto.xlsx",sheet ="Unemployment",skip=10)
UE<-ts(UE[,2],start=c(1987,1),frequency=4)
UE<-window(UE,start=c(1997,1),end=c(2017,4))

RATES<-read_excel("Data_Roberto.xlsx",sheet ="Interest Rates",skip=12)
IRMEX90<-ts(RATES[,2],start=c(1997,1),frequency=4)
IRMEX90<-window(IRMEX90,start=c(1997,1),end=c(2017,4))
IRUSA90<-ts(RATES[,4],start=c(1997,1),frequency=4)
IRUSA90<-window(IRUSA90,start=c(1997,1),end=c(2017,4))


# Plot data before compression and window cut
par(mfrow=c(2,3))
ts.plot(GDP,main="GDP",ylab="Millions of Pesos")
ts.plot(CPI,main="CPI",ylab="Index")
ts.plot(FX,main="Foreign Exchange Rate",ylab="Pesos/Dollars")
ts.plot(MB,main="Money Supply",ylab="Millions of Pesos")
ts.plot(UE, main="Unemployment", ylab="Percentage")
ts.plot(IRMEX90,IRUSA90,col= c("blue","red"),main="3-Month Treasury Securities' Rates",ylab="Percentage")
legend("topright",legend=c("Mexican","American"),col= c("blue","red"),lty=1, cex=.75)
par(mfrow=c(1,1))



##### Lets explore Unemployment a bit more
plot.ts(UE)
lUE <- log(UE)
plot.ts(lUE)

# is it stationary?
q <- round(sqrt(length(lUE)))
CADFtest(lUE,type="trend",criterion="BIC",max.lag.y=q) # stochastic trend

# work in differences
ldUE = diff(lUE)
plot.ts(ldUE)
q <- round(sqrt(length(ldUE)))
CADFtest(ldUE,type="drift",criterion="BIC",max.lag.y=q) # stationary

# seasonal?
monthplot(ldUE)
sld_UE <- diff(ldUE, 4) # seasonal log differences
plot.ts(sld_UE)
q <- round(sqrt(length(sld_UE)))
CADFtest(ldUE,type="drift",criterion="BIC",max.lag.y=q) # stationary

# therefore, sld_UE is stationary, can we modeled with purely seasonal models


#### Dynamic Model 1: GDP ~ CPI ----
par(mfrow=c(2,1))
plot.ts(GDP)
plot.ts(CPI)
par(mfrow=c(1,1))

# take both in log diff and seasonal diff
sld_GDP <- diff(diff(log(GDP), 4))
sld_CPI <- diff(diff(log(CPI), 4))
par(mfrow=c(2,1))
plot.ts(sld_GDP)
plot.ts(sld_CPI)
par(mfrow=c(1,1))

# Try DL(1)
lag = 1
GDP_dat <- embed(sld_GDP, dimension = lag + 1)
CPI_dat <- embed(sld_CPI, dimension = lag + 1)
DL1 <- lm(GDP_dat[,1] ~ CPI_dat)
summary(DL1) # only significant at the 10% level
# stargazer(DL1) # nice latex output of regression table
par(mfrow=c(2,1))
acf(DL1$residuals)
pacf(DL1$residuals)
par(mfrow=c(1,1))
# These show that this model does not fit especially well


# Try DL(2) (though probably not good)
lag = 2
GDP_dat <- embed(sld_GDP, dimension = lag + 1)
CPI_dat <- embed(sld_CPI, dimension = lag + 1)
DL2 <- lm(GDP_dat[,1] ~ CPI_dat)
summary(DL1) # only significant at the 10% level
# stargazer(DL2) # nice latex output of regression table
par(mfrow=c(2,1))
acf(DL2$residuals)
pacf(DL2$residuals)
par(mfrow=c(1,1))
# this really is not helping, we should do a distributed lag model


#### Distributed Lag Model 1:l GDP ~ CPI ----

# looking at the pacf for the residuls from the DL(1) model, we can try just adding
# a single AR term of GDP to the model
lag = 1
GDP_dat <- embed(sld_GDP, dimension = lag + 1)
CPI_dat <- embed(sld_CPI, dimension = lag + 1)
ADL1 <- lm(GDP_dat[,1] ~ GDP_dat[,2] + CPI_dat[,2])
summary(ADL1)
par(mfrow=c(2,1))
acf(ADL1$residuals)
pacf(ADL1$residuals)
par(mfrow=c(1,1))
# here we get better significance

# Lets check for granger causality with this model 
small <- lm(GDP_dat[,1] ~ GDP_dat[,2])
anova(ADL1, small)
# This is showing that it is not good?



#### Dynamic Model 2: GDP ~ MB ----
par(mfrow=c(2,1))
plot.ts(GDP)
plot.ts(MB)
par(mfrow=c(1,1))

# take both in log diff and seasonal diff
sld_GDP <- diff(diff(log(GDP), 4))
sld_MB <- diff(diff(log(MB), 4))
par(mfrow=c(2,1))
plot.ts(sld_GDP)
plot.ts(sld_MB)
par(mfrow=c(1,1))

# Try DL(1)
lag = 1
GDP_dat <- embed(sld_GDP, dimension = lag + 1)
MB_dat <- embed(sld_MB, dimension = lag + 1)
DL1 <- lm(GDP_dat[,1] ~ MB_dat)
summary(DL1) # Nothing is significant
par(mfrow=c(2,1))
acf(DL1$residuals)
pacf(DL1$residuals)
par(mfrow=c(1,1))
# These show that this model does not fit especially well


# Try DL(2) (though probably not good)
lag = 2
GDP_dat <- embed(sld_GDP, dimension = lag + 1)
MB_dat <- embed(sld_MB, dimension = lag + 1)
DL2 <- lm(GDP_dat[,1] ~ MB_dat)
summary(DL2) # Nothing!
par(mfrow=c(2,1))
acf(DL2$residuals)
pacf(DL2$residuals)
par(mfrow=c(1,1))
# this really is not helping, we should do a distributed lag model


#### Distributed Lag Model 2:l GDP ~ MB ----
# looking at the pacf for the residuls from the DL(1) model, we can try just adding
# a single AR term of GDP to the model
lag = 1
GDP_dat <- embed(sld_GDP, dimension = lag + 1)
MB_dat <- embed(sld_MB, dimension = lag + 1)
ADL1 <- lm(GDP_dat[,1] ~ GDP_dat[,2] + MB_dat[,2])
summary(ADL1)
par(mfrow=c(2,1))
acf(ADL1$residuals)
pacf(ADL1$residuals)
par(mfrow=c(1,1))
# Cannot validate this model by acf and pacf

# Lets check for granger causality with this model 
small <- lm(GDP_dat[,1] ~ GDP_dat[,2])
anova(ADL1, small)
# Super no granger causality


#### Try a VAR model with GDP and CPI ----
data <- cbind(sld_GDP, sld_CPI)
VARselect(data)
# by BIC (SC) do a lag of 1

var1 <- vars::VAR(data, p=1)
summary(var1)
# from the looks of it, GDP can influence CPI, but not the other way aroudn
varresid <- resid(var1)
layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))
acf(varresid[,1])
acf(varresid[,2])
ccf(varresid[,1],varresid[,2])
# does not look good

# run the MQ test
mq(varresid, round(sqrt(length(sld_GDP)))) # umm what?

# irf
plot(irf(var1,ortho=FALSE,boot=TRUE))
# we see CPI does not have a later impact on GDP

