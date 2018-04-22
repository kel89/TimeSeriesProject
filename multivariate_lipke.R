# multivariate_lipke.R

rm(list = ls())
setwd("~/Documents/STSCI 4550 -- Time Series/TimeSeriesProject")

# Questions:
# For co-integration, what do we do with seasonal differences?
# Such a thing as seasonal ADL model?
# The VAR model with GDP and CPI, what is that MQ test doing?
#
# Stuggling to interpret MQ test for all 3

# For gdp ~ ue, when I am doing ADL model, how can I get it to fit better, with like
# seasonal terms? I can't valide the model for Granger causaulity without it. Can I just add my 
# own by cherry picking correctly from the `embed` matrix?

# How can we do a seasonal VAR model, which is what I think it wants

# With cointegration, if the series are in seasonal and normal differences
# in the EG test, do we want to run the regression with the series in tottal levels
# or should we leave the seasonal differences?

# What if our series are NOT cointegrated, do we just say that and move on?
# Then we can't do an ECM? Then we need a GARCH model, and we haven't gotten that
# far, but I am not sure if it will be Garch... then what?





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


#### GDP ~ CPI + MB ---- 
par(mfrow = c(3,1))
ts.plot(GDP)
ts.plot(CPI)
ts.plot(MB)

ts.plot(sld_GDP)
ts.plot(sld_CPI)
ts.plot(sld_MB)

# all 3 stationary so lets run normal regression
d <- lm(sld_GDP ~ sld_CPI + sld_MB)
summary(d)
par(mfrow=c(2,1))
acf(d$residuals)
pacf(d$residuals)

l_GDP <- log(GDP)
l_CPI <- log(CPI)
l_MB <- log(MB)
# TRY SAR(1)(1) for erros
arima_d <- Arima(sld_GDP, order=c(1,0,0), seasonal=c(1,0,0), xreg=cbind(sld_CPI, sld_MB))
summary(arima_d) # this shows serious significance!
acf(arima_d$residuals) # lookin good! 
Box.test(arima_d$residuals, lag = round(sqrt(length(arima_d$residuals))), type = "Ljung-Box")
# Super white noise, so that is a good descriptive model
# Looking again at the summary, we see positive in CPI, and negative in MB
# makes sense for CPI, and as MB goes up, this devalues currency, and if the price stays the same
# which is the C.P. assumption, then people are less likely to buy mexican goods abroad?
# which would drive down GDP? I am not sure if that makes sense?


# Now lets look at cross correlations
par(mfrow=c(2,1))
ccf(x=sld_CPI[,1], y=sld_GDP[,1])
ccf(x=sld_MB[,1], y=sld_GDP[,1])
par(mfrow=c(1,1))
# This would lead me to believe that there is no good ADL model to specify?

### NOT SURE WHAT TO DO ABOUT ADL AND GRANGER IF I CAN'T SPEICFY WITH CCF????????




#### VAR model
mydata <- cbind(sld_GDP, sld_CPI, sld_MB)
VARselect(mydata)
# BIC says 1, AIC says 9 (way too many), do 1
var_mod <- vars::VAR(mydata, p=1)
summary(var_mod)
varresid <- resid(var_mod)
par(mfrow=c(3,2)) # there are 7
acf(varresid[,1])
acf(varresid[,2])
acf(varresid[,3])
ccf(varresid[,1], varresid[,2])
ccf(varresid[,1], varresid[,3])
par(mfrow=c(1,1))
# Run a test, because I am honestly not sure?
mq(varresid, lag=floor(sqrt(length(sld_GDP))))
# what do we make of this

# irf
plot(irf(var_mod,ortho=FALSE,boot=TRUE))
# literally noting exciting at all







# What if we did GDP and unemployment?
# is unemployment even stationary?
ts.plot(diff(UE))
monthplot(diff(UE))
ts.plot(diff(diff(UE), 4))
monthplot(diff(diff(UE), 4))
sd_UE <- diff(diff(UE), 4)
q <- round(sqrt(length(sd_UE)))
CADFtest(sd_UE,type="drift",criterion="BIC",max.lag.y=q) # stationary

# Can try dynaim model
gdp_ue <- lm(sld_GDP ~ sd_UE)
summary(gdp_ue) # so in differences it is significant

# Dynamic model
ccf(y = sld_GDP[,1], x = sd_UE[,1]) # there is some significance at one lag, and it appears at a seasonal
# You messed up, it looks like GDP is predicting unemployment... we could spin that
#DL(1)
lag <- 1
gdp.e <- embed(sld_GDP, dimension=lag+1)
ue.e <- embed(sld_UE, dimension = lag+1)
dl1 <- lm(gdp.e[,1] ~ ue.e)
summary(dl1) # marginal significance at the lag
acf(dl1$residuals) # def have some structure in the residuals
Box.test(dl1$residuals, lag = round(sqrt(length(dl1$residuals))), type = "Ljung-Box") # not white noise

# ADL model (I need some seasonal terms)
lag <- 4
gdp.e <- embed(sld_GDP, dimension=lag+1)
ue.e <- embed(sld_UE, dimension = lag+1)
adl1 <- lm(gdp.e[,1] ~ gdp.e[,-1] + ue.e) # no predicitive ability
summary(adl1)
acf(adl1$residuals) # still not good
# take a look at an ADL with predicitive ability
adlp <- lm(gdp.e[,1] ~ gdp.e[,-1] + ue.e[,-1])
summary(adlp) # when controlling for past GDP, UE does not much matter
acf(adlp$residuals) # still not good?

# ADL seasonal model (my attempt) # lets do 1 seasonal terms
# using lag 4 data
adl_s <- lm(gdp.e[,1] ~ gdp.e[,2] + gdp.e[,5] + ue.e[,2] + ue.e[,5]) # lag 1 and 1 seasonal
summary(adl_s)
acf(adl_s$residuals)
Box.test(adl_s$residuals, lag = round(sqrt(length(adl_s$residuals))), type = "Ljung-Box") # Super white noise

# Try granger causality with my seasonal model
small <- lm(gdp.e[,1] ~ gdp.e[,2] + gdp.e[,5])
anova(small, adl_s) # we see marginal significance (at the 8% level) *worth reporting

# Lets try a VAR model
dat <- cbind(sld_GDP, sd_UE)
VARselect(dat) # 4 lags for BIC (i.e. seasonal is what it wants)
var_mod <- vars::VAR(dat, p=4)
summary(var_mod) # significance at ue_{t-4} and gdp_{t-1}, interesting
res <- resid(var_mod)
layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))
acf(res[,1])
acf(res[,2])
ccf(res[,1],res[,2]) # looks super white noise! good
mq(res, lag=round(sqrt(length(res)))) # DEF white noise
plot(irf(var_mod,ortho=FALSE,boot=TRUE)) # not too sure how to interpret this?


# Cointegration with gdp and ue
# both I(1) so check if cointegrated with Engle Granger




### Wait, you read the ccf wrong, lets try it the correct way
ccf(y = sld_GDP[,1], x = sd_UE[,1]) # there is some significance at one lag, and it appears at a seasonal
# GDP has some predictive power for UE

# Lets try dynamic model with one lag
# I am NOT going to include an UE terms, unless I did a seasonal one
# dl_1
lag <- 1
ue.e <- embed(sd_UE, dimension = lag+1)
gdp.e <- embed(sld_GDP, dimension = lag+1)
dl1 <- lm(ue.e[,1] ~ gdp.e)
summary(dl1)
acf(dl1$residuals) # still have a seasonal issue
# we could try a seasonal GDP

# Not ideal, lets try an ADL model to get rid of some extra, specifically adding
# a lagged seasonal term of ue
lag <- 4
ue.e <- embed(sd_UE, dimension = lag+1)
gdp.e <- embed(sld_GDP, dimension = lag+1)
adl <- lm(ue.e[,1] ~ ue.e[,5] + gdp.e[,2])
summary(adl)
acf(adl$residuals) # SO DOPE
Box.test(adl$residuals, lag = round(sqrt(length(adl$residuals))), type = "Ljung-Box") # Super white noise

# lets see if we can get soem granger causality
small <- lm(ue.e[,1] ~ gdp.e[,2])
anova(small, adl) # very much so (gotta think about how to spin this)

# Try var model (same as before, just copied)
dat <- cbind(sld_GDP, sd_UE)
VARselect(dat) # 4 lags for BIC (i.e. seasonal is what it wants)
var_mod <- vars::VAR(dat, p=4)
summary(var_mod) # significance at ue_{t-4} and gdp_{t-1}, interesting
res <- resid(var_mod)
layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))
acf(res[,1])
acf(res[,2])
ccf(res[,1],res[,2]) # looks super white noise! good
mq(res, lag=round(sqrt(length(res)))) # DEF white noise
plot(irf(var_mod,ortho=FALSE,boot=TRUE)) # not too sure how to interpret this?


# We could try a little cointegrationm
co <- lm(diff(UE,4) ~ diff(log(GDP), 4))
summary(co)
q <- round(sqrt(length(co$residuals)))
CADFtest(co$residuals,type="drift",criterion="BIC",max.lag.y=q) # not cointgrated
