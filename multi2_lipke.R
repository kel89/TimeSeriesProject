# mutli2_lipke.R
# Second (cleaner) version of the multivariate analysis

# Questions:
# 1. I am having trouble interpreting coefficients (dynamic model) as our series
# are in differences?

# 2. these series are sort of inversely cointegrated...is there a test for that
# or some sort of transformation we can do then test for cointegration

# 3. The Johnasen approach, I am not failing to reject the top line? what does that mean
# are both of these saying that there is indeed cointegration?

# 4. How do we fit GARCH with a good ARIMA model, i.e. our best model for GDP
# is (1,1,0)(1,1,0), but no seasonal option in Garch?

# 5. GARCH output plots, we need to figure out how to change x-axis to time

rm(list = ls())
setwd("~/Documents/STSCI 4550 -- Time Series/TimeSeriesProject")

library(readxl)
library(CADFtest)
library(forecast)
library(vars)
library(stargazer)
library(MTS)
library(urca)
library(fGarch)


# Read in Data
GDP<- read_excel("Data_Roberto.xlsx",sheet="GDP",skip=17)
GDP<-ts(GDP[,2],start=c(1993,1),frequency=4)
GDP<-window(GDP,start=c(1997,1),end=c(2017,4))

UE<-read_excel("Data_Roberto.xlsx",sheet ="Unemployment",skip=10)
UE<-ts(UE[,2],start=c(1987,1),frequency=4)
UE<-window(UE,start=c(1997,1),end=c(2017,4))

# Transform data (logs, differences, seasonal, etc.)
sld_GDP <- diff(diff(log(GDP), 4))
sd_UE <- diff(diff(UE), 4)

# Following arc of paper, we need to explore UE a bit more
ts.plot(UE)
# No trend, but doesn't look stationary
q <- round(sqrt(length(UE)))
CADFtest(UE,type="drift",criterion="BIC",max.lag.y=q) # not stationary

# Check in differences
ts.plot(diff(UE))
monthplot(diff(UE))
# Check stationarity of seasonal differences
q <- round(sqrt(length(sd_UE)))
CADFtest(sd_UE,type="drift",criterion="BIC",max.lag.y=q) # stationary in first and seasonal differences


##### Dynamic Model -----
# Lets start by investingating predictive power
ccf(y=sd_UE[,1], x=sld_GDP[,1])
# Looking at this, it appears the GDP has predictive power for UE at lag 1, and around seasonality

# Based on this lets try a DL(1) model
lag <- 1
xdat <- embed(sld_GDP, dimension = lag+1)
ydat <- embed(sd_UE, dimension = lag+1)
dl1 <- lm(ydat[,1] ~ xdat) # has contemporaneous and lagged term
summary(dl1) # Very significant terms
acf(dl1$residuals) # still and issue at order of seasonality

# Seasonal DL model
lag <- 4 # order of seasonality
xdat <- embed(sld_GDP, dimension = lag+1)
ydat <- embed(sd_UE, dimension = lag+1)
dls <- lm(ydat[,1] ~ xdat[,1] + xdat[,2] + xdat[,4])
summary(dls)
acf(dls$residuals)
# This did NOT help make the residuals look better
# Not sure what we should do here

# Try DL model with many lags
dl_big <- lm(ydat[,1] ~ xdat)
summary(dl_big)
acf(dl_big$residuals) # no real change

# Autogregressive Distributed Lag Model
# Lets add one y term to the regression
adl <- lm(ydat[,1] ~ ydat[,2] + xdat[,1:2])
summary(adl)
acf(adl$residuals) # This really doesn't help, we have that spike at the order of seasaonlity

# Lets try a big ADL model, with all the y lags?
adl_big <- lm(ydat[,1] ~ ydat[,-1] + xdat)
summary(adl_big)
acf(adl_big$residuals) # looking GOOD
Box.test(adl_big$residuals, lag = q, type = "Ljung-Box") # these are white noise

# Now that we have a validated ADL model, we should test for granger causality
small <- lm(ydat[,1] ~ ydat[,-1])
anova(small, adl_big)
# VERY much granger caused, good!


#### VAR model (for prediction) ----
dat <- cbind(sld_GDP, sd_UE)
VARselect(dat) # 4 (at order of seasonality)
var_mod <- vars::VAR(dat, p=4)
summary(var_mod)
# Notice, the most important terms are GPD at lag 1 (which matches ccf)
# and UE at lag 4, highlighting the seasonal effect we keep seeing
res <- resid(var_mod)
layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))
acf(res[,1])
acf(res[,2])
ccf(res[,1], res[,2])
par(mfrow=c(1,1))
# This looks very good, and we can do an MQ test to verify
mq(res, lag=round(sqrt(length(res)))) # DEF white noise
plot(irf(var_mod,ortho=FALSE,boot=TRUE)) # This is actually interesting
# and it tells a pretty realistic story

#### Cointegration ----
# Now, looking at the graphs I do not see much evidence of cointegration,
# but, apriori, it woul dmake sense for these two series to move together
# Therefore, we will proceed and check for it.

# We already know that in seasonal differences, that these series are of order one
sl_GDP <- diff(log(GDP), 4)
s_UE <- diff(UE,4)
ts.plot(sl_GDP)
ts.plot(s_UE)
# They seem inverted

c_test <- lm(s_UE ~ sl_GDP)
summary(c_test)
q <- round(sqrt(length(c_test$residuals)))
CADFtest(c_test$residuals,type="drift",criterion="BIC",max.lag.y=q)
# Looking at test stat -2.83 is not big enough, this is NOT stationary
# We see that there is no cointegration

# We can try with the Johansen method, but I don't expect any difference
dat <- cbind(sl_GDP, s_UE)
trace_test <- ca.jo(dat, type = "trace", K = 4, ecdet = "const", spec = "transitory")
summary(trace_test) # this would say there IS a cointegration equation???

eig_test <- ca.jo(dat, type = "eigen", K = 4, ecdet = "const", spec = "transitory")
summary(eig_test)


#### GARCH ----
# I do not see much volitlity clustering, but 
# we can just do a quck test to verify this

# I am not even sure what we are supposed to do this on... I am thinking
# Maybe just GDP, to see if it has conditional hetero?

# The (1,1,0)(1,1,0) seemed to be the best model for GDP, lets re-estimate it
# And check for conditional heteroskedasticity
l_GDP <- log(GDP)
sar <- Arima(l_GDP, order=c(1,1,0), seasonal = c(1,1,0))
par(mfrow=c(2,1))
acf(sar$residuals)
acf(sar$residuals^2)
par(mfrow=c(1,1)) # it actually looks like we might have some!
# There is only one significant, so we could probably get away with an ARCH model

# lets do ARCH(1), i.e. GARCH(1,0)
arch <- garchFit(~ arma(4, 0) + garch(1, 0), cond.dist="QMLE", data = sld_GDP) # this is bastardized attempt to work
# with the fGARCH package without a seasonal model
summary(arch)
par(mfrow=c(1,1))
plot(arch, which=10)
plot(arch, which=11)
plot(arch, which=2)
# so we see the most, not surprisingly, at around 2000, and 2008, the greatest 
# periods of financial turmoil
