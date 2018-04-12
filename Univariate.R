# Univariate.R
# Spring 2018

setwd("~/Documents/STSCI 4550 -- Time Series/TimeSeriesProject")
rm(list = ls())

library(forecast)
library(CADFtest)


########################################################
#                     Data
########################################################
# Load Data
gdp <- read.csv("Data/gdp_2013.csv")
cpi <- read.csv("Data/CPI.csv")
fx <- read.csv("Data/FX.csv")
mex90 <- read.csv("Data/90day_mex.csv")
us90 <- read.csv("Data/90day_us.csv")
pesos <- read.csv("Data/pesos.csv") # money supply
emp <- read.csv("Data/unemployment.csv") 


# Adjust data
cpi <- cpi[,1:2]
cpi <- ts(cpi[,2], start=c(1969, 1), frequency=12)

fx <- fx[, 1:2] # first column is peso:us dollar
fx <- ts(fx[,2], start=c(1990, 1), frequency=12)

gdp <- ts(gdp[,2], start=c(1993, 1), frequency=4) # Quarterly!

mex90 <- ts(mex90[,2], start=c(2000, 1), frequency=4)
us90 <- ts(us90[,2], start=c(2000,1), frequency = 12)

pesos <- ts(pesos[,2], start=c(2000, 1), frequency=12)
pesos <- pesos / 1e6 # in millions

emp <- ts(emp[,2], start=c(2000,1), frequency = 4)

# Decide where to cut window
cpi <- window(cpi, start=2000)
fx <- window(fx, start=2000)
gdp <- window(gdp, start=2000)

# Compress Monthly to Quarterly (takes mean of the 3 months)
cpi <- aggregate(cpi, nfrequency = 4)
fx <- aggregate(fx, nfrequency = 4)
us90 <- aggregate(us90, nfrequency = 4)
pesos <- aggregate(pesos, nfrequency = 4)

# Plot data before compression and window cut
par(mfrow=c(4,2))
plot.ts(cpi, main="CPI")
plot.ts(gdp, main="GDP", ylab="Mex")
plot.ts(fx, main="FX", ylab="Pesos/Dollars")
plot.ts(mex90, main="Interest Rate", ylab="%")
plot.ts(us90, main="Interest Rate", ylab="%")
plot.ts(pesos, main="Money Supply", ylab="Pesos (millions)")
plot.ts(emp, main="Unemployment", ylab="%")
par(mfrow=c(1,1))



#######################################################
#                     GDP
#######################################################
# As GDP, we will work in logs, and a trend so log differences
l_gdp <- log(gdp)

# What kind of trend?
q <- floor(sqrt(length(l_gdp)))
CADFtest(l_gdp, type="trend", criterion = "BIC", max.lag.y = q) # fail to reject --> stochastic

# Work in differences
ld_gdp <- diff(log(gdp))
plot.ts(ld_gdp)

# check for seasonality
monthplot(ld_gdp) # very seasonaly
sld_gdp <- diff(ld_gdp, 4) # looks good now
monthplot(sld_gdp)

# Check stationary
q <- floor(sqrt(length(sld_gdp)))
CADFtest(sld_gdp, type="drift", criterion = "BIC", max.lag.y = q)
# reject, therefore no unit root --> startionary

# Check acf, pacf
par(mfrow=c(2,1))
acf(sld_gdp, main="GDP Seasonal Log Differences")
pacf(sld_gdp, main="")
par(mfrow=c(1,1))

# AR(1)
ar <- Arima(l_gdp, order=c(1,1,0), seasonal=c(0,1,1))
summary(ar)
plot.ts(ar$residuals)
acf(ar$residuals)
Box.test(ar$residuals, lag=q, type="Ljung-Box") # Also very good

# SAR(1)
sar <- Arima(l_gdp, order=c(0,1,0), seasonal=c(1,1,0)) # seasonal and normal differences and 1 SAR term
summary(sar) # is significant term
plot.ts(sar$residuals, main="SAR(1) Residuals") # looks pretty good
acf(sar$residuals, main="Autocorrelations of SAR(1) model")
Box.test(sar$residuals, lag=q, type="Ljung-Box") # fail to reject, so good fit

# AR and SAR model
sar_both <- Arima(l_gdp, order=c(1,1,0), seasonal=c(1,1,0))
summary(sar_both)
plot.ts(sar_both$residuals)
acf(sar_both$residuals)
Box.test(sar_both$residuals, lag=q, type="Ljung-Box")


# MA model
ma <- Arima(l_gdp, order=c(0,1,1), seasonal=c(0,1,0)) # seasonal and normal differences and 1 SAR term
summary(ma) # is significant term
plot.ts(ma$residuals, main="MA(1) Residuals") 
acf(ma$residuals, main="Autocorrelations of MA(1) model")
Box.test(ma$residuals, lag=q, type="Ljung-Box") # fail to reject, so good fit

# SMA model
sma <- Arima(l_gdp, order=c(0,1,1), seasonal=c(0,1,1)) # seasonal and normal differences and 1 SAR term
summary(sma) # is significant term
plot.ts(sma$residuals, main="SMA(1) Residuals") 
acf(sma$residuals, main="Autocorrelations of MA(1) model")
Box.test(sma$residuals, lag=q, type="Ljung-Box") # fail to reject, so good fit

# Compare BIC for all models
BIC(ar) # <-- best fit
BIC(sar)
BIC(sar_both)
BIC(ma)
BIC(sma)

# Out of sample error 
window_forecast <- function(y, order, seasonal, h){
  # Takes in a series, an order, and a seasonal order for the model
  # as well as an h for the prediction step
  # Splits the data at the 75% mark
  # Returns a sequence of error
  S = round(0.75*length(y))
  error <- c()
  for (i in S:(length(y)-h))
  {
    submod <- Arima(y[1:i], order=order, seasonal=seasonal)
    pred <- predict(submod, n.ahead=h)$pred[h]
    error <- c(error, y[i+h] - pred)
  }
  return(error)
}

# 1-step errors
error_ar <- window_forecast(l_gdp, c(1,1,0), c(0,1,0), 1)
error_sar <- window_forecast(l_gdp, c(0,1,0), c(1,1,0), 1)
error_sar_both <- window_forecast(l_gdp, c(1,1,0), c(1,1,0), 1)
error_ma <- window_forecast(l_gdp, c(0,1,1), c(0,1,0), 1)
error_sma <- window_forecast(l_gdp, c(0,1,1), c(0,1,1), 1)

# Compare mse
mses <- c(mean(error_ar), mean(error_sar), mean(error_sar_both),
          mean(error_ma), mean(error_sma))
rbind(c("AR", "SAR", "Both", "MA", "SMA"), mses)
# This says that the seasonal AR model is the best

# Run DM test comparing SAR and AR
h <- 1
dm.test(error_sar, error_ar, h=h, power=2) # reject, therefore SAR is indeed better

# Predictions (SAR model)
sar_pred <- predict(sar, n.ahead=8) # 2 years
sar_expected <- sar_pred$pred
sar_lower <- sar_pred$pred-qnorm(0.975)*sar_pred$se;
sar_upper <- sar_pred$pred+qnorm(0.975)*sar_pred$se;
cbind(sar_lower, sar_expected, sar_upper)

# Predictions (AR model)
ar_pred <- predict(ar, n.ahead=8) # 2 years
ar_expected <- ar_pred$pred
ar_lower <- ar_pred$pred-qnorm(0.975)*ar_pred$se;
ar_upper <- ar_pred$pred+qnorm(0.975)*ar_pred$se;
cbind(ar_lower, ar_expected, ar_upper)

# Plot of the two predictions
par(mfrow=c(2,1))
plot.ts(l_gdp, xlim=c(2000, 2020), ylim=c(16.3, 17), main="Log GDP predictions (SAR)", ylab="Log GDP")
lines(sar_expected,col="red")
lines(sar_lower,col="blue")
lines(sar_upper,col="blue")

plot.ts(l_gdp, xlim=c(2000, 2020), ylim=c(16.3, 17), main="Log GDP predictions (AR)", ylab="Log GDP")
lines(ar_expected,col="red")
lines(ar_lower,col="blue")
lines(ar_upper,col="blue")
par(mfrow=c(1,1))
# Those plots are in fact different, just very slightly

# I think a neat way to end the paper would be to make some future economic predictions 
# it will be neat to compare what the univariate prediction vs. the multi are


#######################################################
#                     CPI
#######################################################
# Once again we will work in logs
l_cpi <- log(cpi)

# Check trend type
q <- floor(sqrt(length(l_cpi)))
CADFtest(l_cpi, type="trend", criterion = "BIC", max.lag.y = q) # reject, stochastic

# Work in differences
dl_cpi <- diff(l_cpi)
CADFtest(dl_cpi, type="drift", criterion = "BIC", max.lag.y = q) # stationary

# Check seasonality
monthplot(dl_cpi) # should work in seasonal differences
sdl_cpi <- diff(dl_cpi, 4)
monthplot(sdl_cpi)
CADFtest(sdl_cpi, type="drift", criterion = "BIC", max.lag.y = q) # super stationary

# Now lets start modeling
par(mfrow=c(2,1))
acf(sdl_cpi, main="CPI Seasonal Log Differences")
pacf(sdl_cpi, main="")
par(mfrow=c(1,1))

# Appears 1 Seasonal MA term (or later, 1 or 2 SAR terms)
sma <- Arima(l_cpi, order=c(0,1,0), seasonal=c(0,1,1))
summary(sma)
plot.ts(sma$residuals)
acf(sma$residuals) # on the edge
Box.test(sma$residuals, lag=q, type="Ljung-Box") # fail to reject, so should be fine

# Check SAR models
# SAR(1)
sar1 <- Arima(l_cpi, order=c(0,1,0), seasonal=c(1,1,0))
summary(sar1)
plot.ts(sar1$residuals)
acf(sar1$residuals) # on the edge
Box.test(sar1$residuals, lag=q, type="Ljung-Box") # reject --> not a good model

# SAR(2)
sar2 <- Arima(l_cpi, order=c(0,1,0), seasonal=c(2,1,0))
summary(sar2)
plot.ts(sar2$residuals)
acf(sar2$residuals) # on the edge
Box.test(sar2$residuals, lag=q, type="Ljung-Box") # fail to reject --> kinda OK model

# Compare in-sample metric
BIC(sma)
BIC(sar1) # <-- says best, but I am very skeptical
BIC(sar2)

# out-of-sample error (this is weird, it is saying they are all the same?)
error_sma <- window_forecast(l_cpi, c(0,1,0), c(0,1,1), 1)
error_sar1 <- window_forecast(l_cpi, c(0,1,0), c(1,1,0), 1)
error_sar2 <- window_forecast(l_cpi, c(0,1,0), c(2,1,0), 1)

# Compare
# Compare mse
mses <- c(mean(error_sma^2), mean(error_sar1^2), mean(error_sar2^2))
rbind(c("SMA", "SAR(1)", "SAR(2)"), mses)

# Should do a DM test, but not sure which to even compare?
h <- 1
dm.test(error_sma, error_sar2, h=h, power=2) # error, wtf?

# Lets just make a predict with one for now? SAR(2)
sar2_pred <- predict(sar2, n.ahead=8) # 2 years
sar2_expected <- sar2_pred$pred
sar2_lower <- sar2_pred$pred-qnorm(0.975)*sar2_pred$se;
sar2_upper <- sar2_pred$pred+qnorm(0.975)*sar2_pred$se;
cbind(sar2_lower, sar2_expected, sar2_upper)

# Plot of the two predictions
plot.ts(l_cpi, xlim=c(2000, 2020), ylim=c(5.2, 6.2), main="Log CPI predictions (SAR2)", ylab="Log CPI")
lines(sar2_expected,col="red")
lines(sar2_lower,col="blue")
lines(sar2_upper,col="blue")



#######################################################
#                     FX
#######################################################
# Log it as % changes
l_fx <- log(fx)

# Check trend type
q <- floor(sqrt(length(l_fx)))
CADFtest(l_fx, type="trend", criterion = "BIC", max.lag.y = q) # fail to reject, stochastic 

# Work in differences
dl_fx <- diff(l_fx)
CADFtest(dl_fx, type="drift", criterion = "BIC", max.lag.y = q) # stationary

# Check seasonality
monthplot(dl_fx) # I think it is fine each well within variation

# Now lets start modeling
par(mfrow=c(2,1))
acf(dl_fx, main="FX Log Differences")
pacf(dl_fx, main="")
par(mfrow=c(1,1))
# There is nothing to be modeleded here?




