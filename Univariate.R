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

# Adjust data
cpi <- cpi[,1:2]
cpi <- ts(cpi[,2], start=c(1969, 1), frequency=12)

fx <- fx[, 1:2] # first column is peso:us dollar
fx <- ts(fx[,2], start=c(1990, 1), frequency=12)

gdp <- ts(gdp[,2], start=c(1993, 1), frequency=4) # Quarterly!

# Plot data before compression and window cut
par(mfrow=c(3,1))
plot.ts(cpi, main="CPI")
plot.ts(gdp, main="GDP", ylab="Mex")
plot.ts(fx, main="FX", ylab="Pesos/Dollars")
par(mfrow=c(1,1))

# Decide where to cut window
cpi <- window(cpi, start=2000)
fx <- window(fx, start=2000)
gdp <- window(gdp, start=2000)

# Compress Monthly to Quarterly (takes mean of the 3 months)
cpi <- aggregate(cpi, nfrequency = 4)
fx <- aggregate(fx, nfrequency = 4)

# Plot with all quarterly
par(mfrow=c(3,1))
plot.ts(cpi, main="CPI")
plot.ts(gdp, main="GDP", ylab="Mex")
plot.ts(fx, main="FX", ylab="Pesos/Dollars")
par(mfrow=c(1,1))
# nothing should really change, just slight smoothing
# but now we have the same frequncy for later analysis


#######################################################
#                     GDP
#######################################################
# As GDP, we will work in logs, and a trend so log differences
l_gdp <- log(gdp)
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

# I think a neat way to end the paper would be to make some future economic predictions 
# it will be neat to compare what the univariate prediction vs. the multi are



