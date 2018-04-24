library(CADFtest)
library(forecast)

GDP=read.csv("gdp.csv")
GDP=ts(GDP$SR16573,start=c(1993,1),frequency=4)
ts.plot(GDP)
# There is an obvious trend and try to use log to eliminate the trend
LogGDP=log(GDP)
plot(LogGDP)
# Still a trend here
# Test what kind of trend it is
q=sqrt(length(LogGDP))
CADFtest(LogGDP,type="trend",criterion="BIC",max.lag.y=q)
# P value smaller than 0.05 deterministic trend
# Check for seasonality
SGDP=diff(LogGDP,4)
seasonplot(SGDP)
DSGDP=diff(SGDP)
# There is strong seasonality, work in difference
q=sqrt(length(DSGDP))
CADFtest(DSGDP,type="drift",criterion="BIC",max.lag.y=q)
# p value is less than 0.05, stationary
acf(DSGDP)
pacf(DSGDP)
# Use sarima model to test
SMAGDP=Arima(LogGDP,order=c(0,1,1),seasonal=c(0,1,1),include.constant=TRUE)
abs(SMAGDP$coef/sqrt(diag(SMAGDP$var.coef)))
ts.plot(SMAGDP$residuals)
acf(SMAGDP$residuals)
Box.test(SMAGDP$residuals,lag=q,type="Ljung-Box")
# SAR Model
SARGDP=Arima(LogGDP,order=c(0,1,0),seasonal=c(1,1,0),include.constant=TRUE)
abs(SARGDP$coef/sqrt(diag(SARGDP$var.coef)))
ts.plot(SARGDP$residuals)
acf(SARGDP$residuals)
Box.test(SARGDP$residuals,lag=q,type="Ljung-Box")
# Calculate BIC
BIC(SMAGDP)
# Forecast
h=1
S=round(0.75*length(LogGDP))
error<-c()
for(i in S:(length(LogGDP)-h))
  {
    submod<-Arima(LogGDP[1:i],order=c(0,1,1),seasonal=c(0,1,1))
    pred<-predict(submod, n.ahead=h)$pred[h]
    error<-c(error, LogGDP[i+h]-pred)
  }

SMAGDP_P<-predict(SMAGDP,n.ahead=12)
SMAGDP_E<-SMAGDP_P$pred
SMAGDP_L<-SMAGDP_P$pred-qnorm(0.975)*SMAGDP_P$se;
SMAGDP_U<-SMAGDP_P$pred+qnorm(0.975)*SMAGDP_P$se;
cbind(SMAGDP_L,SMAGDP_E,SMAGDP_U)

h=1
S=round(0.75*length(LogGDP))
error<-c()
for(i in S:(length(LogGDP)-h))
{
  submod<-Arima(LogGDP[1:i],order=c(0,1,0),seasonal=c(1,1,0))
  pred<-predict(submod, n.ahead=h)$pred[h]
  error<-c(error, LogGDP[i+h]-pred)
}

SARGDP_P<-predict(SMAGDP,n.ahead=12)
SARGDP_E<-SARGDP_P$pred
SARGDP_L<-SARGDP_P$pred-qnorm(0.975)*SARGDP_P$se;
SARGDP_U<-SARGDP_P$pred+qnorm(0.975)*SARGDP_P$se;
cbind(SARGDP_L,SARGDP_E,SARGDP_U)

