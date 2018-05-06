#Final Project
#Applied Time Series
#Roberto Galvan

rm(list = ls())
setwd("~/Documents/STSCI 4550 -- Time Series/TimeSeriesProject")


library(readxl)
library(CADFtest)
library(forecast)

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
ts.plot(GDP,main="GDP",ylab="Millions of Pesos")
ts.plot(CPI,main="CPI",ylab="Index")
ts.plot(FX,main="Foreign Exchange Rate",ylab="Pesos/Dollars")
ts.plot(MB,main="Money Supply",ylab="Millions of Pesos")
ts.plot(UE, main="Unemployment", ylab="Percentage")
ts.plot(IRMEX90,IRUSA90,col= c("blue","red"),main="3-Month Treasury Securities' Rates",ylab="Percentage")
legend("topright",legend=c("Mexican","American"),col= c("blue","red"),lty=1)



##############
######GDP#####
##############
#It has a trend
ts.plot(GDP,main="GDP",ylab="Millions of Pesos")
#Use logs to eliminate trend
L_GDP<-log(GDP)
ts.plot(L_GDP,main="Log of GDP")

#The trend was no eliminated using logs. 
#What kind of trend? Fail to reject, Zt is random walk, L_GDP is Stochastic, work in differences
q1<-floor(sqrt(length(L_GDP)))
CADFtest(GDP, type="trend", criterion = "BIC", max.lag.y = q1)
CADFtest(L_GDP,type="trend",criterion="BIC",max.lag.y=q1)

#Work in differences. It seems to be stationary based on the plots, but that is not very clear
DL_GDP<-diff(L_GDP)
ts.plot(DL_GDP,main="Differences of Log GDP")

#Check for seasonality: There is seasonality. After using the seasonal difference, there is no seasonality left
monthplot(DL_GDP, main="GDP Seasonality Plot", ylab="Deviation")
SDL_GDP<-diff(DL_GDP,4)
monthplot(SDL_GDP)

#Check stationary: Reject, SDL_GDP doesn't have a unit root, SDL_GDP is startionary
q2<-floor(sqrt(length(SDL_GDP)))
CADFtest(SDL_GDP,type="drift",criterion="BIC",max.lag.y=q2)

#Suggesting a model
par(mfrow=c(1,2))
#Non seasonal correlation left in the first quarter and seasonal correlation left in year 1.
acf(SDL_GDP, main="", ylim=c(-.5,1))
#Seasonal correlation left in year 1.
pacf(SDL_GDP, main="", ylim=c(-.5,1))
mtext("GDP Correlations", side=1, line=-24, outer=TRUE, font=2)
par(mfrow=c(1,1))

#MA(1) model: MA(1) coefficient is significant. The ACF shows that correlation was left in the residuals
#Ljung box test: reject the model, residuals are not white noise, the model is not validated. Try another one
MAGDP<-Arima(L_GDP,order=c(0,1,1),seasonal=c(0,1,0),include.constant=TRUE)
summary(MAGDP)
abs(MAGDP$coef/sqrt(diag(MAGDP$var.coef)))
ts.plot(MAGDP$residuals,main="MA(1) Residuals")
acf(MAGDP$residuals,main="Autocorrelations of MA(1) model")
Box.test(MAGDP$residuals,lag=q1,type="Ljung-Box")

#MA(1) & SMA(1) model: Both MA(1) and SMA(1) coefficients are significant.
#Fail to reject, resuduals are white noise, the model is validated
SMAGDP<-Arima(L_GDP,order=c(0,1,1),seasonal=c(0,1,1),include.constant=TRUE)
summary(SMAGDP)
abs(SMAGDP$coef/sqrt(diag(SMAGDP$var.coef)))
ts.plot(SMAGDP$residuals,main="MA(1) & SMA(1) Residuals")
acf(SMAGDP$residuals,main="Autocorrelations of MA(1) & SMA(1) model")
Box.test(SMAGDP$residuals,lag=q1,type="Ljung-Box")

#AR(1) model: AR(1) coefficient is significant.The ACF shows that correlation was left in the residuals
#Ljung box test: reject the model, residuals are not white noise, the model is not validated. Try another one
ARGDP<-Arima(L_GDP,order=c(1,1,0),seasonal=c(0,1,0),include.constant=TRUE)
summary(ARGDP)
abs(ARGDP$coef/sqrt(diag(ARGDP$var.coef)))
ts.plot(ARGDP$residuals,main="AR(1) Residuals")
acf(ARGDP$residuals,main="Autocorrelations of AR(1) model")
Box.test(ARGDP$residuals,lag=q1,type="Ljung-Box")

#SAR(1) model: SAR(1) coefficient is significant.
#Fail to reject, resuduals are white noise, the model is validated
SARGDP<-Arima(L_GDP,order=c(0,1,0),seasonal=c(1,1,0),include.constant=TRUE)
summary(SARGDP)
abs(SARGDP$coef/sqrt(diag(SARGDP$var.coef)))
ts.plot(SARGDP$residuals,main="SAR(1) Residuals")
acf(SARGDP$residuals,main="Autocorrelations of SAR(1) model")
Box.test(SARGDP$residuals,lag=q1,type="Ljung-Box")

# ARS(1,1) AR term and Seasonal AR term
ARS <- Arima(L_GDP,order=c(1,1,0),seasonal=c(1,1,0),include.constant=TRUE)
summary(ARS)
acf(ARS$residuals)
Box.test(ARS$residuals,lag=q1,type="Ljung-Box")

#In sample-comparison (BIC comparison)-Only with the two validated models: MA(1) & SMA(1) and SAR(1)
#The SMA(1) was better (lower BIC)
BIC(SMAGDP)
BIC(SARGDP)

# Out of sample error 
window_forecast<-function(y,order,seasonal,h){
  # Takes in a series, an order, and a seasonal order for the model
  # as well as an h for the prediction step
  # Splits the data at the 75% mark
  # Returns a sequence of errors
  S=round(0.75*length(y))
  error<-c()
  for(i in S:(length(y)-h))
  {
    submod<-Arima(y[1:i],order=order,seasonal=seasonal)
    pred<-predict(submod, n.ahead=h)$pred[h]
    error<-c(error, y[i+h]-pred)
  }
  return(error)
}

# 1-Step Errors
error_MAGDP<-window_forecast(L_GDP,c(0,1,1),c(0,1,0),1)
error_SMAGDP<-window_forecast(L_GDP,c(0,1,1),c(0,1,1),1)
error_ARGDP<-window_forecast(L_GDP,c(1,1,0),c(0,1,0),1)
error_SARGDP<-window_forecast(L_GDP,c(0,1,0),c(1,1,0),1)
error_ARS <- window_forecast(L_GDP,c(1,1,0),c(1,1,0),1)

# Compare MAE
#The SAR(1) model was the best, based on the mean absolute error
MAE<-c(mean(abs(error_MAGDP)),mean(abs(error_SMAGDP)),mean(abs(error_ARGDP)),mean(abs(error_SARGDP)))
rbind(c("MA(1)","MA(1) & SMA(1)","AR(1)","SAR(1)"), round(MAE, 5))

mean(abs(error_ARS))

#Diebold-Mariano test with absolute value loss
#Reject the null hypothesis that states both models perform equally well: Therefore, SAR(1) is better.
dm.test(error_SMAGDP,error_SARGDP,h=1,power=1)

#Predictions 3 years (MA(1) & SMA(1) model)
SMAGDP_P<-predict(SMAGDP,n.ahead=12)
SMAGDP_E<-SMAGDP_P$pred
SMAGDP_L<-SMAGDP_P$pred-qnorm(0.975)*SMAGDP_P$se;
SMAGDP_U<-SMAGDP_P$pred+qnorm(0.975)*SMAGDP_P$se;
cbind(SMAGDP_L,SMAGDP_E,SMAGDP_U)

#Predictions 3 years (SAR(1) model)
SARGDP_P<-predict(SMAGDP,n.ahead=12)
SARGDP_E<-SARGDP_P$pred
SARGDP_L<-SARGDP_P$pred-qnorm(0.975)*SARGDP_P$se;
SARGDP_U<-SARGDP_P$pred+qnorm(0.975)*SARGDP_P$se;
cbind(SMAGDP_L,SMAGDP_E,SMAGDP_U)

#Plot of the two predictions
par(mfrow=c(2,1))
plot.ts(L_GDP,xlim=c(2010,2021),ylim=c(16,17.5),main="Log GDP predictions (MA(1) & SMA(1))",ylab="Log GDP")
lines(SMAGDP_E,col="red")
lines(SMAGDP_L,col="blue")
lines(SMAGDP_U,col="blue")

plot.ts(exp(L_GDP),xlim=c(2010,2021), ylim=c(1.2e7, 3.8e7),main="GDP predictions (SAR)",ylab="Millions of Pesos")
lines(exp(SARGDP_E),col="red")
lines(exp(SARGDP_L),col="blue")
lines(exp(SARGDP_U),col="blue")
legend("topleft", legend=c("Prediction", "95% Confidence Bound"), col=c("red", "blue"), lty=1)
par(mfrow=c(1,1))


#############
#####CPI#####
#############
ts.plot(CPI,main="CPI",ylab="Index")
#Use logs to eliminate trend
L_CPI<-log(CPI)
ts.plot(L_CPI,main="Log of CPI")

#The trend was no eliminated using logs. 
#What kind of trend? Fail to reject, Zt is random walk, L_CPI is stochastic, work in differences
q3<-floor(sqrt(length(L_CPI)))
CADFtest(L_CPI,type="trend",criterion="BIC",max.lag.y=q3)

#Work in differences. It doesn't seem to be stationary
DL_CPI<-diff(L_CPI)
ts.plot(DL_CPI,main="Differences of Log Consumer Price Index")

#Check for seasonality: There is seasonality. After using the seasonal difference, there is no seasonality left
monthplot(DL_CPI)
SDL_CPI<-diff(DL_CPI,4)
monthplot(SDL_CPI)

#Check stationary: Reject, SDL_CPI doesn't have a unit root, SDL_CPI is startionary
q4<-floor(sqrt(length(SDL_CPI)))
CADFtest(SDL_CPI,type="drift",criterion="BIC",max.lag.y=q4)

#Suggesting a model
par(mfrow=c(2,1))
#Non seasonal correlation left in the first quarter and seasonal correlation left in year 1.
acf(SDL_CPI, main="Seasonal and Non-Seansonal CPI Log Differences")
#Seasonal correlation left in year 1.
pacf(SDL_CPI, main="")
par(mfrow=c(1,1))

#MA(1) model: MA(1) coefficient is significant. The ACF shows that correlation was left in the residuals in year 1
#Ljung box test: reject the model, residuals are not white noise, the model is not validated. Try another one
#MA(2) didn't work also
MACPI<-Arima(L_CPI,order=c(0,1,1),seasonal=c(0,1,0),include.constant=TRUE)
abs(MACPI$coef/sqrt(diag(MACPI$var.coef)))
ts.plot(MACPI$residuals,main="MA(1) Residuals")
acf(MACPI$residuals,main="Autocorrelations of MA(1) model")
Box.test(MACPI$residuals,lag=q3,type="Ljung-Box")

#MA(1) & SMA(1) model: Both MA(1) and SMA(1) coefficients are significant.
#Fail to reject, resuduals are white noise, the model is validated
#SMA(1) didn't work
SMACPI<-Arima(L_CPI,order=c(0,1,1),seasonal=c(0,1,1),include.constant=TRUE)
abs(SMACPI$coef/sqrt(diag(SMACPI$var.coef)))
ts.plot(SMACPI$residuals,main="MA(1) & SMA(1) Residuals")
acf(SMACPI$residuals,main="Autocorrelations of MA(1) & SMA(1) model")
Box.test(SMACPI$residuals,lag=q3,type="Ljung-Box")

#AR(1) model: AR(1) coefficient is significant.The ACF shows that correlation was left in the residuals
#Ljung box test: reject the model, residuals are not white noise, the model is not validated. Try another one
ARCPI<-Arima(L_CPI,order=c(1,1,0),seasonal=c(0,1,0),include.constant=TRUE)
abs(ARCPI$coef/sqrt(diag(ARCPI$var.coef)))
ts.plot(ARCPI$residuals,main="AR(1) Residuals")
acf(ARCPI$residuals,main="Autocorrelations of AR(1) model")
Box.test(ARCPI$residuals,lag=q3,type="Ljung-Box")

#AR(1) & SAR(1) model: Both AR(1) and SAR(1) coefficients are significant.
#Fail to reject, resuduals are white noise, the model is validated
#Only SAR(1) didn't work
SARCPI<-Arima(L_CPI,order=c(1,1,0), seasonal=c(1,1,0),include.constant=TRUE)
abs(SARCPI$coef/sqrt(diag(SARCPI$var.coef)))
ts.plot(SARCPI$residuals,main="AR(1) & SAR(1) Residuals")
acf(SARCPI$residuals,main="Autocorrelations of MA(1) & SMA(1) model")
Box.test(SARCPI$residuals,lag=q3,type="Ljung-Box")

#In sample-comparison (BIC comparison)-Only with the two validated models: MA(1) & SMA(1) and AR(1) & SAR(1)
#The AR(1) & SAR(1) was better (lower BIC)
BIC(SMACPI)
BIC(SARCPI)

#Out of sample error, 1-Step Errors
error_MACPI<-window_forecast(L_CPI,c(0,1,1),c(0,1,0),1)
error_SMACPI<-window_forecast(L_CPI,c(0,1,1),c(0,1,1),1)
error_ARCPI<-window_forecast(L_CPI,c(1,1,0),c(0,1,0),1)
error_SARCPI<-window_forecast(L_CPI,c(1,1,0),c(1,1,0),1)

#Compare MAE
#The AR(1) & SAR(1) model was the best, based on the mean absolute error
MAE<-c(mean(abs(error_MACPI)),mean(abs(error_SMACPI)),mean(abs(error_ARCPI)),mean(abs(error_SARCPI)))
rbind(c("MA(1)","MA(1) & SMA(1)","AR(1)","AR(1) & SAR(1)"),MAE)

#Diebold-Mariano test with absolute value loss
#Reject the null hypothesis that states both models perform equally well: Therefore, AR(1) & SAR(1) is better.
dm.test(error_SMACPI,error_SARCPI,h=1,power=1)

#Predictions 3 years (MA(1) & SMA(1) model)
SMACPI_P<-predict(SMACPI,n.ahead=12)
SMACPI_E<-SMACPI_P$pred
SMACPI_L<-SMACPI_P$pred-qnorm(0.975)*SMACPI_P$se;
SMACPI_U<-SMACPI_P$pred+qnorm(0.975)*SMACPI_P$se;
cbind(SMACPI_L,SMACPI_E,SMACPI_U)

#Predictions 3 years (AR(1) & SAR(1) model)
SARCPI_P<-predict(SMACPI,n.ahead=12)
SARCPI_E<-SARCPI_P$pred
SARCPI_L<-SARCPI_P$pred-qnorm(0.975)*SARCPI_P$se;
SARCPI_U<-SARCPI_P$pred+qnorm(0.975)*SARCPI_P$se;
cbind(SMACPI_L,SMACPI_E,SMACPI_U)

#Plot of the two predictions
par(mfrow=c(2,1))
plot.ts(L_CPI,xlim=c(2010,2021),ylim=c(4.5,5.5),main="Log CPI predictions (MA(1) & SMA(1))",ylab="Log CPI")
lines(SMACPI_E,col="red")
lines(SMACPI_L,col="blue")
lines(SMACPI_U,col="blue")

plot.ts(L_CPI,xlim=c(2010,2021),ylim=c(4.5,5.5),main="Log CPI predictions AR(1) & SAR(1)",ylab="Log CPI")
lines(SARCPI_E,col="red")
lines(SARCPI_L,col="blue")
lines(SARCPI_U,col="blue")
par(mfrow=c(1,1))


############
#####UE#####
############
ts.plot(UE,main="UE",ylab="Index")

#What kind of trend? Fail to reject, Zt is random walk, UE is stochastic, work in differences
q5<-floor(sqrt(length(UE)))
CADFtest(UE,type="trend",criterion="BIC",max.lag.y=q5)

#Work in differences. It doesn't seem to be stationary
D_UE<-diff(UE)
ts.plot(D_UE,main="Differences of Unemployment")

#Check for seasonality: There is seasonality. After using the seasonal difference, there is no seasonality left
monthplot(D_UE, main="Seasonality Plot:\nUnemployment in First Differences",
            ylab="Seasonal Variation")
SD_UE<-diff(D_UE,4)
monthplot(SD_UE)

#Check stationary: Reject, SD_UE doesn't have a unit root, SD_UE is startionary
q6<-floor(sqrt(length(SD_UE)))
CADFtest(SD_UE,type="drift",criterion="BIC",max.lag.y=q6)

#Suggesting a model
par(mfrow=c(1,2))
#Non seasonal correlation left in the first quarter and seasonal correlation left in year 1.
acf(SD_UE, main="")
#Seasonal correlation left in year 1.
pacf(SD_UE, main="")
mtext("Unemployment Correlations", side=1, line=-24, outer=TRUE, font=2)

par(mfrow=c(1,1))

#SMA(1) model: The SMA(1) coefficient is significant.
#Fail to reject, resuduals are white noise, the model is validated
SMAUE<-Arima(UE,order=c(0,1,0),seasonal=c(0,1,1),include.constant=TRUE)
abs(SMAUE$coef/sqrt(diag(SMAUE$var.coef)))
ts.plot(SMAUE$residuals,main="SMA(1) Residuals")
acf(SMAUE$residuals,main="Autocorrelations of SMA(1) model")
Box.test(SMAUE$residuals,lag=q5,type="Ljung-Box")

#SAR(1) model: The SAR(1) coefficient is significant.
#Fail to reject, resuduals are white noise, the model is validated
SARUE<-Arima(UE,order=c(0,1,0),seasonal=c(1,1,0),include.constant=TRUE)
abs(SARUE$coef/sqrt(diag(SARUE$var.coef)))
ts.plot(SARUE$residuals,main="SAR(1) Residuals")
acf(SARUE$residuals,main="Autocorrelations of SAR(1) model")
Box.test(SARUE$residuals,lag=q1,type="Ljung-Box")

#In sample-comparison (BIC comparison)-Only with the two validated models: SMA(1) and SAR(1)
#The SMA(1) was better (lower BIC)
BIC(SMAUE)
BIC(SARUE)

#Out of sample error, 1-Step Errors
error_SMAUE<-window_forecast(UE,c(0,1,0),c(0,1,1),1)
error_SARUE<-window_forecast(UE,c(0,1,0),c(1,1,0),1)

# Compare MAE
#The AR(1) & SAR(1) model was the best, based on the mean absolute error
MAE<-c(mean(abs(error_SMAUE)),mean(abs(error_SARUE)))
rbind(c("SMA(1)","SAR(1)"),MAE)

#Diebold-Mariano test with absolute value loss
#Reject the null hypothesis that states both models perform equally well: Therefore, AR(1) & SAR(1) is better.
dm.test(error_SMAUE,error_SARUE,h=1,power=1)

#Predictions 3 years (MA(1) & SMA(1) model)
SMAUE_P<-predict(SMAUE,n.ahead=12)
SMAUE_E<-SMAUE_P$pred
SMAUE_L<-SMAUE_P$pred-qnorm(0.975)*SMAUE_P$se;
SMAUE_U<-SMAUE_P$pred+qnorm(0.975)*SMAUE_P$se;
cbind(SMAUE_L,SMAUE_E,SMAUE_U)

#Predictions 3 years (AR(1) & SAR(1) model)
SARUE_P<-predict(SMAUE,n.ahead=12)
SARUE_E<-SARUE_P$pred
SARUE_L<-SARUE_P$pred-qnorm(0.975)*SARUE_P$se;
SARUE_U<-SARUE_P$pred+qnorm(0.975)*SARUE_P$se;
cbind(SMAUE_L,SMAUE_E,SMAUE_U)

#Plot of the two predictions
par(mfrow=c(2,1))
plot.ts(UE,xlim=c(2010,2021),ylim=c(-0.5,6),main="Unemployment predictions SMA(1)",ylab="Unemployment")
lines(SMAUE_E,col="red")
lines(SMAUE_L,col="blue")
lines(SMAUE_U,col="blue")

plot.ts(UE,xlim=c(2010,2021),ylim=c(-0.5,6),main="Unemployment predictions SAR(1)",ylab="Unemployment")
lines(SARUE_E,col="red")
lines(SARUE_L,col="blue")
lines(SARUE_U,col="blue")
par(mfrow=c(1,1))


#############
###IRMEX90###
#############
ts.plot(IRMEX90,main="IRMEX90",ylab="Percentage")

#What kind of trend? Fail to reject, Zt is stationary, IRMEX90 is deterministic, either work in differences or detrend
q7<-floor(sqrt(length(IRMEX90)))
CADFtest(IRMEX90,type="trend",criterion="BIC",max.lag.y=q7)

#Work in differences. It seems to be stationary
D_IRMEX90<-diff(IRMEX90)
ts.plot(D_IRMEX90,main="Differences of Interest Rates")

#Check for seasonality: There is just a bit of seasonality. After using the seasonal difference, it looks perfect
monthplot(D_IRMEX90)
SD_IRMEX90<-diff(D_IRMEX90,4)
monthplot(SD_IRMEX90)

#Check stationary: Reject, SD_IRMEX90 doesn't have a unit root, SD_IRMEX90 is startionary
q8<-floor(sqrt(length(SD_IRMEX90)))
CADFtest(SD_IRMEX90,type="drift",criterion="BIC",max.lag.y=q8)

#Suggesting a model
par(mfrow=c(2,1))
#Non seasonal correlation left in the first quarter and seasonal correlation left in year 1.
acf(SD_IRMEX90, main="Seasonal and Non-Seansonal Differences of Interest Rates")
#Seasonal correlation left in year 1.
pacf(SD_IRMEX90, main="")
par(mfrow=c(1,1))

#MA(1) & SMA(1) model: Both MA(1) and SMA(1) coefficients are significant.
#Fail to reject, resuduals are white noise, the model is validated
#MA(1),MA(2), MA(3) and SMA(1) didn't work also
SMAIRMEX90<-Arima(IRMEX90,order=c(0,1,1),seasonal=c(0,1,1),include.constant=TRUE)
abs(SMAIRMEX90$coef/sqrt(diag(SMAIRMEX90$var.coef)))
ts.plot(SMAIRMEX90$residuals,main="MA(2) & SMA(1) Residuals")
acf(SMAIRMEX90$residuals,main="Autocorrelations of MA(2) & SMA(1) model")
Box.test(SMAIRMEX90$residuals,lag=q7,type="Ljung-Box")

#AR(1),MA(1) & SMA(1) model: Both MA(1) and SMA(1) coefficients are significant.AR(1) is not
#Fail to reject, resuduals are white noise, the model is validated
#AR(1),AR(1) & SAR(1) and AR(2) & SAR(1) didn't work also
SARIMA1_IRMEX90<-Arima(IRMEX90,order=c(1,1,1),seasonal=c(0,1,1),include.constant=TRUE)
abs(SARIMA1_IRMEX90$coef/sqrt(diag(SARIMA1_IRMEX90$var.coef)))
ts.plot(SARIMA1_IRMEX90$residuals,main="AR(1),MA(1) and SMA(1) Residuals")
acf(SARIMA1_IRMEX90$residuals,main="Autocorrelations of AR(1),MA(1) & SMA(1) model")
Box.test(SARIMA1_IRMEX90$residuals,lag=q7,type="Ljung-Box")

#AR(1),MA(1) & SAR(1) model: Both MA(1) and SAR(1) coefficients are significant.AR(1) is not
#Fail to reject, resuduals are white noise, the model is validated
#AR(1),AR(1) & SAR(1) and AR(2) & SAR(1) didn't work also
SARIMA2_IRMEX90<-Arima(IRMEX90,order=c(1,1,1),seasonal=c(1,1,0),include.constant=TRUE)
abs(SARIMA2_IRMEX90$coef/sqrt(diag(SARIMA2_IRMEX90$var.coef)))
ts.plot(SARIMA2_IRMEX90$residuals,main="AR(1),MA(1) and SAR(1) Residuals")
acf(SARIMA2_IRMEX90$residuals,main="Autocorrelations of AR(1),MA(1) & SAR(1) model")
Box.test(SARIMA2_IRMEX90$residuals,lag=q7,type="Ljung-Box")

#In sample-comparison (BIC comparison) 
#AR(1),MA(1) and SMA(1) model is the best one
BIC(SMAIRMEX90)
BIC(SARIMA1_IRMEX90)
BIC(SARIMA2_IRMEX90)

#Out of sample error, 1-Step Errors
error_SMAIRMEX90<-window_forecast(IRMEX90,c(0,1,1),c(0,1,1),1)
error_SARIMA1_IRMEX90<-window_forecast(IRMEX90,c(1,1,1),c(0,1,1),1)
error_SARIMA2_IRMEX90<-window_forecast(IRMEX90,c(1,1,1),c(1,1,0),1)

#Compare MAE
#AR(1) & MA(1) model is the best one
MAE<-c(mean(abs(error_SMAIRMEX90)),mean(abs(error_SARIMA1_IRMEX90)),mean(abs(error_SARIMA2_IRMEX90)))
rbind(c("AR(1) & SMA(1)","AR(1),MA(1) and SMA(1)","AR(1),MA(1) and SAR(1)"),MAE)


#Diebold-Mariano test with absolute value loss
#Fail to reject the null hypothesis. Therefore both models perform as good as each other.
dm.test(error_SMAIRMEX90,error_SARIMA1_IRMEX90,h=1,power=1)

#Predictions 3 years (MA(1) & SMA(1) model)
SMAIRMEX90_P<-predict(SMAIRMEX90,n.ahead=12)
SMAIRMEX90_E<-SMAIRMEX90_P$pred
SMAIRMEX90_L<-SMAIRMEX90_P$pred-qnorm(0.975)*SMAIRMEX90_P$se;
SMAIRMEX90_U<-SMAIRMEX90_P$pred+qnorm(0.975)*SMAIRMEX90_P$se;
cbind(SMAIRMEX90_L,SMAIRMEX90_E,SMAIRMEX90_U)

#Predictions 3 years (AR(1),MA(1) and SMA(1) model)
SARIMA1_IRMEX90_P<-predict(SARIMA1_IRMEX90,n.ahead=12)
SARIMA1_IRMEX90_E<-SARIMA1_IRMEX90_P$pred
SARIMA1_IRMEX90_L<-SARIMA1_IRMEX90_P$pred-qnorm(0.975)*SARIMA1_IRMEX90_P$se;
SARIMA1_IRMEX90_U<-SARIMA1_IRMEX90_P$pred+qnorm(0.975)*SARIMA1_IRMEX90_P$se;
cbind(SARIMA1_IRMEX90_L,SARIMA1_IRMEX90_E,SARIMA1_IRMEX90_U)

#Plot of the predictions
plot.ts(IRMEX90,xlim=c(2010,2021),ylim=c(-10,25),main="3-month interest rates predictions MA(1) & SMA(1) model",ylab="Percentage")
lines(SMAIRMEX90_E,col="red")
lines(SMAIRMEX90_L,col="blue")
lines(SMAIRMEX90_U,col="blue")

#Plot of the predictions
plot.ts(IRMEX90,xlim=c(2010,2021),ylim=c(-10,25),main="3-month interest rates predictions AR(1),MA(1) and SMA(1)",ylab="Percentage")
lines(SARIMA1_IRMEX90_E,col="red")
lines(SARIMA1_IRMEX90_L,col="blue")
lines(SARIMA1_IRMEX90_U,col="blue")
