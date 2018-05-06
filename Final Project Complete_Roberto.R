#Final Project
#Applied Time Series
#Roberto Galvan

rm(list = ls())

install.packages("readxl")
install.packages("CADFtest")
install.packages("forecast")
install.packages("vars")
install.packages("MTS")
install.packages("fGarch")
library(readxl)
library(CADFtest)
library(forecast)
library(vars)
library(MTS)
library(fGarch)


################
##### Data #####
################
GDP<-read_excel("E:/Cornell/Courses/4th semester/ILRST 4550 Applied Time Series Analysis/Project/Data_Roberto.xlsx",sheet ="GDP",skip=17)
GDP<-ts(GDP[,2],start=c(1993,1),frequency=4)
GDP<-window(GDP,start=c(1997,1),end=c(2017,4))

#It has a trend
ts.plot(GDP,main="GDP",ylab="Millions of Pesos")
#Use logs to eliminate trend
L_GDP<-log(GDP)
ts.plot(L_GDP,main="Log of GDP")

#The trend was no eliminated using logs. 
#What kind of trend? Reject, Zt is stationary, L_GDP is deterministic, work in differences
q1<-floor(sqrt(length(L_GDP)))
CADFtest(L_GDP,type="trend",criterion="BIC",max.lag.y=q1)

#Work in differences. It seems to be stationary based on the plots, but that is not very clear
DL_GDP<-diff(L_GDP)
ts.plot(DL_GDP,main="Differences of Log GDP")

#Check for seasonality: There is seasonality. After using the seasonal difference, there is no seasonality left
monthplot(DL_GDP,main="Seasonality Plot: GDP in First Differences",ylab="Deviation")
#Seasonal difference, normal difference and logs of GDP
SDL_GDP<-diff(diff(log(GDP),4))
SL_GDP<-diff(log(GDP),4)
monthplot(SDL_GDP,main="Seasonal and non-seasonal log differences of GDP",ylab="Deviation")

#Check stationary: Reject, SDL_GDP doesn't have a unit root, SDL_GDP is startionary
q2<-floor(sqrt(length(SDL_GDP)))
CADFtest(SDL_GDP,type="drift",criterion="BIC",max.lag.y=q2)

#Suggesting a model
par(mfrow=c(2,1))
#Non-seasonal correlation left in the first quarter and seasonal correlation left in year 1.
acf(SDL_GDP, main="Correlations of the GDP Seasonal and Non-Seansonal Log Differences")
#Seasonal correlation left in year 1.
pacf(SDL_GDP, main="")
par(mfrow=c(1,1))

#MA(1) model: MA(1) coefficient is significant. The ACF shows that correlation was left in the residuals at seasonality order
#Ljung box test: reject the model, residuals are not white noise, the model is not validated. Try another one
MA<-Arima(L_GDP,order=c(0,1,1),seasonal=c(0,1,0),include.constant=TRUE)
abs(MA$coef/sqrt(diag(MA$var.coef)))
ts.plot(MA$residuals,main="MA(1) Residuals")
acf(MA$residuals,main="Autocorrelations of MA(1) model")
Box.test(MA$residuals,lag=q1,type="Ljung-Box")

#MA(1) & SMA(1) model: Both MA(1) and SMA(1) coefficients are significant.
#Fail to reject, resuduals are white noise, the model is validated
MASMA<-Arima(L_GDP,order=c(0,1,1),seasonal=c(0,1,1),include.constant=TRUE)
abs(MASMA$coef/sqrt(diag(MASMA$var.coef)))
ts.plot(MASMA$residuals,main="MA(1) & SMA(1) Residuals")
acf(MASMA$residuals,main="Autocorrelations of MA(1) & SMA(1) model")
Box.test(MASMA$residuals,lag=q1,type="Ljung-Box")

#AR(1) model: AR(1) coefficient is significant.The ACF shows that correlation was left in the residuals at seasonality order
#Ljung box test: reject the model, residuals are not white noise, the model is not validated. Try another one
AR<-Arima(L_GDP,order=c(1,1,0),seasonal=c(0,1,0),include.constant=TRUE)
abs(AR$coef/sqrt(diag(AR$var.coef)))
ts.plot(AR$residuals,main="AR(1) Residuals")
acf(AR$residuals,main="Autocorrelations of AR(1) model")
Box.test(AR$residuals,lag=q1,type="Ljung-Box")

#AR(1) & SAR(1) model: Both AR(1) and SAR(1) coefficients are significant.
#Fail to reject, resuduals are white noise, the model is validated
ARSAR<-Arima(L_GDP,order=c(1,1,0),seasonal=c(1,1,0),include.constant=TRUE)
abs(ARSAR$coef/sqrt(diag(ARSAR$var.coef)))
ts.plot(ARSAR$residuals,main="Ar(1) & SAR(1) Residuals")
acf(ARSAR$residuals,main="Autocorrelations of Ar(1) & SAR(1) model")
Box.test(ARSAR$residuals,lag=q1,type="Ljung-Box")

#SAR(1) model: SAR(1) coefficient is significant.
#Fail to reject, resuduals are white noise, the model is validated
SAR<-Arima(L_GDP,order=c(0,1,0),seasonal=c(1,1,0),include.constant=TRUE)
abs(SAR$coef/sqrt(diag(SAR$var.coef)))
ts.plot(SAR$residuals,main="SAR(1) Residuals")
acf(SAR$residuals,main="Autocorrelations of SAR(1) model")
Box.test(SAR$residuals,lag=q1,type="Ljung-Box")

#In sample-comparison (BIC comparison) along with the four models
#The SMA(1) was better (lower BIC)
BIC(MA)
BIC(MASMA)
BIC(AR)
BIC(ARSAR)
BIC(SAR)
AIC(MA)
AIC(MASMA)
AIC(AR)
AIC(ARSAR)
AIC(SAR)

#Out of sample error 
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

#1-Step Errors
error_MA<-window_forecast(L_GDP,c(0,1,1),c(0,1,0),1)
error_MASMA<-window_forecast(L_GDP,c(0,1,1),c(0,1,1),1)
error_AR<-window_forecast(L_GDP,c(1,1,0),c(0,1,0),1)
error_ARSAR<-window_forecast(L_GDP,c(1,1,0),c(1,1,0),1)
error_SAR<-window_forecast(L_GDP,c(0,1,0),c(1,1,0),1)

#Compare MAE
#The SAR(1) model was the best, based on the mean absolute error
MAE<-c(mean(abs(error_MA)),mean(abs(error_MASMA)),mean(abs(error_AR)),mean(abs(error_ARSAR)),mean(abs(error_SAR)))
rbind(c("MA(1)","MA(1) & SMA(1)","AR(1)","AR(1) & SAR(1)","SAR(1)"),MAE)

#Diebold-Mariano test with absolute value loss
#Reject the null hypothesis that states both models perform equally well: Therefore, SAR(1) is better.
dm.test(error_MASMA,error_SAR,h=1,power=1)

#Predictions 3 years (SAR(1) model)
SAR_P<-predict(SAR,n.ahead=12)
SAR_E<-SAR_P$pred
SAR_L<-SAR_P$pred-qnorm(0.975)*SAR_P$se;
SAR_U<-SAR_P$pred+qnorm(0.975)*SAR_P$se;
cbind(SAR_L,SAR_E,SAR_U)

#Plot of the Predictions
par(mfrow=c(1,1))
plot.ts(L_GDP,xlim=c(2010,2020),ylim=c(16,17.5),main="GDP Predictions, SAR(1) model",ylab="Log of Millons of Pesos")
legend("bottomright",legend=c("Prediction","95% Confidence Bound"),col= c("red","blue"),lty=1)
lines(SAR_E,col="red")
lines(SAR_L,col="blue")
lines(SAR_U,col="blue")


########################
##### Unemployment #####
########################
UE<-read_excel("E:/Cornell/Courses/4th semester/ILRST 4550 Applied Time Series Analysis/Project/Data_Roberto.xlsx",sheet ="Unemployment",skip=10)
UE<-ts(UE[,2],start=c(1987,1),frequency=4)
UE<-window(UE,start=c(1997,1),end=c(2017,4))

#It has a trend
ts.plot(UE,main="Unemployment",ylab="Unemployment Rate (%)")

#What kind of trend? Fail to reject, Zt is random walk, UE is stochastic, work in differences
q3<-floor(sqrt(length(UE)))
CADFtest(UE,type="trend",criterion="BIC",max.lag.y=q3)

#Work in differences. It seems to be stationary, but is not clear
D_UE<-diff(UE)
ts.plot(D_UE,main="Differences of Unemployment")

#Check for seasonality: There is seasonality. After using the seasonal difference, there is no seasonality left
monthplot(D_UE,main="Seasonality Plot: Unemployment in First Differences",ylab="Deviation")
SD_UE<-diff(diff(UE),4)
S_UE<-diff(UE,4)
monthplot(SD_UE)

#Check stationary: Reject, SD_UE doesn't have a unit root, SD_UE is startionary
q4<-floor(sqrt(length(SD_UE)))
CADFtest(SD_UE,type="drift",criterion="BIC",max.lag.y=q4)


####################
##### GDP & UE #####
####################
GDPG<-(GDP/lag(GDP,k=-1))-1
OLR<-lm(GDPG~UE[2:84])
summary(OLR)
plot(UE[2:84],GDPG,main="Okun's Law Regression",xlab="Unemployment Rate (%)",ylab="GDP Growth (%)",pch=16,cex=1.3,col="blue")
abline(OLR,col="red")

#It seems that GDP lags one period UE (which is the opposite of what we wanted)
ccf(x=SDL_GDP[,1],y=SD_UE[,1], main="Cross Correlation: y=UE & x=GDP")


################################################
##### Dynamic Models and Granger Causality #####
################################################
#I will do GDP lagging UE one period
lag<-4
#Embed
SD_UE_l<-embed(SD_UE,dimension=lag+1)
SDL_GDP_l<-embed(SDL_GDP,dimension=lag+1)
Extended_1<-lm(SD_UE_l[,1]~SDL_GDP_l[,1]+SDL_GDP_l[,2])
summary(Extended_1)
acf(Extended_1$residuals,main="ACF Residuals Dynamic Model (1)")
#Reject, the model is not validated, residuals are not white noise
Box.test(Extended_1$residuals,lag=round(sqrt(length(Extended_1$residuals))),type="Ljung-Box")
#All this didn't work

Extended_2<-lm(SD_UE_l[,1]~SDL_GDP_l[,1]+SDL_GDP_l[,2]+SDL_GDP_l[,5])
summary(Extended_2)
acf(Extended_2$residuals,main="ACF Residuals Dynamic Model (2)")
#Reject, the model is not validated, residuals are not white noise
Box.test(Extended_2$residuals,lag=round(sqrt(length(Extended_2$residuals))),type="Ljung-Box")
#All this didn't work

Extended_3<-lm(SD_UE_l[,1]~SDL_GDP_l)
summary(Extended_3)
acf(Extended_3$residuals,main="ACF Residuals Dynamic Model (3)")
#Reject, the model is not validated, residuals are not white noise
Box.test(Extended_3$residuals,lag=round(sqrt(length(Extended_3$residuals))),type="Ljung-Box")
#All this didn't work

Extended_4<-lm(SD_UE_l[,1]~+SD_UE_l[,-1]+SDL_GDP_l[,-1])
summary(Extended_4)
acf(Extended_4$residuals,main="ACF Residuals Dynamic Model (4)")
#Fail to reject, the model is validated, residuals are white noise
Box.test(Extended_4$residuals,lag=round(sqrt(length(Extended_4$residuals))),type="Ljung-Box")
#This model worked

Extended_5<-lm(SD_UE_l[,1]~SD_UE_l[,5]+SDL_GDP_l[,1]+SDL_GDP_l[,2])
summary(Extended_5)
acf(Extended_5$residuals,main="ACF Residuals Dynamic Model (5)")
#Fail to reject, the model is validated, residuals are white noise
Box.test(Extended_5$residuals,lag=round(sqrt(length(Extended_5$residuals))),type="Ljung-Box")
#This model also worked and was the most parsimonious

#Fail to reject the null (Xt does not Granger cause Yt). Therefore, UE does Granger cause GDP
Reduced_4<-lm(SD_UE_l[,1]~SD_UE_l[,-1])
anova(Extended_4,Reduced_4)

###########
### VAR ###
###########
m1<-cbind(SD_UE,SDL_GDP)
#AIC:4, BIC(SC):4
VARselect(m1)
varfit<-vars:::VAR(m1,p=4)

#ACF of UE, ACF of GDP, CCF, I dont think it will be validated
varresid<-resid(varfit)
par(mfrow=c(2,2))
acf(varresid[,1],main="Unemployment residuals")
acf(varresid[,2],main="GDP residuals")
ccf(varresid[,1],varresid[,2],main="Residials Cross Correlation")
par(mfrow=c(1,1))
#The null R1=...=R8=0. The alternative that at least one Ri is different from 0. Fail to reject, the model is validated.
mq(varresid,lag=floor(sqrt(dim(varresid)[1])))

#No significant results of UE over GDP (zero line es inside the interval)
irf_var<-irf(varfit,ortho=FALSE,boot=TRUE)
plot(irf_var)

#######################
#### Cointegration ####
#######################
#Time series plotted
ts.plot(SDL_GDP,SD_UE,main="Time Series Cointegration",col=c("black","red"))
legend("bottomright",legend=c("GDP","Unemployment"),col=c("black","red"),lty=1)

#Engle-Granger approach
coi<-lm(SL_GDP~S_UE)
res_coi<-coi$residuals
summary(coi)
p<-round(sqrt(length(coi$residuals)))
CADFtest(coi$residuals,type="drift",criterion="BIC",max.lag.y=p)
#In this case we reject the null hypothesis that says that et has a unit root-> et is not stationary
#Therefore, et doesn't have a unit root, et is stationary, and therefore GDP and Unnemployment are cointegrated

#Error correcting model
#As UE and GDP were cointegrated, then run DeltaUE~DeltaGDP+residuals_lagged_1
dSL_GDP<-diff(SL_GDP)
dS_UE<-diff(S_UE)
#residuals_lagged_1
ECT<-res_coi[-length(res_coi)]
#DeltaYt~DeltaXt+residuals_lagged_1
fit_ecm<-lm(dSL_GDP~dS_UE+ECT)

#Based on the graphs and the formal test, we cannot vaildate the model (reject the null, residuals are not WN)
p2<-round(sqrt(length(fit_ecm$residuals)))
plot.ts(fit_ecm$residuals,main="Residuals from the Error Correcting Model", ylab="ECM Residuals")
acf(fit_ecm$residuals,main="ACF Residuals from the Error Correcting Model")
Box.test(fit_ecm$residuals,lag=p2,type="Ljung-Box")

#20.19% of the variability of the model is explained by the model (the dependent variables: the seasonal UE and the ECT)
#The error term is negative and significant, as short term effects return to the long run effects
#Therefore, et doesn't have a unit root, et is stationary-> Unemployment and GDP are cointegrated, but the ECM didn't work
summary(fit_ecm)

#Johansen approach
#We reject in all cases, no cointegration
jm<-cbind(S_UE,SL_GDP)
trace_test<-ca.jo(jm,type="trace",K=4,ecdet="const",spec="transitory")
summary(trace_test)
maxeigen_test<-ca.jo(jm,type="eigen",K=4,ecdet="const", spec="transitory")
summary(maxeigen_test)

fit_vecm1<-cajorls(trace_test,r=1)
fit_vecm1

################
#### GARCH #####
################
#SAR(1).The (0,1,0)(1,1,0) seemed to be the best model for GDP, lets re-estimate it
#And check for conditional heteroskedasticity
ARSAR<-Arima(L_GDP,order=c(0,1,0),seasonal=c(1,1,0),include.constant=TRUE)
par(mfrow=c(2,1))
acf(ARSAR$residuals)
#Significant autocorrelations in the squared residuals
acf(ARSAR$residuals^2)
par(mfrow=c(1,1))
Box.test((ARSAR$residuals^2),lag=round(sqrt(length(ARSAR$residuals))),type="Ljung-Box")
#There is only one significant, so we could probably get away with an ARCH model

#lets do ARCH(1), i.e. GARCH(1,1)
ARCH<-garchFit(~arma(4,0)+garch(1,1),cond.dist="QMLE",data=SDL_GDP)
# with the fGARCH package without a seasonal model
summary(ARCH)
par(mfrow=c(1,1))
plot(ARCH, which=10)
plot(ARCH, which=11)
plot(ARCH, which=2)
plot(ARCH, which=13)
#Greatest volatility was observed in 2000, and 2008 (financial crisis)
