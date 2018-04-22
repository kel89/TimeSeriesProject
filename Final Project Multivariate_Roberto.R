#Final Project
#Applied Time Series
#Roberto Galvan

rm(list = ls())

install.packages("readxl")
install.packages("CADFtest")
install.packages("forecast")
install.packages("vars")
install.packages("MTS")
library(readxl)
library(CADFtest)
library(forecast)
library(vars)
library(MTS)

##############
#####Data#####
##############
GDP<-read_excel("E:/Cornell/Courses/4th semester/ILRST 4550 Applied Time Series Analysis/Project/Data_Roberto.xlsx",sheet="GDP",skip=17)
GDP<-ts(GDP[,2],start=c(1993,1),frequency=4)
GDP<-window(GDP,start=c(1997,1),end=c(2017,4))
#Seasonal difference, normal difference and logs of GDP
SDL_GDP<-diff(diff(log(GDP),4))

CPI<-read_excel("E:/Cornell/Courses/4th semester/ILRST 4550 Applied Time Series Analysis/Project/Data_Roberto.xlsx",sheet ="CPI",skip=17)
CPI<-ts(CPI[,2],start=c(1969,1),frequency=12)
CPI<-aggregate(CPI,nfrequency=4)/3
CPI<-window(CPI,start=c(1997,1),end=c(2017,4))
#Seasonal difference, normal difference and logs of CPI
SDL_CPI<-diff(diff(log(CPI),4))

FX<-read_excel("E:/Cornell/Courses/4th semester/ILRST 4550 Applied Time Series Analysis/Project/Data_Roberto.xlsx",sheet ="FX",skip=17)
FX<-ts(FX[,2],start=c(1968,1),frequency=12)
FX<-aggregate(FX,nfrequency=4)/3
FX<-window(FX,start=c(1997,1),end=c(2017,4))

MB<-read_excel("E:/Cornell/Courses/4th semester/ILRST 4550 Applied Time Series Analysis/Project/Data_Roberto.xlsx",sheet ="Monetary Base",skip=17)
MB<-MB[-c(1),]
MB<-ts(MB[,2],start=c(1986,1),frequency=12)
MB<-aggregate(MB,nfrequency=4)/3
MB<-window(MB,start=c(1997,1),end=c(2017,4))

UE<-read_excel("E:/Cornell/Courses/4th semester/ILRST 4550 Applied Time Series Analysis/Project/Data_Roberto.xlsx",sheet ="Unemployment",skip=10)
UE<-ts(UE[,2],start=c(1987,1),frequency=4)
UE<-window(UE,start=c(1997,1),end=c(2017,4))
#Seasonal difference, normal difference of Unemployment
SD_UE<-diff(diff(UE),4)

RATES<-read_excel("E:/Cornell/Courses/4th semester/ILRST 4550 Applied Time Series Analysis/Project/Data_Roberto.xlsx",sheet ="Interest Rates",skip=12)
IRMEX90<-ts(RATES[,2],start=c(1997,1),frequency=4)
IRMEX90<-window(IRMEX90,start=c(1997,1),end=c(2017,4))
IRUSA90<-ts(RATES[,4],start=c(1997,1),frequency=4)
IRUSA90<-window(IRUSA90,start=c(1997,1),end=c(2017,4))
#Seasonal difference, normal difference of 3-months Mexican Interest Rates
SD_IRMEX90<-diff(diff(IRMEX90),4)

###################
#####GDP & CPI#####
###################
#No relation between GDP and CPI
ccf(x=SDL_CPI[,1],y=SDL_GDP[,1])
#Even that a basic lagged model
lag<-1
#Embed
SDL_GDP_l1<-embed(SDL_GDP,dimension=lag+1)
SDL_CPI_l1<-embed(SDL_CPI,dimension=lag+1)
Extended_CPI<-lm(SDL_GDP_l1[,1]~SDL_GDP_l1[,-1]+SDL_CPI_l1[,-1])
acf(Extended_CPI$residuals)
#Reject, the model is not validated, residuals are not white noise
Box.test(Extended_CPI$residuals,lag=round(sqrt(length(Extended_CPI$residuals))),type="Ljung-Box")
#All this didn't work

#Fail to reject the null (Xt does not Granger cause Yt). Therefore, CPI does not Granger cause GDP
Reduced_CPI<-lm(SDL_GDP_l1[,1]~SDL_GDP_l1[,-1])
anova(Extended_CPI,Reduced_CPI)

#VAR
m1<-cbind(SDL_GDP,SDL_CPI)
#AIC:9, BIC(SC):4
VARselect(m1)
varfit1<-vars:::VAR(m1,p=4)

#ACF of GDP ok, ACF of CPI ok, CCF ok (contemporaneous is allowed)
varresid1<-resid(varfit1)
par(mfrow=c(2,2))
acf(varresid1[,1])
acf(varresid1[,2])
ccf(varresid1[,1],varresid1[,2])
par(mfrow=c(1,1))
#The null R1=...=R8=0. The alternative that at least one Ri is different from 0.We fail to reject, the model is validated.
mq(varresid1,lag=floor(sqrt(dim(varresid1)[1])))

#No significant results of CPI over GDP (zero line es inside the interval)
irf_var1<-irf(varfit1,ortho=FALSE,boot=TRUE)
plot(irf_var1)


##################
#####GDP & MB#####
##################
ts.plot(CPI,main="Monetary Base",ylab="Millions of Pesos")
#Use logs to eliminate trend
L_MB<-log(MB)
ts.plot(L_MB,main="Log of MB")

#The trend was no eliminated using logs. 
#What kind of trend? Fail to reject, Zt is random walk, L_MB is stochastic, work in differences
q9<-floor(sqrt(length(L_MB)))
CADFtest(L_MB,type="trend",criterion="BIC",max.lag.y=q9)

#Work in differences. It doesn't seem to be stationary
DL_MB<-diff(L_MB)
ts.plot(DL_MB,main="Differences of Log Monetary Base")

#Check for seasonality: There is seasonality. After using the seasonal difference, there is no seasonality left
monthplot(DL_MB)
SDL_MB<-diff(DL_MB,4)
monthplot(SDL_MB)

#Check stationary: Reject, SDL_MB doesn't have a unit root, SDL_MB is startionary
q10<-floor(sqrt(length(SDL_MB)))
CADFtest(SDL_MB,type="drift",criterion="BIC",max.lag.y=q10)

#It seems that GDP lags one period MB (which is the opposite of what we wanted)
ccf(x=SDL_MB[,1],y=SDL_GDP[,1])
#Even that I will do MB lagging GDP one period
lag<-1
#Embed
SDL_GDP_l1<-embed(SDL_GDP,dimension=lag+1)
SDL_MB_l1<-embed(SDL_MB,dimension=lag+1)
Extended_MB<-lm(SDL_GDP_l1[,1]~SDL_GDP_l1[,-1]+SDL_MB_l1[,-1])
acf(Extended_MB$residuals)
#Reject, the model is not validated, residuals are not white noise
Box.test(Extended_MB$residuals,lag=round(sqrt(length(Extended_MB$residuals))),type="Ljung-Box")
#All this didn't work

#Fail to reject the null (Xt does not Granger cause Yt). Therefore, MB does not Granger cause GDP
Reduced_MB<-lm(SDL_GDP_l1[,1]~SDL_GDP_l1[,-1])
anova(Extended_MB,Reduced_MB)

#VAR
m2<-cbind(SDL_GDP,SDL_MB)
#AIC:9, BIC(SC):1
VARselect(m2)
varfit2<-vars:::VAR(m2,p=1)

#ACF of GDP +-, ACF of MB +-, CCF +- I dont think it will be validated
varresid2<-resid(varfit2)
par(mfrow=c(2,2))
acf(varresid2[,1])
acf(varresid2[,2])
ccf(varresid2[,1],varresid2[,2])
par(mfrow=c(1,1))
#The null R1=...=R8=0. The alternative that at least one Ri is different from 0. Reject, the model is not validated.
mq(varresid2,lag=floor(sqrt(dim(varresid2)[1])))

#No significant results of MB over GDP (zero line es inside the interval)
irf_var2<-irf(varfit2,ortho=FALSE,boot=TRUE)
plot(irf_var2)

##################
#####GDP & UE#####
##################
#It seems that GDP lags one period UE (which is the opposite of what we wanted)
ccf(x=SD_UE[,1],y=SDL_GDP[,1])
#Even that I will do UE lagging GDP one period
lag<-1
#Embed
SDL_GDP_l1<-embed(SDL_GDP,dimension=lag+1)
SD_UE_l1<-embed(SD_UE,dimension=lag+1)
Extended_UE<-lm(SDL_GDP_l1[,1]~SDL_GDP_l1[,-1]+SD_UE_l1[,-1])
acf(Extended_UE$residuals)
#Reject, the model is not validated, residuals are not white noise
Box.test(Extended_UE$residuals,lag=round(sqrt(length(Extended_UE$residuals))),type="Ljung-Box")
#All this didn't work

#Fail to reject the null (Xt does not Granger cause Yt). Therefore, UE does Granger cause GDP
Reduced_UE<-lm(SDL_GDP_l1[,1]~SDL_GDP_l1[,-1])
anova(Extended_UE,Reduced_UE)

#VAR
m3<-cbind(SDL_GDP,SD_UE)
#AIC:9, BIC(SC):1
VARselect(m3)
varfit3<-vars:::VAR(m3,p=1)

#ACF of GDP +-, ACF of UE +-, CCF +- I dont think it will be validated
varresid3<-resid(varfit3)
par(mfrow=c(2,2))
acf(varresid3[,1])
acf(varresid3[,2])
ccf(varresid3[,1],varresid3[,2])
par(mfrow=c(1,1))
#The null R1=...=R8=0. The alternative that at least one Ri is different from 0.Reject, the model is not validated.
mq(varresid3,lag=floor(sqrt(dim(varresid3)[1])))

#No significant results of UE over GDP (zero line es inside the interval)
irf_var3<-irf(varfit3,ortho=FALSE,boot=TRUE)
plot(irf_var3)

#######################
#####GDP & MB + UE#####
#######################
MB_UE<-Arima(SDL_GDP,order=c(1,0,0),seasonal=c(1,0,0),xreg=cbind(SDL_MB,SD_UE))
summary(MB_UE)
acf(MB_UE$residuals)
#Model is validated
#UE vs GDP (-)ok, MB vs GDP (-)?????? In the short run it should be (+)
Box.test(MB_UE$residuals,lag=round(sqrt(length(MB_UE$residuals))),type="Ljung-Box")

par(mfrow=c(2,1))
ccf(x=SDL_CPI[,1],y=SDL_GDP[,1])
ccf(x=SDL_MB[,1],y=SDL_GDP[,1])
par(mfrow=c(1,1))

lag<-1
SDL_GDP_l1<-embed(SDL_GDP,dimension=lag+1)
SDL_MB_l1<-embed(SDL_MB,dimension=lag+1)
SD_UE_l1<-embed(SD_UE,dimension=lag+1)
Extended_MB_UE<-lm(SDL_GDP_l1[,1]~SDL_GDP_l1[,-1]+SDL_MB_l1[,-1]+SD_UE_l1[,-1])
acf(Extended_MB_UE$residuals)
#Reject, the model is not validated, residuals are not white noise
Box.test(Extended_MB_UE$residuals,lag=round(sqrt(length(Extended_MB_UE$residuals))),type="Ljung-Box")
#Didn't work


#VAR
#It seems that GDP lags one period UE (which is the opposite of what we wanted)
m4<-cbind(SDL_GDP,SDL_MB,SD_UE)
#AIC=9, BIC=1
VARselect(m4)
varfit4<-vars::VAR(m4,p=1)
#UE vs GDP (-)ok, MB vs GDP (-)?????? In the short run it should be (+)
summary(varfit4)
varresid4<-resid(varfit4)
par(mfrow=c(3,2))
acf(varresid4[,1])
acf(varresid4[,2])
acf(varresid4[,3])
ccf(varresid4[,1], varresid4[,2])
ccf(varresid4[,1], varresid4[,3])
par(mfrow=c(1,1))
#The null R1=...=R8=0. The alternative that at least one Ri is different from 0.Reject, the model is not validated.
mq(varresid4,lag=floor(sqrt(length(SDL_GDP))))
#Non significan effects of MB over GDP, no signiificant effects of UB over GDP
irf_var4<-irf(varfit4,ortho=FALSE,boot=TRUE)
plot(irf_var4)