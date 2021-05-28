### Aplication for underdisperssion process
##### violence family

y <- c(0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0,
       1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 2, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 0, 1,
       0, 1, 1, 0, 0, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 2,
       2, 0, 0, 0, 0, 0, 0, 1, 0, 2, 1, 0, 0, 1, 3)
ts.plot(y)
ts.plot(diff(y))
acf2(y)
n=length(y)
X=cbind(sin(2*pi*(1:(n))/12),cos(2*pi*(1:(n))/12))

ufit=glm.bg(y~X)
ufit$coefficients
acf2(ufit$resid)
qqPlot(ufit$resid)
shapiro.test(ufit$resid)
X=cbind(sin(2*pi*(1:(n))/12),cos(2*pi*(1:(n))/12))
#X=cbind(cos(2*pi*(1:(n))/12))
#X=as.numeric(y>2)
X_new=cbind(1,X)

#X_new=matrix(1,ncol=1,nrow=length(y))
beta=c(0.1,0.2,0.2);phi=c(0.2);theta=c(0.0);nu=0.8
fitt=bergARMA(beta,phi,theta,nu,y,X_new,link="log")
round(fitt$par,3)
fitt$AIC 
diag(fitt$cov)
table(fitt$par[length(fitt$par)]>abs(fitt$mu-1))
require(astsa)
require(car)
acf2(fitt$resid)
qqPlot(fitt$resid)
shapiro.test(fitt$resid)

#Negative binomial
fit1=garmaFit(y~X,order=c(1,0),family=NBI(log))
fit1$coef
fit1$aic
checkresiduals(fit1)
acf2(fit1$residuals)
qqPlot(fit1$residuals)
shapiro.test(fit1$residuals)

#PO

fit2=garmaFit(y~X,order=c(1,0),family=PO(log))
fit2 
checkresiduals(fit2)
acf2(fit2$residuals)
qqPlot(fit2$residuals)
