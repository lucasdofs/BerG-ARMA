### Aplication for underdisperssion process
##### violence family
require(astsa)
require(car)
require(forecast)
require(ggplot2)
ycomp <- c(0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0,
       1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 2, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 0, 1,
       0, 1, 1, 0, 0, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 2,
       2, 0, 0, 0, 0, 0, 0, 1, 0, 2, 1, 0, 0, 1, 3)
y<- ycomp[1:120]
ts.plot(y,ylab="Monthly counts")
ts.plot(diff(y))
acf2(ycomp[1:120],max.lag=10,main="")
ggseasonplot(x=ts(y,frequency=12))
ggsubseriesplot(x=ts(y,frequency=12))
n=length(y)
X=cbind(sin(2*pi*(1:(n))/12))

#X=cbind(cos(2*pi*(1:(n))/52))

ufit=glm.bg(y~1)
ufit$coefficients
acf2(ufit$resid,main="")
qqPlot(ufit$resid)
shapiro.test(ufit$resid)

#X=cbind(cos(2*pi*(1:(n))/12))
#X=as.numeric(y>2)
X_new=cbind(1,X)

#Com seno e cos
beta=c(0.1);phi=c(0.1,0.2);theta=c(0.0);nu=0.95
fitt=bergARMA(beta,phi,theta,nu,y,rep(1,n),link="log")
est=round(fitt$par,3);est
fitt$AIC;fitt$BIC
sd_est=diag(abs(fitt$cov))^(1/2);sd_est
sd_est[4]=sd_est[4]+0.1
est[1]+c(-1,1)*sd_est[1]*1.96
est[2]+c(-1,1)*sd_est[2]*1.96
est[3]+c(-1,1)*sd_est[3]*1.96
est[4]+c(-1,1)*(sd_est[4]+0.1)*1.96
est[5]+c(-1,1)*sd_est[5]*1.96

#Teste Hipotese
round(2*min(pnorm(est[1]/sd_est[1]), 1-pnorm(est[1]/sd_est[1])),3)
round(2*min(pnorm(est[2]/sd_est[2]), 1-pnorm(est[2]/sd_est[2])),3)
round(2*min(pnorm(est[3]/sd_est[3]), 1-pnorm(est[3]/sd_est[3])),3)
round(2*min(pnorm(est[4]/sd_est[4]), 1-pnorm(est[4]/sd_est[4])),3)
round(2*min(pnorm((est[4]-1)/sd_est[4]), 1-pnorm((est[4]-1)/sd_est[4])),3)


table(fitt$par[length(fitt$par)]>abs(fitt$mu-1))

acf2(fitt$resid)
qqPlot(fitt$resid)
shapiro.test(fitt$resid)


#Residuals chart
var(fitt$resid);mean(fitt$res)
plot(fitt$resid,ylab="Quantile Residuals")
abline(h=-3)
aa=seq(-5,n+6,by=1)
par(mfrow=c(1,1))
par(mar=c(2.8, 2.7, 1.2, 1)) # margens c(baixo,esq,cima,direia)
par(mgp=c(1.7, 0.45, 0))
plot(fitt$resid,main=" ",xlab="Index",ylab="Quantile residuals", pch = "+",ylim=c(-4,4))
lines(aa,rep(-3,n+12),lty=2,col=1)
lines(aa,rep(3,n+12),lty=2,col=1)
lines(aa,rep(-2,n+12),lty=3,col=1)
lines(aa,rep(2,n+12),lty=3,col=1)


acf2(fitt$resid,main="",max.lag = 15)
qqPlot(fitt$resid,ylab="Quantile residuals")
shapiro.test(fitt$resid)
Box.test(fitt$resid,lag=15)

#Fitting the other models
#Negative binomial
fit1=garmaFit(y~1,order=c(2,0),family=NBI(log))
fit1$coef
fit1$aic;fit1$sbc
round(diag(fit1$vcov)^(1/2),4)
checkresiduals(fit1)
acf2(fit1$residuals)
qqPlot(fit1$residuals)
shapiro.test(fit1$residuals)
round(2*min(pnorm(fit1$coef[1]/sqrt(diag(fit1$vcov))[1]), 1-pnorm(fit1$coef[1]/sqrt(diag(fit1$vcov))[1])),4)
round(2*min(pnorm(fit1$coef[2]/sqrt(diag(fit1$vcov))[2]), 1-pnorm(fit1$coef[2]/sqrt(diag(fit1$vcov))[2])),4)
round(2*min(pnorm(fit1$coef[3]/sqrt(diag(fit1$vcov))[3]), 1-pnorm(fit1$coef[3]/sqrt(diag(fit1$vcov))[3])),4)
round(2*min(pnorm(fit1$coef[4]/sqrt(diag(fit1$vcov))[4]), 1-pnorm(fit1$coef[4]/sqrt(diag(fit1$vcov))[4])),4)






#ZIP
fit2=garmaFit(y~1,order=c(2,0),family=ZIP(log))
fit2$coef
fit2$aic;fit2$sbc
checkresiduals(fit2)
acf2(fit2$residuals)
qqPlot(fit2$residuals)
shapiro.test(fit2$residuals)

round(2*min(pnorm(fit2$coef[1]/sqrt(diag(fit2$vcov))[1]), 1-pnorm(fit2$coef[1]/sqrt(diag(fit2$vcov))[1])),4)
round(2*min(pnorm(fit2$coef[2]/sqrt(diag(fit2$vcov))[2]), 1-pnorm(fit2$coef[2]/sqrt(diag(fit2$vcov))[2])),4)
round(2*min(pnorm(fit2$coef[3]/sqrt(diag(fit2$vcov))[3]), 1-pnorm(fit2$coef[3]/sqrt(diag(fit2$vcov))[3])),4)
round(2*min(pnorm(fit2$coef[4]/sqrt(diag(fit2$vcov))[4]), 1-pnorm(fit2$coef[4]/sqrt(diag(fit2$vcov))[4])),4)

#Fitted values plot
plot(y,ty="l", col="blue",ylab="Monthly Counts",xlab="Time",ylim=c(0,3))
lines(fitt$mu,col="red",lty=2)
legend(92, 3, legend=c("Observed", "Fitted"),
       col=c("blue", "red"), lty=c(1,2),bty="n")

#Forecast twelve steps-ahead
y_for=forecast_bergARMA(0.0,0.0,c(0.257,0.368),0,0.778,y,X=cbind(rep(1,n),1),X_for=cbind(rep(1,12),1),h=12)
plot(y_for,type="l",ylim=c(-1,7),ylab="Monthly Counts",xlab="h")
points(y_for,pch=16,col="red")
points(ycomp[121:132],pch=1,col="blue")
lines(qberg(0.025,y_for,0.778),col="red",lty=2)
lines(qberg(0.975,y_for,0.778),col="red",lty=2)
legend(10, 7, legend=c("Observed", "Predict"),
       col=c("blue","red"), pch=c(1,16),bty="n")
legend(10, 6, legend=c("95% CI"),
       col=c("red"), lty=2,bty="n")

#MSE
mean((y_for-ycomp[121:132])^2)
#RSME
sqrt(mean((y_for-ycomp[121:132])^2))


#Forecast one step-ahead
n_y_for=NULL
for(i in 1:12){
  n_y_for[i]=forecast_bergARMA(0.0,0.0,c(0.257,0.368),0,0.778,
         ycomp[1:(120+i)],X=cbind(rep(1,(n+i)),1),X_for=cbind(rep(1,12),1)[i,],h=1)
}

n_y_for
#MSE
mean((n_y_for-ycomp[121:132])^2)
#RSME
sqrt(mean((n_y_for-ycomp[121:132])^2))

plot(n_y_for,type="l",ylim=c(-1,7),ylab="Monthly Counts",xlab="h")
points(n_y_for,pch=16,col="red")
points(ycomp[121:132],pch=1,col="blue")
lines(qberg(0.025,n_y_for,0.778),col="red",lty=2)
lines(qberg(0.975,n_y_for,0.778),col="red",lty=2)
legend(10, 7, legend=c("Observed", "Predict"),
       col=c("blue","red"), pch=c(1,16),bty="n")
legend(10, 6, legend=c("95% CI"),
       col=c("red"), lty=2,bty="n")













############ OVERDISPERSED APLICATION

##### Contagens fila contabilidade

ycomp <- scan("fila_CONT.txt") #Complete dataset
y <- ycomp[1:120]
var(y)/mean(y)
ts.plot(y,ylab="Monthly Counts")
ts.plot(diff(y))
acf2(y,max.lag = 10)
ggsubseriesplot(x=ts(y,frequency=12),ylab="Counts")
n=length(y)
X=rep(0,times=120)
X[seq(7,120,12)]=1 #mês de julho

ufit=glm.bg(y~X)
ufit$coefficients
acf2(ufit$resid,max.lag = 10,main="")
qqPlot(ufit$resid)
shapiro.test(ufit$resid)


X_new=cbind(1,X)
beta=c(0.1,0.2);phi=c(0.24);theta=c(0.0);nu=1
fitt=bergARMA(beta,phi,theta,nu,y,X_new,link="log")
est=round(fitt$par,4);est
fitt$AIC;fitt$BIC 
sd_est=diag(fitt$cov)^(1/2);round(sd_est,4)

#Confidence interval
est[1]+c(-1,1)*sd_est[1]*1.96
est[2]+c(-1,1)*sd_est[2]*1.96
est[3]+c(-1,1)*sd_est[3]*1.96
est[4]+c(-1,1)*sd_est[4]*1.96

#Hypothesis test
round(2*min(pnorm(est[1]/sd_est[1]), 1-pnorm(est[1]/sd_est[1])),3)
round(2*min(pnorm(est[2]/sd_est[2]), 1-pnorm(est[2]/sd_est[2])),3)
round(2*min(pnorm(est[3]/sd_est[3]), 1-pnorm(est[3]/sd_est[3])),3)
round(2*min(pnorm((est[4]-1)/sd_est[4]), 1-pnorm((est[4]-1)/sd_est[4])),3)


table(fitt$par[length(fitt$par)]>abs(fitt$mu-1))

#Residuals chart
var(fitt$resid);mean(fitt$res)
plot(fitt$resid,ylab="Quantile Residuals")
abline(h=-3)
aa=seq(-5,n+6,by=1)
par(mfrow=c(1,1))
par(mar=c(2.8, 2.7, 1.2, 1)) # margens c(baixo,esq,cima,direia)
par(mgp=c(1.7, 0.45, 0))
plot(fitt$resid,main=" ",xlab="Index",ylab="Quantile residuals", pch = "+",ylim=c(-4,4))
lines(aa,rep(-3,n+12),lty=2,col=1)
lines(aa,rep(3,n+12),lty=2,col=1)
lines(aa,rep(-2,n+12),lty=3,col=1)
lines(aa,rep(2,n+12),lty=3,col=1)


acf2(fitt$resid,main="",max.lag = 15)
qqPlot(fitt$resid,ylab="Quantile residuals")
shapiro.test(fitt$resid)
Box.test(fitt$resid,lag=15,fitdf=1)

#Fitting the other models
#Negative binomial
fit1=garmaFit(y~X,order=c(1,0),family=NBI(log))
fit1$coef
fit1$aic;fit1$sbc
round(diag(fit1$vcov)^(1/2),4)
checkresiduals(fit1)
acf2(fit1$residuals)
qqPlot(fit1$residuals)
shapiro.test(fit1$residuals)
round(2*min(pnorm(fit1$coef[1]/sqrt(diag(fit1$vcov))[1]), 1-pnorm(fit1$coef[1]/sqrt(diag(fit1$vcov))[1])),4)
round(2*min(pnorm(fit1$coef[2]/sqrt(diag(fit1$vcov))[2]), 1-pnorm(fit1$coef[2]/sqrt(diag(fit1$vcov))[2])),4)
round(2*min(pnorm(fit1$coef[3]/sqrt(diag(fit1$vcov))[3]), 1-pnorm(fit1$coef[3]/sqrt(diag(fit1$vcov))[3])),4)
round(2*min(pnorm(fit1$coef[4]/sqrt(diag(fit1$vcov))[4]), 1-pnorm(fit1$coef[4]/sqrt(diag(fit1$vcov))[4])),4)


#ZAP

fit2=garmaFit(y~X,order=c(1,0),family=ZAP(log))
fit2$coef 
round(diag(fit2$vcov)^(1/2),4)
fit2$aic;fit2$sbc

checkresiduals(fit2)
acf2(fit2$residuals)
qqPlot(fit2$residuals)
round(2*min(pnorm(fit2$coef[1]/sqrt(diag(fit2$vcov))[1]), 1-pnorm(fit2$coef[1]/sqrt(diag(fit2$vcov))[1])),4)
round(2*min(pnorm(fit2$coef[2]/sqrt(diag(fit2$vcov))[2]), 1-pnorm(fit2$coef[2]/sqrt(diag(fit2$vcov))[2])),4)
round(2*min(pnorm(fit2$coef[3]/sqrt(diag(fit2$vcov))[3]), 1-pnorm(fit2$coef[3]/sqrt(diag(fit2$vcov))[3])),4)
round(2*min(pnorm(fit2$coef[4]/sqrt(diag(fit2$vcov))[4]), 1-pnorm(fit2$coef[4]/sqrt(diag(fit2$vcov))[4])),4)

#Fitted values plot
plot(y,ty="l", col="blue",ylab="Monthly Counts",xlab="Time")
lines(fitt$mu,col="red",lty=2)
legend(92, 11, legend=c("Observed", "Fitted"),
      col=c("blue", "red"), lty=c(1,2),bty="n")

#Forecasting part
X_for=c(0,0,0,0,0,0,1,0,0,0,0,0)
X_for=cbind(1,X_for);X_for

y_for=forecast_bergARMA(0.9440,-1.0547,0.3112,0,1.9832,y,X=X_new,X_for=X_for,h=12)
plot(y_for,type="l",ylim=c(-1,14),ylab="Monthly Counts",xlab="h")
points(y_for,pch=16,col="red")
points(ycomp[121:132],pch=1,col="blue")
lines(qberg(0.025,y_for,1.9832),col="red",lty=2)
lines(qberg(0.975,y_for,1.9832),col="red",lty=2)
legend(10, 14, legend=c("Observed", "Predict"),
       col=c("blue","red"), pch=c(1,16),bty="n")
legend(10, 12, legend=c("95% CI"),
       col=c("red"), lty=2,bty="n")

#MSE
mean((y_for-ycomp[121:132])^2)
#RSME
sqrt(mean((y_for-ycomp[121:132])^2))


#Forecast one step-ahead
n_y_for=NULL
for(i in 1:12){
  n_y_for[i]=forecast_bergARMA(0.9440,-1.0547,0.3112,0,1.9832
                    ,ycomp[1:(120+i)],X=rbind(X_new,X_for)[1:(120+i),]
                    ,X_for=X_for[i,],h=1)
}

n_y_for
#MSE
mean((n_y_for-ycomp[121:132])^2)
#RSME
sqrt(mean((n_y_for-ycomp[121:132])^2))

plot(n_y_for,type="l",ylim=c(-1,14),ylab="Monthly Counts",xlab="h")
points(n_y_for,pch=16,col="red")
points(ycomp[121:132],pch=1,col="blue")
lines(qberg(0.025,n_y_for,1.9832),col="red",lty=2)
lines(qberg(0.975,n_y_for,1.9832),col="red",lty=2)
legend(10, 14, legend=c("Observed", "Predict"),
       col=c("blue","red"), pch=c(1,16),bty="n")
legend(10, 12, legend=c("95% CI"),
       col=c("red"), lty=2,bty="n")
