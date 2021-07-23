# Beta Regression Analysis 

# set working directory
setwd("~/Desktop/Beta Regression Paper/Beta Regression R analysis")

# load data
SoilNMR<-read.csv("NMR data.csv")

# load packages
library(betareg)
library(car)
library(gamlss)
library(emmeans)
library(dplyr)
library(lmtest)
library(glmmTMB)
library(DHARMa)

------------------------------------CARBONYL CARBON ANALYSES------------------------------------

# ----- Soil NMR -----
# explore data
head(SoilNMR)
summary(SoilNMR)

# General Linear Model --- Carbonyl Carbon
m1.SoilNMR<-glm(CC~Region+Layer+Region*Layer,data=na.omit(SoilNMR),
             family = gaussian(link="identity"))
plot(m1.SoilNMR)
hist(resid(m1.SoilNMR))
# plots for manuscript
plot(x=fitted(m1.SoilNMR),y=resid(m1.SoilNMR),main=NULL,xlab="Fitted Values",ylab="Residuals")
hist(resid(m1.SoilNMR),main=NULL,xlab="Residuals")
qqnorm(resid(m1.SoilNMR),main=NULL)
qqline(resid(m1.SoilNMR),col='red')
plot(m1.SoilNMR,which=4,main=NULL)
# model summary/diagnostics
summary(m1.SoilNMR)
exp(logLik(m1.SoilNMR))
anova(m1.SoilNMR)
#ANODEV
m1.intercept<-glm(CC~1,data=na.omit(SoilNMR),
                  family=gaussian(link = "identity"))
m1.OHCC<-glm(CC~1+Region,data=na.omit(SoilNMR),
              family=gaussian(link = "identity"))
m1.CR<-glm(CC~1+Region+Layer,data=na.omit(SoilNMR),
           family=gaussian(link = "identity"))
m1.OHCCCR<-glm(CC~1+Region+Layer+Region*Layer,data=na.omit(SoilNMR),
              family=gaussian(link = "identity"))
m1ANODEV<-lrtest(m1.intercept,m1.OHCC,m1.CR,m1.OHCCCR)
m1ANODEV

# Regular beta regressions {betareg}
library(betareg)
m3<-betareg(CC~Region+Layer+Region*Layer,data=SoilNMR)
plot(m3)
# diagnostic plots for manuscript
plot(m3,which=1,type="pearson")
plot(m3,which=4,type="pearson")
plot(m3,which=5,type="deviance")
plot(m3,which=2,type="pearson")
# model summary
summary(m3)
#ANODEV
m3.intercept<-betareg(CC~1,data=SoilNMR)
m3.OHCC<-betareg(CC~1+Region,data=SoilNMR)
m3.CR<-betareg(CC~1+Region+Layer,data=SoilNMR)
m3.OHCCCR<-betareg(CC~1+Region+Layer+Layer*Region,data=SoilNMR)

install.packages("lrtest")
m3ANODEV<-lrtest(m3.intercept,m3.OHCC,m3.CR,m3.OHCCCR)
m3ANODEV # same result as the GLMMTMB

# ---- Model comparison: SoilNMR ----

glm.model.SoilNMR.CC<-as.data.frame(m1ANODEV)%>%
  rename(numDf='#Df')

CCglm.SoilNMR<-glm.model.SoilNMR.CC%>%
  mutate(Deviance=2*LogLik)%>%
  mutate(dDeviance=Deviance-lag(Deviance))%>%
  mutate(LR=exp(dDeviance/2))%>%
  mutate(AIC=dDeviance-2*numDf)%>%
  mutate(LR_Df=LR/Df)%>%
  mutate(dAIC=AIC-lag(AIC))

beta.model.SoilNMR.CC<-as.data.frame(m3ANODEV)%>%
  rename(numDf='#Df')

CCbeta.SoilNMR<-beta.model.SoilNMR.CC%>%
  mutate(Deviance=2*LogLik)%>%
  mutate(dDeviance=Deviance-lag(Deviance))%>%
  mutate(LR=exp(dDeviance/2))%>%
  mutate(AIC=dDeviance-2*numDf)%>%
  mutate(LR_Df=LR/Df)%>%
  mutate(dAIC=AIC-lag(AIC))

# save tables
write.csv(CCglm.SoilNMR,file ="CCglm.model.csv",row.names=FALSE)
write.csv(CCbeta.SoilNMR,file ="CCbeta.model.csv",row.names = FALSE)

# ---- Model parameters: SoilNMR ----
m1.SoilNMR$coefficients
summary(m1.SoilNMR)
m1.SoilNMR$coefficients

summary(m3)
betacoef<-m3$mu.coefficients
odds<-exp(coef(m3))
beta.prob<-odds/(1+odds)

parameter.summary.CC<-data.frame(glm=m1.SoilNMR$coefficients,
                              beta=beta.prob[-c(13)])

# save table
write.csv(parameter.summary.CC,"CCparameters.csv",row.names = FALSE)

# ---- Fit statistics: SoilNMR -----
fit.SoilNMR<-data.frame(Fit.statistics=c('AIC','LR','Residual Deviance'),
                     glm=c(AIC(m1.SoilNMR),exp(logLik(m1.SoilNMR)),
                           sum(m1.SoilNMR$residuals^2)/m1.SoilNMR$df.residual),
                     beta=c(AIC(m3),exp(logLik(m3)),
                            sum(m3$residuals^2)/m3$df.residual))

# save table
write.csv(fit.SoilNMR,"Carbonylcarbonfit.stats.csv",row.names = FALSE)
