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
-------------------------------------ALKYL CARBON ANALYSES------------------------------------

  # ----- Soil NMR -----
# explore data
head(SoilNMR)
summary(SoilNMR)

# General Linear Model --- Alkyl Carbon
m4.SoilNMR<-glm(AC~Region+Layer+Region*Layer,data=na.omit(SoilNMR),
                family = gaussian(link="identity"))
plot(m4.SoilNMR)
hist(resid(m4.SoilNMR))
# plots for manuscript
plot(x=fitted(m4.SoilNMR),y=resid(m4.SoilNMR),main=NULL,xlab="Fitted Values",ylab="Residuals")
hist(resid(m4.SoilNMR),main=NULL,xlab="Residuals")
qqnorm(resid(m4.SoilNMR),main=NULL)
qqline(resid(m4.SoilNMR),col='red')
plot(m4.SoilNMR,which=4,main=NULL)
# model summary/diagnostics
summary(m4.SoilNMR)
exp(logLik(m4.SoilNMR))
anova(m4.SoilNMR)
#ANODEV
m4.intercept<-glm(AC~1,data=na.omit(SoilNMR),
                  family=gaussian(link = "identity"))
m4.OHCC<-glm(AC~1+Region,data=na.omit(SoilNMR),
             family=gaussian(link = "identity"))
m4.CR<-glm(AC~1+Region+Layer,data=na.omit(SoilNMR),
           family=gaussian(link = "identity"))
m4.OHCCCR<-glm(AC~1+Region+Layer+Region*Layer,data=na.omit(SoilNMR),
               family=gaussian(link = "identity"))
m4ANODEV<-lrtest(m4.intercept,m4.OHCC,m4.CR,m4.OHCCCR)
m4ANODEV

# Regular beta regressions {betareg}
library(betareg)
m6.SoilNMR<-betareg(AC~Region+Layer+Region*Layer,data=SoilNMR)
plot(m6.SoilNMR)
# diagnostic plots for manuscript
plot(m6.SoilNMR,which=1,type="pearson")
plot(m6.SoilNMR,which=4,type="pearson")
plot(m6.SoilNMR,which=5,type="deviance")
plot(m6.SoilNMR,which=2,type="pearson")
# model summary
summary(m6.SoilNMR)
#ANODEV
m6.intercept<-betareg(AC~1,data=SoilNMR)
m6.OHCC<-betareg(AC~1+Region,data=SoilNMR)
m6.CR<-betareg(AC~1+Region+Layer,data=SoilNMR)
m6.OHCCCR<-betareg(AC~1+Region+Layer+Region*Layer,data=SoilNMR)

install.packages("lrtest")
m6ANODEV<-lrtest(m6.intercept,m6.OHCC,m6.CR,m6.OHCCCR)
m6ANODEV # same result as the GLMMTMB


# ---- Model comparison: SoilNMR ----
glm.model.SoilNMR.AC<-as.data.frame(m4ANODEV)%>%
  rename(numDf='#Df')

glm.SoilNMR<-glm.model.SoilNMR.AC%>%
  mutate(Deviance=2*LogLik)%>%
  mutate(dDeviance=Deviance-lag(Deviance))%>%
  mutate(LR=exp(dDeviance/2))%>%
  mutate(AIC=dDeviance-2*numDf)%>%
  mutate(LR_Df=LR/Df)%>%
  mutate(dAIC=AIC-lag(AIC))

beta.model.SoilNMR.AC<-as.data.frame(m6ANODEV)%>%
  rename(numDf='#Df')
beta.SoilNMR<-beta.model.SoilNMR.AC%>%
  mutate(Deviance=2*LogLik)%>%
  mutate(dDeviance=Deviance-lag(Deviance))%>%
  mutate(LR=exp(dDeviance/2))%>%
  mutate(AIC=dDeviance-2*numDf)%>%
  mutate(LR_Df=LR/Df)%>%
  mutate(dAIC=AIC-lag(AIC))

# save tables
write.csv(glm.SoilNMR,file ="AC.glm.model.csv",row.names=FALSE)
write.csv(beta.SoilNMR,file ="AC.beta.model.csv",row.names = FALSE)

# ---- Model parameters: SoilNMR ----
m4.SoilNMR$coefficients
summary(m4.SoilNMR)
m4.SoilNMR$coefficients

summary(m6.SoilNMR)
betacoef.AC<-m6.SoilNMR$mu.coefficients
odds.AC<-exp(coef(m6.SoilNMR))
beta.prob.AC<-odds/(1+odds)

parameter.summary.AC<-data.frame(glm=m4.SoilNMR$coefficients,
                              beta=beta.prob.AC[-c(13)])

# save table
write.csv(parameter.summary,"AC.parameters.csv",row.names = FALSE)

# ---- Fit statistics: SoilNMR -----
fit.SoilNMR<-data.frame(Fit.statistics=c('AIC','LR','Residual Deviance'),
                        glm=c(AIC(m4.SoilNMR),exp(logLik(m4.SoilNMR)),
                              sum(m4.SoilNMR$residuals^2)/m4.SoilNMR$df.residual),
                        beta=c(AIC(m6.SoilNMR),exp(logLik(m6.SoilNMR)),
                               sum(m6.SoilNMR$residuals^2)/m6.SoilNMR$df.residual))

# save table
write.csv(fit.SoilNMR,"Alkylcarbonfit.stats.csv",row.names = FALSE)



