# Beta Regression Analysis 

# set working directory
setwd("C:/Users/user/Documents/Research/Beta/")

# load data
cod<-read.csv("beta-fish.csv")

# load packages
library(betareg)
library(car)
library(gamlss)
library(emmeans)
library(dplyr)
library(lmtest)
library(glmmTMB)
library(DHARMa)

# ----- Fish -----
# explore data
head(cod)
summary(cod)

# General Linear Model
m1.fish<-glm(pi~factor(pulse)+year+month,data=na.omit(cod),
        family = gaussian(link="identity"))
plot(m1.fish)
hist(resid(m1))
# plots for manuscript
plot(x=fitted(m1),y=resid(m1),main=NULL,xlab="Fitted Values",ylab="Residuals")
hist(resid(m1),main=NULL,xlab="Residuals")
qqnorm(resid(m1),main=NULL)
qqline(resid(m1),col='red')
plot(m1,which=4,main=NULL)
# model summary/diagnostics
summary(m1.fish)
exp(logLik(m1.fish))
anova(m1.fish)
#ANODEV
m1.intercept<-glm(pi~1,data=na.omit(cod),
                  family=gaussian(link = "identity"))
m1.pulse<-glm(pi~1+factor(pulse),data=na.omit(cod),
              family=gaussian(link = "identity"))
m1.yr<-glm(pi~1+factor(pulse)+year,data=na.omit(cod),
           family=gaussian(link = "identity"))
m1.month<-glm(pi~1+factor(pulse)+year+month,data=na.omit(cod),
              family=gaussian(link = "identity"))
m1ANODEV<-lrtest(m1.intercept,m1.pulse,m1.yr,m1.month)
m1ANODEV


# Zero-inflated Beta Regression
m2.fish<-gamlss(pi~factor(pulse)+year+month,family=BEZI,data=na.omit(cod),trace=F)
summary(m2.fish)
#model diagnostics
plot(m2.fish)

exp(logLik(m2.fish))

# ANODEV
m2.intercept<-gamlss(pi~1,family=BEZI,data=na.omit(cod0),trace=F)
m2.pulse<-gamlss(pi~1+factor(pulse),family=BEZI,data=na.omit(cod0),trace=F)
m2.yr<-gamlss(pi~1+factor(pulse)+year,family=BEZI,data=na.omit(cod0),trace=F)
m2.month<-gamlss(pi~1+factor(pulse)+year+month,family=BEZI,data=na.omit(cod0),trace=F)

m2ANODEV<-lrtest(m2.intercept,m2.pulse,m2.yr,m2.month)
m2ANODEV

# ---- Regular Beta Regression
# glmmTMB with type III ANOVA 
cod0<-cod%>%
  filter(pi!=0)
m3<-glmmTMB(pi~factor(pulse),
            beta_family(link = "logit"),
            data=cod0)
res3<-simulateResiduals(m3)
plot(res3)
summary(m3)
Anova.glmmTMB(m3,type = "III")

m3.intercept<-glmmTMB(pi~1,beta_family(link="logit"),data=cod0)
m3.pulse<-glmmTMB(pi~1+factor(pulse),beta_family(link="logit"),data=cod0)

m3ANODEV<-lrtest(m3.intercept,m3.pulse)
m3ANODEV # off from anova just a little bit

# Advantages: ANOVA table
# Disadvantages: no diagnostic plots

# Regular beta regressions {betareg}
m4<-betareg(pi~factor(pulse),data=cod0)
plot(m4)
# diagnostic plots for manuscript
plot(m2,which=1,type="pearson")
plot(m2,which=4,type="pearson")
plot(m2,which=5,type="deviance")
plot(m2,which=2,type="pearson")
# model summary
summary(m4)
#ANODEV
m4.intercept<-betareg(pi~1,data=cod0)
m4.pulse<-betareg(pi~1+factor(pulse),data=cod0)

m4ANODEV<-lrtest(m4.intercept,m4.pulse)
m4ANODEV # same result as the GLMMTMB

# might be able to use ANOVA table from GLMMTMB
# and use diagnostics from betareg
# Summaries produce the same estimates
# will use combo of glmmTMB and betareg for two analyses
# then use zero-inflated for the third
