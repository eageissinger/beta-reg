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

# ----- Fish models -----
# explore data
head(cod)
summary(cod)

# General Linear Model
m1.fish<-glm(pi~factor(pulse)+year+month,data=na.omit(cod),
        family = gaussian(link="identity"))
plot(m1.fish)
hist(resid(m1.fish))
# plots for manuscript
plot(x=fitted(m1.fish),y=resid(m1.fish),main=NULL,xlab="Fitted Values",ylab="Residuals")
hist(resid(m1.fish),main=NULL,xlab="Residuals")
qqnorm(resid(m1.fish),main=NULL)
qqline(resid(m1.fish),col='red')
plot(m1.fish,which=4,main=NULL)

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
hist(resid(m2.fish))


#ANODEV
m2.intercept<-gamlss(pi~0,family=BEZI,data=na.omit(cod),trace=F)
m2.pulse<-gamlss(pi~1+factor(pulse),family=BEZI,data=na.omit(cod),trace=F)
m2.yr<-gamlss(pi~1+factor(pulse)+year,family=BEZI,data=na.omit(cod),trace=F)
m2.month<-gamlss(pi~1+factor(pulse)+year+month,family=BEZI,data=na.omit(cod),trace=F)

m2ANODEV<-lrtest(m2.intercept,m2.pulse,m2.yr,m2.month)
m2ANODEV
# ---- Model comparison: fish ----
glm.model.fish<-as.data.frame(m1ANODEV)%>%
  rename(numDf='#Df')

glm.fish<-glm.model.fish%>%
  mutate(Deviance=2*LogLik)%>%
  mutate(dDeviance=Deviance-lag(Deviance))%>%
  mutate(LR=exp(dDeviance/2))%>%
  mutate(AIC=dDeviance-2*numDf)%>%
  mutate(LR_Df=LR/Df)%>%
  mutate(dAIC=AIC-lag(AIC))

beta.model.fish<-as.data.frame(m2ANODEV)%>%
  rename(numDf='#Df')
beta.fish<-beta.model.fish%>%
  mutate(Deviance=2*LogLik)%>%
  mutate(dDeviance=Deviance-lag(Deviance))%>%
  mutate(LR=exp(dDeviance/2))%>%
  mutate(AIC=dDeviance-2*numDf)%>%
  mutate(LR_Df=LR/Df)%>%
  mutate(dAIC=AIC-lag(AIC))

# save tables
write.csv(glm.fish,"./output/glm.model.csv",row.names=FALSE)
write.csv(beta.fish,"./output/beta.model.csv",row.names = FALSE)

# ---- Model parameters: fish ----
m1.fish$coefficients
summary(m1.fish)
m1.fish$coefficients

summary(m2.fish)
betacoef<-m2.fish$mu.coefficients
odds<-exp(betacoef)
beta.prob<-odds/(1+odds)

parameter.summary<-data.frame(glm=m1.fish$coefficients,
                              beta=beta.prob)
# save table
write.csv(parameter.summary,"./output/parameters.csv",row.names = FALSE)
# ---- Fit statistics: fish -----
fit.fish<-data.frame(Fit.statistics=c('AIC','LR','Residual Deviance'),
                     glm=c(AIC(m1.fish),exp(logLik(m1.fish)),
                           sum(m1.fish$residuals^2)/m1.fish$df.residual),
                     beta=c(AIC(m2.fish),exp(logLik(m2.fish)),
                            sum(m2.fish$residuals^2)/m2.fish$df.residual))
# save table
write.csv(fit.fish,"./output/fit.stats.csv",row.names = FALSE)


# ---- Regular Beta Regression -----
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
m2<-betareg(pi~factor(pulse),data=cod0)
plot(m2)
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
