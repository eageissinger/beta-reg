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
library(DHARMa) # for glmmTMB
library(rcompanion)
library(hnp)
library(ggpubr)
library(grDevices)
library(gamlssdiag)

source("revised-gamlss-plot-fx.R")
# explore data
head(cod)
summary(cod)

# General Linear Model
m1.fish<-glm(pi~factor(pulse)+year+month,data=na.omit(cod),
        family = gaussian(link="identity"))
plot(m1.fish)
hist(resid(m1.fish),cex.lab=1.25)
mtext("A",side=2,line=1,at=50,col='black',las=1,font=2,size=2)
hnp(m1.fish,pch=16) #half normal with envelope (for dispersion)
# plots for manuscript
# plots for manuscript
png(width = 160, height = 160, units = "mm",res =600)
par(mfrow=c(2,2))
plot(x=fitted(m1.fish),y=resid(m1.fish),main=NULL,
     xlab="Fitted Values",ylab="Residuals",cex.lab=1.15)
mtext("A",side=2,line=2,at=1,col='black',font=2,las=1,size=1.75)
hist(resid(m1.fish),main=NULL,xlab="Residuals",
     cex.lab=1.15)
mtext("B",side=2,line=2,at=50,col='black',font=2,las=1,size=1.75)
qqnorm(resid(m1.fish),main=NULL,cex.lab=1.15)
mtext("C",side=2,line=2,at=1.1,col='black',font=2,las=1,size=1.75)
qqline(resid(m1.fish),col='red')
plot(m1.fish,which=4,caption = NULL,main = NULL,cex.lab=1.15) 
mtext("D",side=2,line=2,at=.26,col='black',font=2,las=1,size=1.75)
dev.off()
# will need to white out the equation
# save as 800 x 650 PNG

# ---- Manuscript Figures ----
fit<-fitted(m1.fish)
res<-resid(m1.fish)

residuals<-data.frame(cbind(fit,res))

fig1a<-ggplot(residuals,aes(x=fit,y=res))+
  geom_point(fill='white')+
  theme_classic()+
  ylab("Residuals")+xlab("Fitted Values")
fig1b<-ggplot(residuals)+geom_histogram(aes(x=res),bins=10,colour='black',fill='white')+
  theme_classic()+
  ylab("Frequency")+xlab("Residuals")
fig1data=data.frame(qqnorm(resid(m1.fish),main=NULL))
fig1c<-ggplot(fig1data)+geom_point(aes(x=x,y=y))+
  geom_smooth(method="lm",aes(x=x,y=y),se=FALSE,color='red')+
  theme_classic()+
  xlab("Theoretical Quantiles")+ylab("Sample Quantiles")
fig1d<-ggplot(m1.fish, aes(seq_along(.cooksd), .cooksd))+geom_bar(stat="identity", position="identity",colour='black',fill='black')+
  xlab("Obs. Number")+ylab("Cook's distance")+
  theme_classic()

Fig1<-ggarrange(fig1a,fig1b,fig1c,fig1d,nrow=2,ncol=2,
                labels=c("a","b","c","d"))
hnp(m1.fish,pch=16)->test
x<-test["x"]
u<-test["upper"]
l<-test["lower"]
m<-test["median"]
r<-test["residuals"]
envelope<-bind_cols(x,u,l,m,r)
ggplot(envelope)+geom_point(aes(x=x,y=residuals))+
  geom_line(aes(x=x,y=upper))+
  geom_line(aes(x=x,y=lower))+
  geom_line(aes(x=x,y=median),linetype='dashed')+
  theme_classic()+
  xlab("Theoretical quantiles")+
  ylab("Residuals")
glimpse(envelope)

# model summary/diagnostics
summary(m1.fish)
exp(logLik(m1.fish))
anova(m1.fish)

#ANODEV
m1.intercept<-glm(pi~1,data=na.omit(cod),
                  family=gaussian(link = "identity"))
m1.pulse<-glm(pi~factor(pulse),data=na.omit(cod),
              family=gaussian(link = "identity"))
m1.yr<-glm(pi~factor(pulse)+year,data=na.omit(cod),
           family=gaussian(link = "identity"))
m1.month<-glm(pi~factor(pulse)+year+month,data=na.omit(cod),
              family=gaussian(link = "identity"))
m1ANODEV<-lrtest(m1.intercept,m1.pulse,m1.yr,m1.month)
m1ANODEV


# Zero-inflated Beta Regression
m2.fish<-gamlss(pi~factor(pulse)+year+month,family=BEZI,data=na.omit(cod),trace=F)
summary(m2.fish)
#model diagnostics
png(width = 160, height = 160, units = "mm",res =600)
plot(m2.fish)
dev.off()
hist(resid(m2.fish))
plot(cooksd(gamlss,pi~factor(pulse)+year+month,family=BEZI,data=na.omit(cod)),type="h",
     ylab="Cook's distance",xlab="Obs. Number")
dev.off()

#ANODEV
m2.intercept<-gamlss(pi~0,family=BEZI,data=na.omit(cod),trace=F)
m2.pulse<-gamlss(pi~1+factor(pulse),family=BEZI,data=na.omit(cod),trace=F)
m2.yr<-gamlss(pi~1+factor(pulse)+year,family=BEZI,data=na.omit(cod),trace=F)
m2.month<-gamlss(pi~1+factor(pulse)+year+month,family=BEZI,data=na.omit(cod),trace=F)

m2ANODEV<-lrtest(m2.intercept,m2.pulse,m2.yr,m2.month)
m2ANODEV

# ---- Regular Beta Regression -----
# using {betareg}
cod0<-cod%>%
  filter(pi!=0) # take out 0's for example set

m1<-betareg(pi~as.factor(pulse),data=cod0)
summary(m1)
plot(m1)
# diagnostic plots for manuscript
par(mfrow=c(2,2))
plot(m1,which=1,type="pearson",caption=NULL,pch=16)
plot(m1,which=4,type="pearson",caption=NULL)
plot(m1,which=5,type="deviance",caption = NULL)
plot(m1,which=2,type="pearson",caption = NULL)
# save as 800 x 650 PNG
par(mfrow=c(1,1))

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


# ---- Model comparison: fish ----
# ANOVA Table GLM
glm.model.fish<-as.data.frame(m1ANODEV)%>%
  rename(numDf='#Df') # adjust name of #DF so that R can run
# ANOVA with parameter comparison for GLM
glm.fish<-glm.model.fish%>%
  mutate(Deviance=2*LogLik)%>%
  mutate(dDeviance=Deviance-lag(Deviance))%>%
  mutate(LR=exp(dDeviance/2))%>%
  mutate(AIC=dDeviance-2*numDf)%>%
  mutate(LR_Df=LR/Df)%>%
  mutate(dAIC=AIC-lag(AIC))

# ANOVA Table Beta
beta.model.fish<-as.data.frame(m2ANODEV)%>%
  rename(numDf='#Df') # adjust name of #DF so that R can run
# ANOVA with parameter comparison for beta
beta.fish<-beta.model.fish%>%
  mutate(Deviance=2*LogLik)%>%
  mutate(dDeviance=Deviance-lag(Deviance))%>%
  mutate(LR=exp(dDeviance/2))%>%
  mutate(AIC=dDeviance-2*numDf)%>%
  mutate(dAIC=AIC-lag(AIC))

# save tables
write.csv(glm.fish,"./output/glm.model.csv",row.names=FALSE)
write.csv(beta.fish,"./output/beta.model.csv",row.names = FALSE)

# ---- Model parameters: fish ----
# Calculate model parameters for glm and beta
m1.fish$coefficients
summary(m1.fish)
m1.fish$coefficients

summary(m2.fish)
betacoef<-m2.fish$mu.coefficients
odds<-exp(betacoef)
beta.prob<-odds/(1+odds)

# combine model parameters into single table for comparison
parameter.summary<-data.frame(glm=m1.fish$coefficients,
                              beta=beta.prob)
# save table
write.csv(parameter.summary,"./output/parameters.csv",row.names = FALSE)
# ---- Fit statistics: fish -----

# Calculate r-squared for LR
# LR = (1-R^2)^(-n/2), where n is the number of observations
# Use Cox and Snell method (Cox and Snell 1989)
# GLM LR
nagelkerke(m1.fish) # r-squared for glm, extract cox and snell
(1-0.350309)^(-101/2) # cox snell
# BetaReg LR
(1-m1$pseudo.r.squared)^(-101/2) # LR for {betareg}
# zero-inflated beta LR
(1-Rsq(m2.fish,type=c("Cox Snell")))^(-101/2) # cox snell for beta regression

fit.fish<-data.frame(Fit.statistics=c('AIC','LR','Residual Deviance'),
                     glm=c(AIC(m1.fish),(1-0.350309)^(-101/2),
                           sum(m1.fish$residuals^2)/m1.fish$df.residual),
                     beta=c(AIC(m2.fish),(1-Rsq(m2.fish,type=c("Cox Snell")))^(-101/2),
                            sum(m2.fish$residuals^2)/m2.fish$df.residual))

# Alternative R-squared using Nagelkerke (Nagelkerke 1991) 
(1-Rsq(m2.fish))^(-101/2) # nagelkerke
nagelkerke(m1.fish)
(1-.818843)^(-101/2) # nagelkerke


# save table
write.csv(fit.fish,"./output/fit.stats.csv",row.names = FALSE)

# LR = exp(G/2)
