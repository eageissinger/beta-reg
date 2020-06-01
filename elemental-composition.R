# Author: Isabella Richmond & Emilie Geissinger (https://github.com/eageissinger)
# Date: June 1, 2020
# code for beta regression in the natural sciences paper. Comparing beta regression to
# general linear model. Using %N for Abies balsamea (balsam fir) 0 modelling site and year.
# comparing models with likelihood ratios, residual deviance/residual, and AICc.

#Load required packages 
library("easypackages")
libraries("betareg", "car", "lme4", "lmtest", "DescTools", "tidyverse")

# load data
stoich <- read_csv("C:/Users/Isabella Richmond/Documents/M.Sc/General Linear Models/Final Exam/Data/TotalStoich_2016_2017_Biomass.csv")
# subset balsam fir (ABBA)
abba <- subset(stoich, Species == "ABBA")
head(abba)
summary(abba)

# ---- Elemental composition models ----
# GzLM with normal error and identity link 
# nitrogen (%)
gzlmN <- glm(N_dec~Year+Site, family = gaussian(link = identity), data = abba)
gzlmNFrameNormal<-cbind(abba,residuals(gzlmN),fitted(gzlmN))
# diagnostic plots for manuscript 
png("C:/Users/Isabella Richmond/Documents/M.Sc/General Linear Models/Beta Regression/gzlmN_diagnostics.png", width = 800, height = 650)
par(mfrow=c(2,2))
plot(x=fitted(gzlmN),y=resid(gzlmN),main=NULL,xlab="Fitted Values",ylab="Residuals")
hist(resid(gzlmN),main=NULL,xlab="Residuals")
qqnorm(resid(gzlmN),main=NULL)
qqline(resid(gzlmN),col='red')
plot(gzlmN,which=4,main=NULL)
dev.off()
# model summary and fit statistics 
summary(gzlmN)
exp(logLik(gzlmN))
anova(gzlmN)
# build ANODEV table sequentially, same method as used for beta regression 
gzlmNintercept<-glm(N_dec~1, family = gaussian(link = identity),data=abba)
gzlmNyear<-glm(N_dec~1+Year, family = gaussian(link = identity),data=abba)
gzlmNyearsite<-glm(N_dec~1+Year+Site, family = gaussian(link = identity),data=abba)
gzlmNANODEV <- lrtest(gzlmNintercept,gzlmNyear,gzlmNyearsite)
gzlmNANODEV

# now use beta regression for the same model 
# nitrogen (%)
brN<-betareg(N_dec~Year+Site, link="logit",data=abba)
brNFrameBeta<-cbind(abba,residuals(brN),fitted(brN))
# diagnostic plots for manuscript
png("C:/Users/Isabella Richmond/Documents/M.Sc/General Linear Models/Beta Regression/betaN_diagnostics.png", width = 800, height = 650)
par(mfrow=c(2,2))
plot(brN,which=1,type="pearson")
plot(brN,which=4,type="pearson")
plot(brN,which=5,type="deviance")
plot(brN,which=2,type="pearson")
dev.off()
# model summary and fit statistics
summary(brN)
AIC(brN)
exp(logLik(brN))
res.devN<-residuals(brN,type = "deviance")
(sum(res.devN^2))/109
# ANODEV
brNintercept<-betareg(N_dec~1, link="logit",data=abba)
brNyear<-betareg(N_dec~1+Year, link="logit",data=abba)
brNyearsite<-betareg(N_dec~1+Year+Site, link="logit",data=abba)
brNANODEV <- lrtest(brNintercept,brNyear,brNyearsite)
brNANODEV #gives the same result as GLMMTMB


# ---- Model comparison: elemental composition ----
# ANOVA Table GLM
glm.model.abbaN<-as.data.frame(gzlmNANODEV)%>%
  rename(numDf='#Df') # adjust name of #DF so that R can run
# ANOVA with parameter comparison for GLM
glm.abbaN<-glm.model.abbaN%>%
  mutate(Deviance=2*LogLik)%>%
  mutate(dDeviance=Deviance-lag(Deviance))%>%
  mutate(LR=exp(dDeviance/2))%>%
  mutate(AIC=dDeviance-2*numDf)%>%
  mutate(LR_Df=LR/Df)%>%
  mutate(dAIC=AIC-lag(AIC))

# ANOVA Table Beta
beta.model.abbaN<-as.data.frame(brNANODEV)%>%
  rename(numDf='#Df') # adjust name of #DF so that R can run
# ANOVA with parameter comparison for beta
beta.abbaN<-beta.model.abbaN%>%
  mutate(Deviance=2*LogLik)%>%
  mutate(dDeviance=Deviance-lag(Deviance))%>%
  mutate(LR=exp(dDeviance/2))%>%
  mutate(AIC=dDeviance-2*numDf)%>%
  mutate(LR_Df=LR/Df)%>%
  mutate(dAIC=AIC-lag(AIC))

# save tables
write.csv(glm.abbaN,"C:/Users/Isabella Richmond/Documents/M.Sc/General Linear Models/Beta Regression/glm.model.csv",row.names=FALSE)
write.csv(beta.abbaN,"C:/Users/Isabella Richmond/Documents/M.Sc/General Linear Models/Beta Regression/beta.model.csv",row.names = FALSE)

# ---- Model parameters: elemental composition ----
# calculate model parameters for glm and beta
gzlmN$coefficients
summary(gzlmN)

summary(brN)
betacoef<-brN$coefficients
odds<-exp(betacoef$mean)
beta.prob<-odds/(1+odds)

# combine model parameters into single table for comparison
parameter.summary<-data.frame(glm=gzlmN$coefficients,
                              beta=beta.prob)
# save table
write.csv(parameter.summary,"C:/Users/Isabella Richmond/Documents/M.Sc/General Linear Models/Beta Regression/parameters.csv",row.names = FALSE)

# ---- Fit statistics: elemental composition -----

# Calculate r-squared for LR
# LR = (1-R^2)^(-n/2), where n is the number of observations
# Use Cox and Snell method (Cox and Snell 1989)
# GLM LR
PseudoR2(gzlmN, which = "CoxSnell") # r-squared for glm, extract cox and snell
(1-PseudoR2(gzlmN, which = "CoxSnell"))^(-101/2) 
# BetaReg LR
(1-brN$pseudo.r.squared)^(-101/2) # LR for {betareg}

fit.abbaN<-data.frame(Fit.statistics=c('AIC','LR','Residual Deviance'),
                     glm=c(AIC(gzlmN),(1-PseudoR2(gzlmN, which = "CoxSnell"))^(-101/2),
                           sum(gzlmN$residuals^2)/gzlmN$df.residual),
                     beta=c(AIC(brN),(1-brN$pseudo.r.squared)^(-101/2),
                            sum(brN$residuals^2)/brN$df.residual))
# save table
write.csv(fit.abbaN,"C:/Users/Isabella Richmond/Documents/M.Sc/General Linear Models/Beta Regression/fit.stats.csv",row.names = FALSE)
