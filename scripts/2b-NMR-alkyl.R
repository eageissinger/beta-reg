# Alkyl Carbon Analyses
# Comparing beta regression to general linear model with arcsin transformation. 
# Using alkyl carbon in soil - modelling region and soil layer
# comparing models with likelihood ratios, residual deviance/residual, and AIac.

#### Packages ####
# load required packages 
p <- c("betareg", "lmtest", "DescTools", "dplyr", "Metrics")
lapply(p, library, character.only = T)

#### Data #### 
# load data
soil <- read.csv("input/Beta_NMRdata.csv")
# select relevant columns
soilac <- select(soil, c(Region, Layer, AC))
# add carbonyl carbon data with arcsin square root transformation
soilac <- mutate(soilac, AC_arc = asin(sqrt(soilac$AC)))

#### General Linear Model with Arcsin Transformation ####
m1 <- glm(AC_arc ~ Region+Layer+Region*Layer,data=na.omit(soilac))
m1FrameNormal<-cbind(soilac,residuals(m1),fitted(m1))
# diagnostic plots for manuscript
png("figures/glm_alkyl-carbon_diagnostics.png", width = 160, height = 160, units = "mm",res = 600)
par(mfrow=c(2,2))
plot(x=fitted(m1),y=resid(m1),main=NULL,
     xlab="Fitted Values",ylab="Residuals",cex.lab=1.15)
mtext("A",side=2,line=2,at=0.06,col='black',font=2,las=1,size=1.75)
hist(resid(m1),main=NULL,xlab="Residuals",
     cex.lab=1.15)
mtext("B",side=2,line=2,at=9,col='black',font=2,las=1,size=1.75)
qqnorm(resid(m1),main=NULL,cex.lab=1.15)
mtext("C",side=2,line=2,at = 0.06,col='black',font=2,las=1,size=1.75)
qqline(resid(m1),col='red')
plot(m1,which=4,caption = NULL,main = NULL,cex.lab=1.15)
mtext("D",side=2,line=2,at=0.3,col='black',font=2,las=1,size=1.75)
dev.off()
# model summary and fit statistics
summary(m1)
exp(logLik(m1))
anova(m1)
#ANODEV
m1.intercept<-glm(AC_arc~1,data=na.omit(soilac))
m1.OHac<-glm(AC_arc~1+Region,data=na.omit(soilac))
m1.CR<-glm(AC_arc~1+Region+Layer,data=na.omit(soilac))
m1.OHacCR<-glm(AC_arc~1+Region+Layer+Region*Layer,data=na.omit(soilac))
m1ANODEV<-lrtest(m1.intercept,m1.OHac,m1.CR,m1.OHacCR)
m1ANODEV

#### Beta Regression ####
m2 <- betareg(AC~Region+Layer+Region*Layer,data=soilac)
m2FrameBeta<-cbind(soilac,residuals(m2),fitted(m2))
# diagnostic plots for manuscript
png("figures/beta_alkyl-carbon_diagnostics.png",  width = 160, height = 160, units = "mm",res = 600)
par(mfrow=c(2,2))
plot(x=fitted(m2),y=resid(m2),main=NULL,
     xlab="Fitted Values",ylab="Residuals",cex.lab=1.15, type="pearson")
mtext("A",side=2,line=2,at=4,col='black',font=2,las=1,size=1.75)
plot(m2,which = 4, caption=NULL, cex.lab=1.15, type="pearson") 
mtext("B",side=2,line=2,at=3.5,col='black',font=2,las=1,size=1.75)
plot(m2, which = 5, caption = NULL, type="deviance", cex.lab=1.15)
mtext("C",side=2,line=2,at=3.5,col='black',font=2,las=1,size=1.75)
plot(m2,which=2, type="pearson",caption=NULL, cex.lab=1.15)
mtext("D",side=2,line=2,at=0.43,col='black',font=2,las=1,size=1.75)
dev.off()
# model summary and fit statistics
summary(m2)
AIC(m2)
exp(logLik(m2))
res.dev<-residuals(m2,type = "deviance")
(sum(res.dev^2))/109
sum(res.dev)/109
#ANODEV
m2.intercept<-betareg(AC~1,data=soilac)
m2.OHac<-betareg(AC~1+Region,data=soilac)
m2.CR<-betareg(AC~1+Region+Layer,data=soilac)
m2.OHacCR<-betareg(AC~1+Region+Layer+Layer*Region,data=soilac)
m2ANODEV<-lrtest(m2.intercept,m2.OHac,m2.CR,m2.OHacCR)
m2ANODEV

#### Model Comparison ####
# ANOVA Table GLM
glm.model.ac<-as.data.frame(m1ANODEV)%>%
  rename(numDf='#Df')
# ANOVA with parameter comparison for GLM
acglm.soilac<-glm.model.ac%>%
  mutate(Deviance=2*LogLik)%>%
  mutate(dDeviance=Deviance-lag(Deviance))%>%
  mutate(LR=exp(dDeviance/2))%>%
  mutate(AIC=dDeviance-2*numDf)%>%
  mutate(LR_Df=LR/Df)%>%
  mutate(dAIC=AIC-lag(AIC))

# ANOVA Table Beta
beta.model.ac<-as.data.frame(m2ANODEV)%>%
  rename(numDf='#Df')
# ANOVA with parameter comparison for beta
acbeta.soilac <-beta.model.ac%>%
  mutate(Deviance=2*LogLik)%>%
  mutate(dDeviance=Deviance-lag(Deviance))%>%
  mutate(LR=exp(dDeviance/2))%>%
  mutate(AIC=dDeviance-2*numDf)%>%
  mutate(LR_Df=LR/Df)%>%
  mutate(dAIC=AIC-lag(AIC))

# save tables
write.csv(acglm.soilac,file ="output/glm_alkyl-carbon_model.csv",row.names=FALSE)
write.csv(acbeta.soilac,file ="output/beta_alkyl-carbon_model.csv",row.names = FALSE)

#### Model Parameters ####
m1$coefficients
summary(m1)

summary(m2)
betacoef<-m2$coefficients
odds<-exp(betacoef$mean)
beta.prob<-odds/(1+odds)

parameter.summary<-data.frame(glm=m1$coefficients,
                              beta.coef=betacoef,
                              beta.prob=beta.prob)

# save table
write.csv(parameter.summary,"output/alkyl-carbon_parameters.csv",row.names = FALSE)

#### Fit Statistics ####
# Calculate r-squared for LR
# LR = (1-R^2)^(-n/2), where n is the number of observations
# Use Cox and Snell method (Cox and Snell 1989)
# GLM LR
PseudoR2(m1, which = "CoxSnell") # r-squared for glm, extract cox and snell
(1-PseudoR2(m1, which = "CoxSnell"))^(-36/2) 
# BetaReg LR
(1-m2$pseudo.r.squared)^(-36/2) # LR for {betareg}

#predicted values for RSME calculation
data.predicted.lm<-c(predict(m1,na.omit(soilac)))
data.predicted.beta<-(c(predict(m2,na.omit(soilac))))

fit.soilac<-data.frame(Fit.statistics=c('AIC','LR','R^2','RMSE','Residual Deviance'),
                       glm=c(AIC(m1), # AIC
                             (1-PseudoR2(m1, which = "CoxSnell"))^(-36/2), # LR
                             PseudoR2(m1, which = "CoxSnell"), #R^2
                             rmse(na.omit(soilac$AC_arc), data.predicted.lm), #RSME
                             sum(m1$residuals^2)/m1$df.residual), # Residual Deviance
                       beta=c(AIC(m2), # AIC
                              (1-m2$pseudo.r.squared)^(-36/2), # LR
                              m2$pseudo.r.squared, # R^2
                              rmse(na.omit(soilac$AC),data.predicted.beta), #RSME
                              sum(m2$residuals^2)/m2$df.residual)) # Residual Deviance

# save table
write.csv(fit.soilac,"output/alkyl-carbon_fit.csv",row.names = FALSE)
