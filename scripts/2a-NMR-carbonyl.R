# Carbonyl Carbon Analyses
# Comparing beta regression to general linear model with arcsin transformation. 
# Using carbonyl carbon in soil - modelling region and soil layer
# comparing models with likelihood ratios, residual deviance/residual, and AICc.

#### Packages ####
# load required packages 
p <- c("betareg", "lmtest", "DescTools", "dplyr","Metrics")
lapply(p, library, character.only = T)

#### Data #### 
# load data
soil <- read.csv("input/Beta_NMRdata.csv")
# select relevant columns
soilcc <- select(soil, c(Region, Layer, CC))
# add carbonyl carbon data with arcsin square root transformation
soilcc <- mutate(soilcc, CC_arc = asin(sqrt(soilcc$CC)))

#### General Linear Model with Arcsin Transformation ####
m1.soilcc <- glm(CC_arc ~ Region+Layer+Region*Layer,data=na.omit(soilcc))
m1.soilccFrameNormal<-cbind(soilcc,residuals(m1.soilcc),fitted(m1.soilcc))
# diagnostic plots for manuscript
png("figures/glm_carbonyl-carbon_diagnostics.png", width = 160, height = 160, units = "mm",res = 600)
par(mfrow=c(2,2))
plot(x=fitted(m1.soilcc),y=resid(m1.soilcc),main=NULL,
     xlab="Fitted Values",ylab="Residuals",cex.lab=1.15)
mtext("A",side=2,line=2,at=0.02,col='black',font=2,las=1,size=1.75)
hist(resid(m1.soilcc),main=NULL,xlab="Residuals",
     cex.lab=1.15)
mtext("B",side=2,line=2,at=12,col='black',font=2,las=1,size=1.75)
qqnorm(resid(m1.soilcc),main=NULL,cex.lab=1.15)
mtext("C",side=2,line=2,at = 0.025,col='black',font=2,las=1,size=1.75)
qqline(resid(m1.soilcc),col='red')
plot(m1.soilcc,which=4,caption = NULL,main = NULL,cex.lab=1.15)
mtext("D",side=2,line=2,at=0.2,col='black',font=2,las=1,size=1.75)
dev.off()
# model summary and fit statistics
summary(m1.soilcc)
exp(logLik(m1.soilcc))
anova(m1.soilcc)
#ANODEV
m1.intercept<-glm(CC_arc~1,data=na.omit(soilcc))
m1.OHCC<-glm(CC_arc~1+Region,data=na.omit(soilcc))
m1.CR<-glm(CC_arc~1+Region+Layer,data=na.omit(soilcc))
m1.OHCCCR<-glm(CC_arc~1+Region+Layer+Region*Layer,data=na.omit(soilcc))
m1ANODEV<-lrtest(m1.intercept,m1.OHCC,m1.CR,m1.OHCCCR)
m1ANODEV

#### Beta Regression ####
m2 <- betareg(CC~Region+Layer+Region*Layer,
              link="logit",link.phi = "identity",
              data=soilcc)
m2FrameBeta<-cbind(soilcc,residuals(m2),fitted(m2))
# diagnostic plots for manuscript
png("figures/beta_carboynl-carbon_diagnostics.png",  width = 160, height = 160, units = "mm",res = 600)
par(mfrow=c(2,2))
plot(x=fitted(m2),y=resid(m2),main=NULL,
     xlab="Fitted Values",ylab="Residuals",cex.lab=1.15)
mtext("A",side=2,line=2,at=3,col='black',font=2,las=1,size=1.75)
plot(m2,which = 4, caption=NULL, cex.lab=1.15, type="deviance") 
mtext("B",side=2,line=2,at=2.6,col='black',font=2,las=1,size=1.75)
plot(m2, which = 5, caption = NULL, type="deviance", cex.lab=1.15)
mtext("C",side=2,line=2,at=3.5,col='black',font=2,las=1,size=1.75)
plot(m2,which=2, type="deviance",caption=NULL, cex.lab=1.15)
mtext("D",side=2,line=2,at=0.26,col='black',font=2,las=1,size=1.75)
dev.off()
# model summary and fit statistics
summary(m2)
AIC(m2)
exp(logLik(m2))
res.dev<-residuals(m2,type = "deviance")
(sum(res.dev^2))/109
sum(res.dev)/109
#ANODEV
m2.intercept<-betareg(CC~1,data=soilcc)
m2.OHCC<-betareg(CC~1+Region,data=soilcc)
m2.CR<-betareg(CC~1+Region+Layer,data=soilcc)
m2.OHCCCR<-betareg(CC~1+Region+Layer+Layer*Region,data=soilcc)
m2ANODEV<-lrtest(m2.intercept,m2.OHCC,m2.CR,m2.OHCCCR)
m2ANODEV

#### Model Comparison ####
# ANOVA Table GLM
glm.model.CC<-as.data.frame(m1ANODEV)%>%
  rename(numDf='#Df')
# ANOVA with parameter comparison for GLM
CCglm.soilcc<-glm.model.CC%>%
  mutate(Deviance=2*LogLik)%>%
  mutate(dDeviance=Deviance-lag(Deviance))%>%
  mutate(LR=exp(dDeviance/2))%>%
  mutate(AIC=dDeviance-2*numDf)%>%
  mutate(LR_Df=LR/Df)%>%
  mutate(dAIC=AIC-lag(AIC))

# ANOVA Table Beta
beta.model.CC<-as.data.frame(m2ANODEV)%>%
  rename(numDf='#Df')
# ANOVA with parameter comparison for beta
CCbeta.soilcc <-beta.model.CC%>%
  mutate(Deviance=2*LogLik)%>%
  mutate(dDeviance=Deviance-lag(Deviance))%>%
  mutate(LR=exp(dDeviance/2))%>%
  mutate(AIC=dDeviance-2*numDf)%>%
  mutate(LR_Df=LR/Df)%>%
  mutate(dAIC=AIC-lag(AIC))

# save tables
write.csv(CCglm.soilcc,file ="output/glm_carbonyl-carbon_model.csv",row.names=FALSE)
write.csv(CCbeta.soilcc,file ="output/beta_carbonyl-carbon_model.csv",row.names = FALSE)

#### Model Parameters ####
m1.soilcc$coefficients
m1.squared<-(sin(m1.soilcc$coefficients))^2
summary(m1.soilcc)

summary(m2)
betacoef<-m2$coefficients
odds<-exp(betacoef$mean)
beta.prob<-odds/(1+odds)

parameter.summary<-data.frame(glm.coef=m1.soilcc$coefficients,
                              glm2=m1.squared,
                              beta.coef=betacoef,
                              beta.prob=beta.prob)

# save table
write.csv(parameter.summary,"output/carbonyl-carbon_parameters.csv",row.names = FALSE)

#### Fit Statistics ####
# Calculate r-squared for LR
# LR = (1-R^2)^(-n/2), where n is the number of observations
# Use Cox and Snell method (Cox and Snell 1989)
# GLM LR
PseudoR2(m1.soilcc, which = "CoxSnell") # r-squared for glm, extract cox and snell
(1-PseudoR2(m1.soilcc, which = "CoxSnell"))^(-36/2) 
# BetaReg LR
(1-m2$pseudo.r.squared)^(-36/2) # LR for {betareg}

#predicted values for RSME calculation
data.predicted.lm<-c(predict(m1.soilcc,na.omit(soilcc)))
data.predicted.beta<-(c(predict(m2,na.omit(soilcc))))

fit.soilcc<-data.frame(Fit.statistics=c('AIC','LR','R^2','RMSE','Residual Deviance'),
                     glm=c(AIC(m1.soilcc), # AIC
                           (1-PseudoR2(m1.soilcc, which = "CoxSnell"))^(-36/2), # LR
                           PseudoR2(m1.soilcc, which = "CoxSnell"), #R^2
                           rmse(na.omit(soilcc$CC_arc), data.predicted.lm), #RSME
                           sum(m1.soilcc$residuals^2)/m1.soilcc$df.residual), # Residual Deviance
                     beta=c(AIC(m2), # AIC
                            (1-m2$pseudo.r.squared)^(-36/2), # LR
                            m2$pseudo.r.squared, # R^2
                            rmse(na.omit(soilcc$CC),data.predicted.beta), #RSME
                            sum(m2$residuals^2)/m2$df.residual)) # Residual Deviance

# save table
write.csv(fit.soilcc,"output/carbonyl-carbon_fit.csv",row.names = FALSE)

