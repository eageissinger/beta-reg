# Comparing beta regression to general linear model with arcsin transformation. 
# Using %N for Abies balsamea (balsam fir) - modelling site and year.
# comparing models with likelihood ratios, residual deviance/residual, and AICc.

#### Packages ####
# load required packages 
p <- c("betareg", "lmtest", "DescTools", "dplyr", "Metrics")
lapply(p, library, character.only = T)

#### Data #### 
# load data
stoich <- read.csv("input/elemental-composition.csv")
# select relevant columns and subset Balsam fir (ABBA)
abba <- stoich %>%
  select(c(Year, Site, Species, N_dec)) %>%
  subset(Species == "ABBA")
# add nitrogen data with arcsin square root transformation
abba <- mutate(abba, N_arc = asin(sqrt(abba$N_dec)))

#### General Linear Model with Arcsin Transformation ####
glmN <- glm(N_arc ~ Year + Site, data = abba)
glmNFrameNormal<-cbind(abba,residuals(glmN),fitted(glmN))
# diagnostic plots for manuscript 
png("figures/glm_elemental-composition_diagnostics.png", width = 160, height = 160, units = "mm",res = 600)
par(mfrow=c(2,2))
plot(x=fitted(glmN),y=resid(glmN),main=NULL,
     xlab="Fitted Values",ylab="Residuals",cex.lab=1.15)
mtext("A",side=2,line=2,at=0.035,col='black',font=2,las=1,size=1.75)
hist(resid(glmN),main=NULL,xlab="Residuals",
     cex.lab=1.15)
mtext("B",side=2,line=2,at=33,col='black',font=2,las=1,size=1.75)
qqnorm(resid(glmN),main=NULL,cex.lab=1.15)
mtext("C",side=2,line=2,at = 0.03,col='black',font=2,las=1,size=1.75)
qqline(resid(glmN),col='red')
plot(glmN,which=4,caption = NULL,main = NULL,cex.lab=1.15)
mtext("D",side=2,line=2,at=0.7,col='black',font=2,las=1,size=1.75)
dev.off()
# model summary and fit statistics 
summary(glmN)
exp(logLik(glmN))
anova(glmN)
# build ANODEV table sequentially, same method as used for beta regression 
glmNintercept<-glm(N_arc~1, data=abba)
glmNyear<-glm(N_arc ~ 1 + Year, data=abba)
glmNyearsite<-glm(N_arc ~ 1 + Year + Site, data=abba)
glmNANODEV <- lrtest(glmNintercept, glmNyear, glmNyearsite)
glmNANODEV

#### Beta Regression ####
brN<-betareg(N_dec ~ Year + Site, link="logit", data=abba)
brNFrameBeta<-cbind(abba,residuals(brN),fitted(brN))
# diagnostic plots for manuscript
png("figures/beta_elemental-composition_diagnostics.png",  width = 160, height = 160, units = "mm",res = 600)
par(mfrow=c(2,2))
plot(x=fitted(brN),y=resid(brN),main=NULL,
     xlab="Fitted Values",ylab="Residuals",cex.lab=1.15, type="pearson")
mtext("A",side=2,line=2,at=4,col='black',font=2,las=1,size=1.75)
plot(brN,which = 4, caption=NULL, cex.lab=1.15, type="pearson") 
mtext("B",side=2,line=2,at=5,col='black',font=2,las=1,size=1.75)
plot(brN, which = 5, caption = NULL, type="deviance", cex.lab=1.15)
mtext("C",side=2,line=2,at=4.5,col='black',font=2,las=1,size=1.75)
plot(brN,which=2, type="pearson",caption=NULL, cex.lab=1.15)
mtext("D",side=2,line=2,at=0.9,col='black',font=2,las=1,size=1.75)
dev.off()

# model summary and fit statistics
summary(brN)
AIC(brN)
exp(logLik(brN))
res.devN<-residuals(brN,type = "deviance")
(sum(res.devN^2))/109
sum(res.devN)/109

# ANODEV
brNintercept<-betareg(N_dec~1, link="logit",data=abba)
brNyear<-betareg(N_dec~1+Year, link="logit",data=abba)
brNyearsite<-betareg(N_dec~1+Year+Site, link="logit",data=abba)
brNANODEV <- lrtest(brNintercept,brNyear,brNyearsite)
brNANODEV #gives the same result as GLMMTMB


#### Model Comparison ####
# ANOVA Table GLM
glm.model.abbaN<-as.data.frame(glmNANODEV)%>%
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
write.csv(glm.abbaN,"output/glm_elemental-composition_model.csv",row.names=FALSE)
write.csv(beta.abbaN,"output/beta_elemental-composition_model.csv",row.names = FALSE)

#### Model Parameters ####
# calculate model parameters for glm and beta
glmN$coefficients
summary(glmN)

summary(brN)
betacoef<-brN$coefficients
odds<-exp(betacoef$mean)
beta.prob<-odds/(1+odds)

# combine model parameters into single table for comparison
parameter.summary<-data.frame(glm=glmN$coefficients,
                              beta.coef=betacoef,
                              beta.prob=beta.prob)
# save table
write.csv(parameter.summary,"output/elemental-composition_parameters.csv",row.names = FALSE)

#### Fit Statistics ####
# Calculate r-squared for LR
# LR = (1-R^2)^(-n/2), where n is the number of observations
# Use Cox and Snell method (Cox and Snell 1989)
# GLM LR
PseudoR2(glmN, which = "CoxSnell") # r-squared for glm, extract cox and snell
(1-PseudoR2(glmN, which = "CoxSnell"))^(-114/2) 
# BetaReg LR
(1-brN$pseudo.r.squared)^(-114/2) # LR for {betareg}

# predicted values for RSME calculation
data.predicted.lm<- c(predict(glmN,abba))
data.predicted.beta<- c(predict(brN,abba))

fit.abbaN<-data.frame(Fit.statistics=c('AIC','LR','R^2',"RMSE",'Residual Deviance'),
                     glm=c(AIC(glmN), # AIC
                           (1-PseudoR2(glmN, which = "CoxSnell"))^(-114/2), #LR
                           PseudoR2(glmN, which = "CoxSnell"), #R^2
                           rmse(abba$N_arc,data.predicted.lm), #RMSE
                           sum(glmN$residuals^2)/glmN$df.residual), #Residual Deviance
                     beta=c(AIC(brN), # AIC
                            (1-brN$pseudo.r.squared)^(-114/2), #LR
                            brN$pseudo.r.squared, # R^2
                            rmse(abba$N_dec,data.predicted.beta), # RMSE
                            sum(brN$residuals^2)/brN$df.residual)) # Residual Deviance
# save table
write.csv(fit.abbaN,"output/elemental-composition_fit.csv",row.names = FALSE)
