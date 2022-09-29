library(predtools)
library(magrittr)
library(dplyr)
library(ggplot2)
library(riskRegression)
library(prodlim)
library(glmtoolbox)
library(aplore3)
library(rms)
library(rfUtilities)

ASE <- data.frame(read.csv(
  file = '/Users/quincy/Documents/Research/HVI/Radiomics - Automated/FINAL/Submission/JACC/Editor_Reviewer_Comments/Resubmission/Appeal/Revised/ASE-REWARD2.csv'))

CHOICE <- data.frame(read.csv(
  file = '/Users/quincy/Documents/Research/HVI/Radiomics - Automated/FINAL/Submission/JACC/Editor_Reviewer_Comments/Resubmission/Appeal/Revised/CHOICE.csv'))

WVU <- data.frame(read.csv(
  file = '/Users/quincy/Documents/Research/HVI/Radiomics - Automated/FINAL/Submission/JACC/Editor_Reviewer_Comments/Resubmission/Appeal/Revised/WVU.csv'))

reg <- glm(Outcome~.,data=ASE,family=binomial(link="logit"))
summary(reg)

ASE$pred<- predict.glm(reg, type = 'response')
CHOICE$pred <- predict.glm(reg, newdata = CHOICE, type = 'response')
WVU$pred <- predict.glm(reg, newdata = WVU, type = 'response')

calibration_plot(data = ASE, obs = "Outcome", pred = "pred", title = "Calibration plot for development data")

calibration_plot(data = CHOICE, obs = "Outcome", pred = "pred", y_lim = c(0, 0.6),
                 title = "Calibration plot for validation data", group = "Outcome")



# binary
fb1=glm(Outcome~pred,data=ASE,family="binomial")
fb2=glm(Outcome~calibrated,data=ASE,family="binomial")
ASE$pred <- predict.glm(fb1, type = 'response')
y.hat <- predict(fb1, type = 'response')
y.hat <- format(round(y.hat, 3), nsmall = 3)
y.hat <- as.numeric(y.hat)
y <- as.numeric(ASE$Outcome)

calibrated <- probability.calibration(y, y.hat, regularization = TRUE)
ASE$calibrated <- calibrated

write.csv(ASE, file = "/Users/quincy/Documents/Research/HVI/Radiomics - Automated/FINAL/Submission/JACC/Editor_Reviewer_Comments/Resubmission/Appeal/Revised/ASEPred.csv")

xb=Score(list(model1=fb1, model2=fb2),Outcome~1,data=ASE, plots=c("calibration","ROC"))
xb
plotROC(xb)
plotCalibration(xb)
plotCalibration(xb,bars=TRUE,model="model1")
plotCalibration(xb,models=1,bars=TRUE,names.cex=1.3)
hltest(fb2)

