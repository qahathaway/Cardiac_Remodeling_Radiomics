
####Regression Predictions####
library(PredictABEL)
library(prodlim)
library(riskRegression)

##Load Dataset##
MESA <- data.frame(read.csv(
  file = 'path/to/file.csv'))

##Impute Median Values##
library(mlr)
imputed = impute(MESA, target = character(0), classes = list(numeric = imputeMedian(), integer = imputeMedian()))
final <- as.data.frame(imputed$data)

# detach mlr package because of 'plotCalibration' conflict
detach("package:mlr", unload = TRUE)

fgrARIC<-FGR(Hist(MACE_days, MACE)~ARIC, data=final, cause=1)
risk_score_ARIC <- predictRisk(fgrARIC,times=2,newdata=final)

fgrLVH<-FGR(Hist(MACE_days, MACE)~ARIC+LV.Mass.Index, data=final, cause=1)
risk_score_LVH <- predictRisk(fgrLVH,times=2,newdata=final)

fgrLVEF<-FGR(Hist(MACE_days, MACE)~ARIC+LV.Mass.Index+LVEF, data=final, cause=1)
risk_score_LVEF <- predictRisk(fgrLVEF,times=2,newdata=final)

fgrRAD<-FGR(Hist(MACE_days, MACE)~ARIC+LV.Mass.Index+LVEF+Ultrasomics, data=final, cause=1)
risk_score_RAD <- predictRisk(fgrRAD,times=2,newdata=final)


df_list <- list(risk_score_LVEF, risk_score_RAD)
df <- data.frame(df_list)

write.csv(df, file = "path/to/file.csv")


###Sample 1###
# specify column number of outcome variable
cOutcome <- 1
# specify column numbers of non-genetic predictors
cNonGenPred <- c(13)
# specify column numbers of non-genetic predictors that are categorical
cNonGenPredCat <- c(0)
# specify column numbers of genetic predictors
cGenPred <- c(0)
# specify column numbers of genetic predictors that are categorical
cGenPredCat <- c(0)
# fit logistic regression model
riskmodel1 <- fitLogRegModel(data=final, cOutcome=cOutcome,
                            cNonGenPreds=cNonGenPred, cNonGenPredsCat=cNonGenPredCat,
                            cGenPreds=cGenPred, cGenPredsCat=cGenPredCat)
summary(riskmodel1)

###Sample 2###
# specify column numbers of non-genetic predictors
cNonGenPred <- c(14)
# fit logistic regression model
riskmodel2 <- fitLogRegModel(data=final, cOutcome=cOutcome,
                            cNonGenPreds=cNonGenPred, cNonGenPredsCat=cNonGenPredCat,
                            cGenPreds=cGenPred, cGenPredsCat=cGenPredCat)
summary(riskmodel2)


# obtain multivariate OR(95% CI) for all predictors of the fitted model
ORmultivariate(riskModel=riskmodel1)
ORmultivariate(riskModel=riskmodel2)

# obtain predicted risks
predRisk1 <- predRisk(riskmodel1)
predRisk2 <- predRisk(riskmodel2)

#compute reclassification measures
reclassification(data=final, cOutcome=cOutcome,
                 predrisk1=predRisk1, predrisk2=predRisk2)


# compute calibration measures and produce calibration plot
plotCalibration(data=final, cOutcome=cOutcome, predRisk=predRisk1, groups = 100)
plotCalibration(data=final, cOutcome=cOutcome, predRisk=predRisk2, groups = 100)
