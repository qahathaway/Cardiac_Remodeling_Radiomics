
####Regression Predictions####
library(PredictABEL)
library(prodlim)
library(riskRegression)

##Load Dataset##
MESA <- data.frame(read.csv(
  file = '/Users/quincy/Documents/Research/HVI/Radiomics - Automated/FINAL/Submission/Nature Cardiovascular Research/C-Statistic/NEW_PredScore.csv'))

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

fgrRAD<-FGR(Hist(MACE_days, MACE)~ARIC+LV.Mass.Index+LVEF+Probability.of.Myocardial.Remodeling, data=final, cause=1)
risk_score_RAD <- predictRisk(fgrRAD,times=2,newdata=final)


df_list <- list(risk_score_LVEF, risk_score_RAD)
df <- data.frame(df_list)

write.csv(df, file = "/Users/quincy/Documents/Research/HVI/Radiomics - Automated/FINAL/Submission/Nature Cardiovascular Research/C-Statistic/NEW_PredScore2.csv")


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

# specify cutoff values for risk categories
cutoff <- c(0.0,0.15,1.0)
#cutoff <- c(0, 0.27, 0.36, 0.78, 1)
#compute reclassification measures
reclassification(data=final, cOutcome=cOutcome,
                 predrisk1=predRisk1, predrisk2=predRisk2, cutoff=cutoff)




# compute calibration measures and produce calibration plot
plotCalibration(data=final, cOutcome=cOutcome, predRisk=predRisk1, groups = 100)
plotCalibration(data=final, cOutcome=cOutcome, predRisk=predRisk2, groups = 100)


# specify labels for the groups without and with the outcome of interest
labels <- c("No Mortality", "Mortality")
# produce discrimination box plot
plotDiscriminationBox(data=final, cOutcome=cOutcome, predrisk=predRisk1,
                      labels=labels)
plotDiscriminationBox(data=final, cOutcome=cOutcome, predrisk=predRisk2,
                      labels=labels)

# specify range of y-axis
rangeyaxis <- c(0,1)
# specify labels of the predictiveness curves
labels <- c("B", "Framingham")
# produce predictiveness curves
plotPredictivenessCurve(predrisk=cbind(predRisk1,predRisk2),
                        rangeyaxis=rangeyaxis, labels=labels)


# specify the size of each interval
interval <- .05
# specify label of x-axis
xlabel <- "Predicted risk"
# specify label of y-axis
ylabel <- "Percentage"
# specify range of x-axis
xrange <- c(0,1)
# specify range of y-axis
yrange <- c(0,40)
# specify title for the plot
maintitle <- "Distribution of predicted risks"
# specify labels
labels <- c("No Mortality", "Mortality")
# produce risk distribution plot
plotRiskDistribution(data=final, cOutcome=cOutcome,
                     risks=predRisk1, interval=interval, plottitle=maintitle, rangexaxis=xrange,
                     rangeyaxis=yrange, xlabel=xlabel, ylabel=ylabel, labels=labels)

# specify range of x-axis
rangexaxis <- c(0,12)
# specify range of y-axis
rangeyaxis <- c(0,1)
# specify label of x-axis
xlabel <- "Risk score"
# specify label of y-axis
ylabel <- "Predicted risk"
# specify title for the plot
plottitle <- "Risk score versus predicted risk"
# function to compute unweighted genetic risk scores
riskScore1 <- riskScore(weights=riskmodel1, data=final,
                       cGenPreds=cGenPred, Type="unweighted")

# produce risk score-predicted risk plot
plotRiskscorePredrisk(data=final, riskScore=riskScore1, predRisk=predRisk1,
                      plottitle=plottitle, xlabel=xlabel, ylabel=ylabel, rangexaxis=rangexaxis,
                      rangeyaxis=rangeyaxis)


# specify cutoff values for risk categories
cutoff <- c(0,.10,.20,.30,.40,.50,.60,.70,.80,.90,1)
# compute reclassification measures
reclassification(data=final, cOutcome=cOutcome,
                 predrisk1=predRisk1, predrisk2=predRisk2, cutoff=cutoff)


#--- sample data (pbc in survival package) ---
D=subset(pbc, select=c("time","status","age","albumin","edema","protime","bili"))
D$status=as.numeric(D$status==2)
D=D[!is.na(apply(D,1,mean)),] ; dim(D)
mydata=D[1:100,]





####IDI Survival####
library(survIDINRI)

final <- data.frame(read.csv(
  file = '/Users/quincy/Documents/Research/HVI/Radiomics - Automated/FINAL/Submission/Nature Cardiovascular Research/C-Statistic/NEW_PredScore.csv'))

library(mlr)
imputed = impute(final, target = character(0), classes = list(numeric = imputeMedian(), integer = imputeMedian()))
final1 <- as.data.frame(imputed$data)
final1$event=as.integer(final1$event)
final1$duration=as.integer(final1$duration)
detach("package:mlr", unload = TRUE)

covs0<-as.matrix(final1[5])
covs1<-as.matrix(final1[6])
#--- inference ---
t0=1140

x<-IDI.INF(final1[,2:1], covs0, covs1, t0, npert=10)
#--- results ---
IDI.INF.OUT(x)
#--- Graphical presentaion of the estimates###
IDI.INF.GRAPH(x)
