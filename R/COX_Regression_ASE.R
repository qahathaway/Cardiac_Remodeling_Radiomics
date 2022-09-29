
##Load Packages##
library(tidyverse)
library(reshape2)
library(xgboost)
library(randomForest)
library(rfUtilities)
library(caret)
library(survival)
library(survminer)
library(randomForestSRC)

##Load Dataset##
MESA <- data.frame(read.csv(
  file = '/Users/quincy/Documents/Research/HVI/Radiomics - Automated/FINAL/Submission/Nature Cardiovascular Research/C-Statistic/NEW_Univariate.csv'))

##Impute Median Values##
library(mlr)
imputed = impute(MESA, target = character(0), classes = list(numeric = imputeMedian(), integer = imputeMedian()))
final <- as.data.frame(imputed$data)

colnames(final)

##Univariate##
res.cox <- coxph(Surv(MACE_days, MACE) ~ARIC3+LV.Mass.Index+LVEF+Probability.of.Myocardial.Remodeling, data = final)
res.cox
summary(res.cox)

res.cox2 <- coxph(Surv(MACE_months, MACE) ~ IVS_original_gldm_LowGrayLevelEmphasis, data = final)
res.cox2
summary(res.cox2)
        
covariates <- c("Moderate.Severe.Valve", "LV.Mass.Index",
                "LV.EDV", "LA.ESVI", "LVEF", "Average.e.prime", "Average.Ee.prime.ratio",
                "MV.EA.ratio", "WMSI", "Pul.HTN", "Probability.of.Myocardial.Remodeling")

covariates <- c("ARIC3")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(MACE_days, MACE)~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = final)})
# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
as.data.frame(res)


##Multivariate##
library(survcomp)

cARIC <- concordance.index(x=final$ARIC3,
                          method="noether", surv.time=final$MACE_days, surv.event=final$MACE)

cLVH <- concordance.index(x=final$ARIC3+final$LV.Mass.Index,
                           method="noether", surv.time=final$MACE_days, surv.event=final$MACE)

cLVEF <- concordance.index(x=final$ARIC3+final$LVEF,
                          method="noether", surv.time=final$MACE_days, surv.event=final$MACE)

cRAD <- concordance.index(x=final$ARIC3+final$Probability.of.Myocardial.Remodeling,
                           method="noether", surv.time=final$MACE_days, surv.event=final$MACE)

cindex.comp(cARIC, cRAD)



res.ARIC <- coxph(Surv(MACE_days, MACE) ~ARIC, data = final)
res.Valve <- coxph(Surv(MACE_days, MACE) ~ARIC+Moderate.Severe.Valve, data = final)
res.LVEF <- coxph(Surv(MACE_days, MACE) ~ARIC+LVEF, data = final)
res.LVH <- coxph(Surv(MACE_days, MACE) ~LV.Mass.Index, data = final)
res.LVHLVEF <- coxph(Surv(MACE_days, MACE) ~LV.Mass.Index+LVEF, data = final)
res.LVHLVEFValve <- coxph(Surv(MACE_days, MACE) ~LV.Mass.Index+LVEF+Moderate.Severe.Valve, data = final)
res.ALL <- coxph(Surv(MACE_days, MACE) ~LV.Mass.Index+LVEF+Moderate.Severe.Valve+Probability.of.Myocardial.Remodeling, data = final)

summary(res.ARIC)
summary(res.Valve)
summary(res.LVEF)
summary(res.LVH)
summary(res.LVHLVEF)
summary(res.LVHLVEFValve)
summary(res.ALL)
concordance(res.ARIC)
concordance(res.LVH)
concordance(res.LVEF)
concordance(res.RAD)
anova(res.ARIC, res.LVH)
anova(res.ARIC, res.LVEF)
anova(res.ARIC, res.RAD)


##Plot the baseline survival function##
fit <- surv_fit(Surv(Cardiovascular_days, Cardiovascular) ~Prob,
                data = final)

ggsurvplot(fit, data = final, pval = TRUE, break.time.by = 200)

ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for 
  fun="event",
  conf.int.style = "step",  # customize style of confidence intervals
  xlab = "Time in days",
  ylab = "MACE Probability",# customize X axis label.
  fontsize = 2,
  break.time.by = 30,     # break X axis in time intervals by 200.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table = "abs_pct",  # absolute number and percentage at risk.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.y.text = TRUE,
  risk.table.font = 3,
  risk.table.title = "Number at Risk (%)",
  risk.table.pos = "out",
  cumevents.title = "Cumulative Events",
  cumevents = FALSE,
  cumevents.font = 3,
  font.x = 12,
  font.y = 12,
  ncensor.plot = FALSE,      # plot the number of censored subjects at time t
  surv.median.line = "hv",  # add the median survival pointer.
  ylim = c(0, .25),
  legend.labs = 
    c("Normal", "MACE"),    # change legend labels.
  palette = 
    c("black", "red") # custom color palettes.
)

# Fit a Cox proportional hazards model
surv_object <- Surv(time = final$MACE_days, event = final$MACE)
fit.coxph <- coxph(surv_object ~Prob_Bin, 
                   data = final)
ggforest(fit.coxph, data = final)


####RandomForestSRC####
train <- sample(1:nrow(final), round(nrow(final) * 0.80))

final.grow <- rfsrc(Surv(ToD, CardiacMortality) ~ Age+Gender+Race+Hypertension+Mean_SBP+Mean_DBP+Diabetes+Hyperlipidemia+
                      LDL+HDL+Chol+Trig+Smoking+Pack_Years+Metabolic_Syndrome+BMI+
                      Homocysteine+Tumor_Necrosis_Factor+IL6+Plasmin_Antiplasmin+Fibrinogen_Antigen+
                      C_Reactive_Protein+D_Dimer+Factor_VIII+Aortic_Valve_Calcium+Mitral_Valve_Calcium+
                      Ascending_Thoracic_Aortic_Calcium+Descending_Thoracic_Aortic_Calcium+Aortic_Valve_Ring_Calcium+
                      Left_Ventricular_Size+CalciumExam1+AnnualizedProgressCalcium+SQRTProgressCalcium, final[train, ], ntree = 1000,  importance = TRUE)

final.pred <- predict(final.grow, final[-train , ])
print(final.grow)
print(final.pred)

plot(final.grow)
plot.survival(final.grow)
plot.variable(final.grow, xvar.names = c("Age", "CalciumExam1"), surv.type = "surv")

## plot survival curves for first 10 individuals -- direct way
matplot(final.grow$time.interest, 100 * t(final.grow$survival.oob[1:10, ]),
        xlab = "Time", ylab = "Survival", type = "l", lty = 1)
## plot survival curves for first 10 individuals
plot.survival(final.grow, subset = 1:10)



####Compare RF-SRC to Cox Regression####
library(survival)
library(pec)
library(prodlim)
library(riskRegression)

if (library("survival", logical.return = TRUE)
    & library("pec", logical.return = TRUE)
    & library("prodlim", logical.return = TRUE))
{
  ##prediction function required for pec
  predictSurvProb.rfsrc <- function(object, newdata, times, ...){
    ptemp <- predict(object,newdata=newdata,...)$survival
    pos <- sindex(jump.times = object$time.interest, eval.times = times)
    p <- cbind(1,ptemp)[, pos + 1]
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
      stop("Prediction failed")
    p
  }
  ## data, formula specifications
  surv.f <- as.formula(Surv(ToD, CardiacMortality) ~ Age+Gender+Race+Hypertension+Mean_SBP+Mean_DBP+Diabetes+Hyperlipidemia+
                         LDL+HDL+Chol+Trig+Smoking+Pack_Years+Metabolic_Syndrome+BMI)
  
  pec.f <- as.formula(Hist(ToD, CardiacMortality) ~ 1)
  ## run cox/rfsrc models
  cox.obj <- coxph(surv.f, data = final, x = TRUE)
  rfsrc.obj <- rfsrc(surv.f, final, ntree = 1000)
  ## compute bootstrap cross-validation estimate of expected Brier score
  set.seed(100)
  prederror.pbc <- pec(list(cox.obj,rfsrc.obj), data = final, formula = pec.f,
                       splitMethod = "bootcv", B = 2)
  print(prederror.pbc)
  plot(prederror.pbc)
  ## compute out-of-bag C-index for cox regression and compare to rfsrc
  rfsrc.obj <- rfsrc(surv.f, final)
  cat("out-of-bag Cox Analysis ...", "\n")
  cox.err <- sapply(1:100, function(b) {
    if (b%%10 == 0) cat("cox bootstrap:", b, "\n")
    train <- sample(1:nrow(final), nrow(final), replace = TRUE)
    cox.obj <- tryCatch({coxph(surv.f, final[train, ])}, error=function(ex){NULL})
    if (!is.null(cox.obj)) {
      get.cindex(final$ToD[-train], final$CardiacMortality[-train], predict(cox.obj, final[-train, ]))
    } else NA
  })
  cat("\n\tOOB error rates\n\n")
  cat("\tRSF : ", rfsrc.obj$err.rate[rfsrc.obj$ntree], "\n")
  cat("\tCox regression : ", mean(cox.err, na.rm = TRUE), "\n")
}


####Compare 2 or more Survival Curves####
#library(FHtest)
#data(bcos)
#FHtesticp(Surv(left, right, type = "interval2")~treatment, data = bcos)

