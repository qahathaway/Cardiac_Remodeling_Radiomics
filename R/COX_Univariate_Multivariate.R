
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
  file = 'path/to/file.csv'))

##Impute Median Values##
library(mlr)
imputed = impute(MESA, target = character(0), classes = list(numeric = imputeMedian(), integer = imputeMedian()))
final <- as.data.frame(imputed$data)

colnames(final)

##Univariate##     
covariates <- c("Moderate.Severe.Valve", "LV.Mass.Index",
                "LV.EDV", "LA.ESVI", "LVEF", "Average.e.prime", "Average.Ee.prime.ratio",
                "MV.EA.ratio", "WMSI", "Pul.HTN", "Ultrasomics")

covariates <- c("ARIC")

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

cARIC <- concordance.index(x=final$ARIC,
                          method="noether", surv.time=final$MACE_days, surv.event=final$MACE)

cLVH <- concordance.index(x=final$ARIC+final$LV.Mass.Index,
                           method="noether", surv.time=final$MACE_days, surv.event=final$MACE)

cLVEF <- concordance.index(x=final$ARIC+final$LVEF,
                          method="noether", surv.time=final$MACE_days, surv.event=final$MACE)

cRAD <- concordance.index(x=final$ARIC+final$Ultrasomics,
                           method="noether", surv.time=final$MACE_days, surv.event=final$MACE)

cindex.comp(cARIC, cRAD)


res.ARIC <- coxph(Surv(MACE_days, MACE) ~ARIC, data = final)
res.Valve <- coxph(Surv(MACE_days, MACE) ~ARIC+Moderate.Severe.Valve, data = final)
res.LVEF <- coxph(Surv(MACE_days, MACE) ~ARIC+LVEF, data = final)
res.LVH <- coxph(Surv(MACE_days, MACE) ~LV.Mass.Index, data = final)
res.LVHLVEF <- coxph(Surv(MACE_days, MACE) ~LV.Mass.Index+LVEF, data = final)
res.LVHLVEFValve <- coxph(Surv(MACE_days, MACE) ~LV.Mass.Index+LVEF+Moderate.Severe.Valve, data = final)
res.ALL <- coxph(Surv(MACE_days, MACE) ~LV.Mass.Index+LVEF+Moderate.Severe.Valve+Ultrasomics, data = final)

summary(res.ARIC)
concordance(res.ARIC)
anova(res.ARIC, res.ALL)


##Plot the baseline survival function##
fit <- surv_fit(Surv(MACE_days, MACE) ~Prob,
                data = final)

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
