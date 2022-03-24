
##Load Packages##
library(reshape2)
library(survival)
library(survminer)

##Load Dataset##
final <- data.frame(read.csv(
  file = '/path/to/file.csv'))

colnames(final)

##Univariate##
res.cox <- coxph(Surv(MACE_months, MACE) ~., data = final)
res.cox
summary(res.cox)
        
covariates <- c("IVS_original_gldm_LargeDependenceLowGrayLevelEmphasis",
                "IVS_original_gldm_LowGrayLevelEmphasis","Myocardial.Remodeling","Age_x",
                "Gender","BMI","Systolic.BP","Diastolic.BP","NYHA.class","CHF","CAD","A.Fib","CVA","DM","HLD","COPD",
                "Pul.HTN","Moderate.Severe.Valve","LV.Mass.Index","IVSd","PWTd","RWT",
                "LVIDd","LV.EDV","LA.ESVI","LVEF","Average.e.prime","Average.Ee.prime.ratio",
                "MV.EA.ratio","WMSI","RWMA")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(MACE_months, MACE)~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = final)})
                        
##Extract data##
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
res.cox <- coxph(Surv(MACE_months, MACE) ~.,
                 data = final)

res.cox1 <- coxph(Surv(MACE_months, MACE) ~Myocardial.Remodeling+
                   Age_x+CAD+
                   LV.Mass.Index+Average.e.prime,
                 data = final)

res.cox3 <- coxph(Surv(MACE_months, MACE) ~IVS_original_gldm_LowGrayLevelEmphasis+
                    Age_x+CAD+
                    LV.Mass.Index+Average.e.prime,
                 data = final)
                        
res.cox2 <- coxph(Surv(MACE_months, MACE) ~IVS_original_gldm_LargeDependenceLowGrayLevelEmphasis+
                    Age_x+CAD+
                    LV.Mass.Index+Average.e.prime,
                 data = final)

summary(res.cox)
summary(res.cox1)
summary(res.cox2)
summary(res.cox3)
