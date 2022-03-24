
##Load Dataset##
MESA <- data.frame(read.csv(file = '/Users/quincy/Documents/Research/HVI/Radiomics - Automated/FINAL/Submission/Nature Cardiovascular Research/C-Statistic/CHOICE_2.csv'))

##Impute Median Values##
detach("package:randomForestSRC", unload = TRUE)
library(mlr)
imputed = impute(MESA, target = character(0), classes = list(numeric = imputeMedian(), integer = imputeMedian()))
dat <- as.data.frame(imputed$data)

#Load Packages
library(survival)
library(randomForestSRC)
library(pec)
library(prodlim)

#COXPH Regression
coxProb <- coxph(Surv(duration,event)~Prob,data=dat,x=TRUE,y=TRUE)

coxALL <- coxph(Surv(duration,event)~.,data=dat,x=TRUE,y=TRUE)

#Random Survival Forest
rsfProb <- rfsrc(Surv(duration,event)~Prob,data=dat,ntree=1000,forest=TRUE)

rsfALL <- rfsrc(Surv(duration,event)~.,data=dat,ntree=1000,forest=TRUE)

# compute the apparent estimate of the C-index at different time points
ApparrentCindex <- pec::cindex(list("COXPH Prob"=coxProb,
                                     "COXPH ALL"=coxALL),
              formula=Surv(duration,event)~.,data=dat,
              eval.times=seq(0,189,1), pred.times=seq(0,189,1))

ApparrentCindex2 <- pec::cindex(list("COXPH Prob"=rsfProb,
                                    "COXPH ALL"=rsfALL),
                               formula=Surv(duration,event)~.,data=dat,
                               eval.times=seq(0,189,1), pred.times=seq(0,189,1))

plot(ApparrentCindex, legend = FALSE, xlim=c(0,200))
plot(ApparrentCindex2, legend = FALSE, xlim=c(0,200))
write.csv(ApparrentCindex2$AppCindex, file = "/Users/quincy/Documents/Research/HVI/MESA/FINAL/Testerss.csv")


##cumulative prediction error##

Models <-  list('CoxProb'=coxph(Surv(duration,event)~Prob,data=dat,x=TRUE,y=TRUE),
               'CoxALL'=coxph(Surv(duration,event)~.,data=dat,x=TRUE,y=TRUE))

# compute the apparent prediction error
PredError <- pec(object=Models,
                 formula=Surv(duration,event)~.,
                 data=dat,
                 exact=TRUE,
                 cens.model="marginal",
                 splitMethod="none",
                 B=0,
                 verbose=TRUE)
print(PredError,times=seq(189,189,1))
summary(PredError)
plot(PredError,xlim=c(0,200), legend = FALSE)
