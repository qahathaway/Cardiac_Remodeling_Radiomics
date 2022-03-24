
##Load Dataset##
dat <- data.frame(read.csv(file = '/path/to/file.csv'))

#Load Packages
library(survival)
library(pec)
library(prodlim)

#COXPH Regression
coxProb <- coxph(Surv(duration,event)~Prob,data=dat,x=TRUE,y=TRUE)

coxALL <- coxph(Surv(duration,event)~.,data=dat,x=TRUE,y=TRUE)

#Random Survival Forest
rsfProb <- rfsrc(Surv(duration,event)~Prob,data=dat,ntree=1000,forest=TRUE)

rsfALL <- rfsrc(Surv(duration,event)~.,data=dat,ntree=1000,forest=TRUE)

#Compute the apparent estimate of the C-index at different time points - External Validation (POCUS)
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
write.csv(ApparrentCindex2$AppCindex, file = "/path/to/file.csv")

#Compute the apparent estimate of the C-index at different time points - External Validation (High-End)
ApparrentCindex <- pec::cindex(list("COXPH Prob"=coxProb,
                                     "COXPH ALL"=coxALL),
              formula=Surv(duration,event)~.,data=dat,
              eval.times=seq(0,266,1), pred.times=seq(0,266,1))

ApparrentCindex2 <- pec::cindex(list("COXPH Prob"=rsfProb,
                                    "COXPH ALL"=rsfALL),
                               formula=Surv(duration,event)~.,data=dat,
                               eval.times=seq(0,266,1), pred.times=seq(0,266,1))

plot(ApparrentCindex, legend = FALSE, xlim=c(0,300))
plot(ApparrentCindex2, legend = FALSE, xlim=c(0,300))
write.csv(ApparrentCindex2$AppCindex, file = "/path/to/file.csv")

