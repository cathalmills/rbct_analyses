#R Script for Within Culling Areas

#Load in Libraries
library(glmmTMB)
library(MuMIn)
library(rstanarm)
library(bridgesampling)
library(DHARMa)
require(AER)
require(vcd)
require(data.table)
require(performance)
require(cowplot)
require(rstan)
require(boot)
require(olsrr)
require(bc)
require(lme4)
require(caret)
require(loo)
require(scales)
require(bayesplot)
require(extraDistr)
require(readxl)


#Set directories to your local directories which will contain the input data
  # and output data

# data.dir2 <- ""
# out.dir2 <- ""

#Setting colours for plots
plot_val_cols <- c("#2196F3", "red") #Setting colours for plots
nb_quant_cols <- c("white", "forestgreen", "white")
pois_quant_cols <- c("white", "gold", "white")

#Incidence data from initial cull until September 2005;
  # Data as per Donnelly et Al, Nature (2006): http://dx.doi.org/10.1038/nature04454
rbctconf <- data.table(read_xlsx(file.path(data.dir2,"confirmed_vetnet.xlsx")))
#Incidence data from first follow-up cull until September 2005;
  # Data as per Donnelly et Al, Nature (2006): http://dx.doi.org/10.1038/nature04454
rbctconf_follow <- data.table(read_xlsx(file.path(data.dir2,"confirmed_vetnet_FOLLOW_UP.xlsx")))

#Incidence data from post-trial period (1 year after last proactive cull until 2013); 
rbctconf_after <- data.table(read_xlsx(file.path(data.dir2,"confirmed_rbctconf_after_trial.xlsx")))


#Recode triplets and treatment (proactive culling vs survey-only)
  # to correspond to Torgerson et Al and thus, enable model review
rbctconf$Triplet<-relevel(factor(rbctconf$Triplet), ref="J")
rbctconf$Treatment<-relevel(factor(rbctconf$Treatment), ref="Survey-only")
rbctconf$Treatment

rbctconf_follow$Triplet<-relevel(factor(rbctconf_follow$Triplet), ref="J")
rbctconf_follow$Treatment<-relevel(factor(rbctconf_follow$Treatment), ref="Survey-only")


rbctconf_after$Triplet<-relevel(factor(rbctconf_follow$Triplet), ref="J")
rbctconf_after$Treatment<-relevel(factor(rbctconf_follow$Treatment), ref="Survey-only")


#Model 1 (as per Nature 2006 paper)----
#Fit to initial period
rs1<-glm(Incidence~Treatment+Triplet+log(Hist3yr)+log(Baseline), family = poisson, data=rbctconf)
summary(rs1)
dispersiontest(rs1) #Fail to reject null of equidispersion
performance::check_model(rs1)
?dispersiontest
#Follow-Up Period
rs1_follow<-glm(Incidence~Treatment+Triplet+log(Hist3yr)+log(Baseline), family = poisson, data=rbctconf_follow)
dispersiontest(rs1_follow) #Fail to reject null of equidispersion
performance::check_model(rs1_follow)

rs1_after<-glm(Incidence~Treatment+Triplet+log(Hist3yr)+log(Baseline), family = poisson, data=rbctconf_after)
summary(rs1_after)
dispersiontest(rs1_after) #Fail to reject null of equidispersion
performance::check_model(rs1_after)


#Extracting SE's for estimates of effect size
rs1_se <- sqrt(diag(vcov(rs1)))
rs1_follow_se <- sqrt(diag(vcov(rs1_follow)))
rs1_after_se <- sqrt(diag(vcov(rs1_after)))*sqrt(rs1_after$deviance/rs1_after$df.residual)
rs1_after_se
#Overdispersion-adjusted CI's (for consistency with previous analyses)
exp(rs1$coefficients[2])-1
exp(rs1$coefficients[2]-qnorm(0.975) * rs1_se[2])-1
exp(rs1$coefficients[2]+qnorm(0.975) * rs1_se[2])-1

exp(rs1_follow$coefficients[2])-1
exp(rs1_follow$coefficients[2]-qnorm(0.975) * rs1_follow_se[2])-1
exp(rs1_follow$coefficients[2]+qnorm(0.975) * rs1_follow_se[2])-1

exp(rs1_after$coefficients[2])-1
exp(rs1_after$coefficients[2]-qnorm(0.975) * rs1_after_se[2])-1
exp(rs1_after$coefficients[2]+qnorm(0.975) * rs1_after_se[2])-1



# Setup for leave-one-out cross-validation (LOOCV), 
  # to be employed for each model and each period of analysis (3)
train.control <- trainControl(method = "LOOCV")

#Model names 
mod_vector <- c("rs1", "rs2", "rs3", "rs3a", "rs4", "rs4a", "rs4c",
                "rs5", "rs5a")
rs4
#List of datasets for each period of analysis
dataset_vector <- list(rbctconf, rbctconf_follow, rbctconf_after)

#Data.table will store model names and LOOCV accuracy metrics
cv_mods <- data.table(mod = mod_vector, 
                                     RMSE = rep(0, length(mod_vector)),
                                     MAE = rep(0, length(mod_vector)), 
                                     Rsquared = rep(0, length(mod_vector)))

#Using caret::train for LOOCV procedure
rs1_cv <- train(Incidence~Treatment+Triplet+log(Hist3yr)+log(Baseline), family = poisson, data=rbctconf, method = "glm",
                trControl = train.control)
cv_mods[which(mod == "rs1"), RMSE:= rs1_cv$results$RMSE]
cv_mods[which(mod == "rs1"), MAE:= rs1_cv$results$MAE]

#As above, with definition of all models for period from follow-up cull until 2005.
mod_vector_follow <- c("rs1_follow", "rs2_follow", "rs3_follow", "rs3a_follow", "rs4_follow", "rs4a_follow",
                       "rs4c_follow", "rs5_follow", "rs5a_follow")

#Data.table will store model names and LOOCV accuracy metrics
cv_mods_follow <- data.table(mod = mod_vector_follow, 
                                            RMSE = rep(0, length(mod_vector_follow)),
                                            MAE = rep(0, length(mod_vector_follow)), 
                                            Rsquared = rep(0, length(mod_vector_follow)))

rs1_follow_cv <- train(Incidence~Treatment+Triplet+log(Hist3yr)+log(Baseline), family = poisson, data=rbctconf_follow, method = "glm",
                       trControl = train.control)
cv_mods_follow[which(mod == "rs1_follow"), RMSE:= rs1_follow_cv$results$RMSE]
cv_mods_follow[which(mod == "rs1_follow"), MAE:= rs1_follow_cv$results$MAE]


mod_vector_after <- c("rs1_after", "rs2_after", "rs3_after", "rs3a_after", "rs4_after",
                      "rs4a_after", "rs4c_after", "rs5_after", "rs5a_after")

cv_mods_after <- data.table(mod = mod_vector_after, 
                                           RMSE = rep(0, length(mod_vector_after)),
                                           MAE = rep(0, length(mod_vector_after)), 
                                           Rsquared = rep(0, length(mod_vector_after)))

rs1_after_cv <- train(Incidence~Treatment+Triplet+log(Hist3yr)+log(Baseline), family = poisson, data=rbctconf_after, method = "glm",
                      trControl = train.control)
cv_mods_after[which(mod == "rs1_after"), RMSE:= rs1_after_cv$results$RMSE]
cv_mods_after[which(mod == "rs1_after"), MAE:= rs1_after_cv$results$MAE]


#DHARMA Simulated Residuals Plots
plot(simulateResiduals(rs1))
plot(simulateResiduals(rs1_follow))
plot(simulateResiduals(rs1_after))


#Information Criteria
BIC_rs1 <- BIC(rs1)
AICc_rs1 <- AICc(rs1)

BIC_rs1_follow <- BIC(rs1_follow)
AICc_rs1_follow <- AICc(rs1_follow)

BIC_rs1_after <- BIC(rs1_after)
AICc_rs1_after <- AICc(rs1_after)

#Model 2-----
# As before, with quasi-Poisson model structure
rs2<-glm(Incidence~Treatment+Triplet+log(Hist3yr)+offset(log(Baseline)), family = quasipoisson, data=rbctconf)
check_model(rs2)

rs2_follow<-glm(Incidence~Treatment+Triplet+log(Hist3yr)+offset(log(Baseline)), family = quasipoisson, data=rbctconf_follow)
summary(rs2_follow)
check_model(rs2_follow)

rs2_after<-glm(Incidence~Treatment+Triplet+log(Hist3yr)+offset(log(Baseline)), family = quasipoisson, data=rbctconf_after)
summary(rs2_after)
check_model(rs2_after)

rs2_ovr<-glm(Incidence~Treatment+Triplet+log(Hist3yr)+offset(log(Baseline)), family = quasipoisson, data=rbctconf_ovr)
check_model(rs2_ovr)
summary(rs2_ovr)

#LOOCV procedure
rs2_cv <- train(Incidence~Treatment+Triplet+log(Hist3yr)+offset(log(Baseline)), family = quasipoisson, data=rbctconf, method = "glm",
                trControl = train.control)
cv_mods[which(mod == "rs2"), RMSE:= rs2_cv$results$RMSE]
cv_mods[which(mod == "rs2"), MAE:= rs2_cv$results$MAE]


rs2_follow_cv <- train(Incidence~Treatment+Triplet+log(Hist3yr)+offset(log(Baseline)), family = quasipoisson, data=rbctconf_follow, method = "glm",
                       trControl = train.control)
cv_mods_follow[which(mod == "rs2_follow"), RMSE:= rs2_follow_cv$results$RMSE]
cv_mods_follow[which(mod == "rs2_follow"), MAE:= rs2_follow_cv$results$MAE]


rs2_after_cv <- train(Incidence~Treatment+Triplet+log(Hist3yr)+offset(log(Baseline)), family = quasipoisson, data=rbctconf_after, method = "glm",
                      trControl = train.control)
cv_mods_after[which(mod == "rs2_after"), RMSE:= rs2_after_cv$results$RMSE]
cv_mods_after[which(mod == "rs2_after"), MAE:= rs2_after_cv$results$MAE]

#SE's for estimating 95% CI's of treatment/culling effects
rs2_se <- sqrt(diag(vcov(rs2)))
exp(rs2$coefficients[2])-1
exp(rs2$coefficients[2]-qnorm(0.975) * rs2_se[2])-1
exp(rs2$coefficients[2]+qnorm(0.975) * rs2_se[2])-1

rs2_follow_se <- sqrt(diag(vcov(rs2_follow)))
exp(rs2_follow$coefficients[2])-1
exp(rs2_follow$coefficients[2]-qnorm(0.975) * rs2_follow_se[2])-1
exp(rs2_follow$coefficients[2]+qnorm(0.975) * rs2_follow_se[2])-1

rs2_after_se <- sqrt(diag(vcov(rs2_after)))
exp(rs2_after$coefficients[2])-1
exp(rs2_after$coefficients[2]-qnorm(0.975) * rs2_after_se[2])-1
exp(rs2_after$coefficients[2]+qnorm(0.975) * rs2_after_se[2])-1

#Influential Observations Investigation
#1 Observation identified as influential using Cook's distance
influence.measures(rs2)$is.inf
colSums( influence.measures(rs2)$is.inf )
colSums( influence.measures(rs2_follow)$is.inf )
colSums( influence.measures(rs2_after)$is.inf )
colSums( influence.measures(rs1)$is.inf ) #None for "original" Poisson model of Nature (2006)

pdf(paste0(file.path(out.dir2, "rs2_influence_plot.pdf")))
check_model(rs2, check = "outliers", main = "Influence Plots")
influencePlot(rs2)
dev.off()

check_model(rs2_follow, check = "outliers", main = "Influence Plots")
influencePlot(rs2_follow)



#Model 3----
#Original model of Nature (2006) in generalised Poisson form
rs3 <- glmmTMB(Incidence~Treatment+Triplet+log(Hist3yr)+log(Baseline),family=genpois, data=rbctconf)
check_model(rs3) 

rs3_follow<-glmmTMB(Incidence~Treatment+Triplet+log(Hist3yr)+log(Baseline),family=genpois, data=rbctconf_follow)
check_model(rs3_follow)
rs3_after<-glmmTMB(Incidence~Treatment+Triplet+log(Hist3yr)+log(Baseline),family=genpois, data=rbctconf_after)
summary(rs3_after)
check_model(rs3_after)

pdf("dharma_example.pdf")
plot(simulateResiduals(rs3_follow))
dev.off()


pdf(paste0(file.path(out.dir2, "rs3_dharma_plot.pdf")), height = 12, width = 20,
    pointsize = 30)
plot(simulateResiduals(rs3))

dev.off()

pdf(paste0(file.path(out.dir2, "rs3_after_dharma_plot.pdf")), height = 12, width = 20,
    pointsize = 30)
plot(simulateResiduals(rs3_after))
dev.off()

pdf(paste0(file.path(out.dir2, "rs3_follow_dharma_plot.pdf")), height = 12, width = 20,
    pointsize = 30)
plot(simulateResiduals(rs3_follow))
dev.off()

#SE's for 95% CI's for treatment effects in each period of analysis 
rs3_se <- summary(rs3$sdr, "fixed")[2,2]
exp(summary(rs3$sdr, "fixed")[2,1])-1
exp(summary(rs3$sdr, "fixed")[2,1]+qnorm(0.975) * rs3_se)-1
exp(summary(rs3$sdr, "fixed")[2,1]-qnorm(0.975) * rs3_se)-1

rs3_follow_se <- summary(rs3_follow$sdr, "fixed")[2,2]
exp(summary(rs3_follow$sdr, "fixed")[2,1])-1
exp(summary(rs3_follow$sdr, "fixed")[2,1]-qnorm(0.975) * rs3_follow_se)-1
exp(summary(rs3_follow$sdr, "fixed")[2,1]+qnorm(0.975) * rs3_follow_se)-1

rs3_after_se <- summary(rs3_after$sdr, "fixed")[2,2]
exp(summary(rs3_after$sdr, "fixed")[2,1])-1
exp(summary(rs3_after$sdr, "fixed")[2,1]-qnorm(0.975) * rs3_after_se)-1
exp(summary(rs3_after$sdr, "fixed")[2,1]+qnorm(0.975) * rs3_after_se)-1

#Information Criteria
BIC_rs3 <- BIC(rs3)
AICc_rs3 <- AICc(rs3)

BIC_rs3_follow <- BIC(rs3_follow)
AICc_rs3_follow <- AICc(rs3_follow)

BIC_rs3_after <- BIC(rs3_after)
AICc_rs3_after <- AICc(rs3_after)

#Need model-specific LOOCV as caret::train doesn't work with glmmTMB
LOOCV <- lapply(1:nrow(rbctconf), function(x){
  m1 <- glmmTMB(Incidence~Treatment+Triplet+log(Hist3yr)+log(Baseline),
                family=genpois, data=rbctconf[-x,])
  return(predict(m1, rbctconf[x,], type = "response"))
})

#Squared Errors 
sq_errors <- rep(0, length(LOOCV))
abs_errors <- rep(0, length(LOOCV))

for(i in 1:length(LOOCV)){
  sq_errors[i] <- (LOOCV[[i]]-rbctconf$Incidence[i])^2
  abs_errors[i] <- abs(LOOCV[[i]]-rbctconf$Incidence[i])
}
cv_mods[which(mod == "rs3"), RMSE:= sqrt(mean(sq_errors))]
cv_mods[which(mod == "rs3"), MAE:= mean(abs_errors)]



LOOCV_3_follow <- lapply(1:nrow(rbctconf_follow), function(x){
  m1_3_follow <- glmmTMB(Incidence~Treatment+Triplet+log(Hist3yr)+log(Baseline),
                         family=genpois, data=rbctconf_follow[-x,])
  return(predict(m1_3_follow, rbctconf_follow[x,], type = "response"))
})
sq_errors_3_follow <- rep(0, length(LOOCV_3_follow))
abs_errors_3_follow <- rep(0, length(LOOCV_3_follow))
for(i in 1:length(LOOCV_3_follow)){
  sq_errors_3_follow[i] <- (LOOCV_3_follow[[i]]-rbctconf_follow$Incidence[i])^2
  abs_errors_3_follow[i] <- abs(LOOCV_3_follow[[i]]-rbctconf_follow$Incidence[i])
}
cv_mods_follow[which(mod == "rs3_follow"), RMSE:= sqrt(mean(sq_errors_3_follow))]
cv_mods_follow[which(mod == "rs3_follow"), MAE:= mean(abs_errors_3_follow)]


LOOCV_3_after <- lapply(1:nrow(rbctconf_after), function(x){
  m1_3_after <- glmmTMB(Incidence~Treatment+Triplet+log(Hist3yr)+ log(Baseline),
                        family=genpois, data=rbctconf_after[-x,])
  return(predict(m1_3_after, rbctconf_after[x,], type = "response"))
})
sq_errors_3_after <- rep(0, length(LOOCV_3_after))
abs_errors_3_after <- rep(0, length(LOOCV_3_after))
for(i in 1:length(LOOCV_3_after)){
  sq_errors_3_after[i] <- (LOOCV_3_after[[i]]-rbctconf_after$Incidence[i])^2
  abs_errors_3_after[i] <- abs(LOOCV_3_after[[i]]-rbctconf_after$Incidence[i])
}
cv_mods_after[which(mod == "rs3_after"), RMSE:= sqrt(mean(sq_errors_3_after))]
cv_mods_after[which(mod == "rs3_after"), MAE:= mean(abs_errors_3_after)]
cv_mods_after




#Model 3a----
rs3a<- glmmTMB(Incidence~log(Hist3yr)+log(Baseline),family=genpois, data=rbctconf)
summary(rs3a)
plot(simulateResiduals(rs3a))
check_model(rs3a) 


#Posterior Predictive check indicates some potential misfit
rs3a_follow<-glmmTMB(Incidence~log(Hist3yr)+log(Baseline),family=genpois, data=rbctconf_follow)
summary(rs3a_follow)
plot(simulateResiduals(rs3a_follow))
check_model(rs3a_follow)

rs3a_after<-glmmTMB(Incidence~log(Hist3yr)+log(Baseline),family=genpois, data=rbctconf_after)
rs3a_ovr<-glmmTMB(Incidence~log(Hist3yr)+log(Baseline),family=genpois, data=rbctconf_ovr)

summary(rs3a_after)
plot(simulateResiduals(rs3a_after))
check_model(rs3a_after)
check_model(rs3a_ovr)


BIC_rs3a <- BIC(rs3a)
BIC_rs3a
AICc_rs3a <- AICc(rs3a)
AICc_rs3a

BIC_rs3a_follow <- BIC(rs3a_follow)
BIC_rs3a_follow
AICc_rs3a_follow <- AICc(rs3a_follow)
AICc_rs3a_follow

BIC_rs3a_after <- BIC(rs3a_after)
BIC_rs3a_after
AICc_rs3a_after <- AICc(rs3a_after)
AICc_rs3a_after


LOOCV_3a <- lapply(1:nrow(rbctconf), function(x){
  m3a <- glmmTMB(Incidence~log(Hist3yr)+log(Baseline),
                 family=genpois, data=rbctconf[-x,])
  return(predict(m3a, rbctconf[x,], type = "response"))
})
sq_errors_3a <- rep(0, length(LOOCV_3a))
abs_errors_3a <- rep(0, length(LOOCV_3a))

for(i in 1:length(LOOCV_3a)){
  sq_errors_3a[i] <- (LOOCV_3a[[i]]-rbctconf$Incidence[i])^2
  abs_errors_3a[i] <- abs(LOOCV_3a[[i]]-rbctconf$Incidence[i])
}
cv_mods[which(mod == "rs3a"), RMSE:= sqrt(mean(sq_errors_3a))]
cv_mods[which(mod == "rs3a"), MAE:= mean(abs_errors_3a)]


LOOCV_3a_follow <- lapply(1:nrow(rbctconf_follow), function(x){
  m3a_follow <- glmmTMB(Incidence~log(Hist3yr)+log(Baseline),
                        family=genpois, data=rbctconf_follow[-x,])
  return(predict(m3a_follow, rbctconf_follow[x,], type = "response"))
})
sq_errors_3a_follow <- rep(0, length(LOOCV_3a_follow))
abs_errors_3a_follow <- rep(0, length(LOOCV_3a_follow))

for(i in 1:length(LOOCV_3a_follow)){
  sq_errors_3a_follow[i] <- (LOOCV_3a_follow[[i]]-rbctconf_follow$Incidence[i])^2
  abs_errors_3a_follow[i] <- abs(LOOCV_3a_follow[[i]]-rbctconf_follow$Incidence[i])
}
cv_mods_follow[which(mod == "rs3a_follow"), RMSE:= sqrt(mean(sq_errors_3a_follow))]
cv_mods_follow[which(mod == "rs3a_follow"), MAE:= mean(abs_errors_3a_follow)]
cv_mods_follow






LOOCV_3a_after <- lapply(1:nrow(rbctconf_after), function(x){
  m3a_after <- glmmTMB(Incidence~log(Hist3yr)+log(Baseline),
                       family=genpois, data=rbctconf_after[-x,])
  return(predict(m3a_after, rbctconf_after[x,], type = "response"))
})
sq_errors_3a_after <- rep(0, length(LOOCV_3a_after))
abs_errors_3a_after <- rep(0, length(LOOCV_3a_after))

for(i in 1:length(LOOCV_3a_after)){
  sq_errors_3a_after[i] <- (LOOCV_3a_after[[i]]-rbctconf_after$Incidence[i])^2
  abs_errors_3a_after[i] <- abs(LOOCV_3a_after[[i]]-rbctconf_after$Incidence[i])
}
cv_mods_after[which(mod == "rs3a_after"), RMSE:= sqrt(mean(sq_errors_3a_after))]
cv_mods_after[which(mod == "rs3a_after"), MAE:= mean(abs_errors_3a_after)]
cv_mods_after



LOOCV_3a_ovr <- lapply(1:nrow(rbctconf_ovr), function(x){
  m3a_ovr <- glmmTMB(Incidence~log(Hist3yr)+log(Baseline),
                     family=genpois, data=rbctconf_ovr[-x,])
  return(predict(m3a_ovr, rbctconf_ovr[x,], type = "response"))
})
sq_errors_3a_ovr <- rep(0, length(LOOCV_3a_ovr))
abs_errors_3a_ovr <- rep(0, length(LOOCV_3a_ovr))

for(i in 1:length(LOOCV_3a_ovr)){
  sq_errors_3a_ovr[i] <- (LOOCV_3a_ovr[[i]]-rbctconf_ovr$Incidence[i])^2
  abs_errors_3a_ovr[i] <- abs(LOOCV_3a_ovr[[i]]-rbctconf_ovr$Incidence[i])
}
cv_mods_ovr[which(mod == "rs3a_ovr"), RMSE:= sqrt(mean(sq_errors_3a_ovr))]
cv_mods_ovr[which(mod == "rs3a_ovr"), MAE:= mean(abs_errors_3a_ovr)]
cv_mods_ovr

rs3c_follow<-glmmTMB(Incidence~Treatment + log(Hist3yr)+offset(log(Baseline))+Triplet,family=genpois, data=rbctconf_follow)
summary(rs3c_follow)
AICc(rs3a_follow)
summary(rs2_follow)
?glmmTMB

rs3b<-glmmTMB(Incidence~Treatment + Triplet + offset(log(Baseline)),family=genpois, data=rbctconf)
summary(rs3b)

testResiduals(rs3a)
plot(simulateResiduals(rs3c_follow, n = 1000))
plot(simulateResiduals(rs3, n = 1000))

rs3a_follow<-glmmTMB(Incidence~log(Hist3yr)+log(Baseline),family=genpois, data=rbctconf_follow)
summary(rs3a_follow)
testResiduals(rs3a_follow)
plot(simulateResiduals(rs3a_follow, n = 1000))
AIC(rs3a_follow)
AIC(rs3_follow)
glmmTMB::diagnose(rs3a)
glmmTMB::diagnose(rs3)



#Model 4----
rs4<-glmmTMB(Incidence~Treatment+Triplet+
               log(Hist3yr)+offset(log(hdyrsrisk)), family=genpois,data=rbctconf)
summary(rs4)
AICc(rs4)
check_model(rs4)
rs4_se <- summary(rs4$sdr, "fixed")[2,2]
summary(rs4$sdr, "fixed")[2,1]
exp(summary(rs4$sdr, "fixed")[2,1])-1

exp(summary(rs4$sdr, "fixed")[2,1]+qnorm(0.975) * rs4_se)-1
exp(summary(rs4$sdr, "fixed")[2,1]-qnorm(0.975) * rs4_se)-1

LOOCV_4 <- lapply(1:nrow(rbctconf), function(x){
  m4 <- glmmTMB(Incidence~Treatment+Triplet+
                  log(Hist3yr)+offset(log(hdyrsrisk)), 
                family=genpois, data=rbctconf[-x,])
  return(predict(m4, rbctconf[x,], type = "response"))
})
sq_errors_4 <- rep(0, length(LOOCV_4))
abs_errors_4 <- rep(0, length(LOOCV_4))

for(i in 1:length(LOOCV_4)){
  sq_errors_4[i] <- (LOOCV_4[[i]]-rbctconf$Incidence[i])^2
  abs_errors_4[i] <- abs(LOOCV_4[[i]]-rbctconf$Incidence[i])
}
cv_mods[which(mod == "rs4"), RMSE:= sqrt(mean(sq_errors_4))]
cv_mods[which(mod == "rs4"), MAE:= mean(abs_errors_4)]
cv_mods



plot(simulateResiduals(rs4, n = 2500, quantreg = TRUE))



rs4_follow<-glmmTMB(Incidence~Treatment+Triplet+
                      log(Hist3yr)+offset(log(hdyrsrisk)), family=genpois,data=rbctconf_follow)
summary(rs4_follow)
AICc(rs4_follow)
rs4_follow_se <- summary(rs4_follow$sdr, "fixed")[2,2]
exp(summary(rs4_follow$sdr, "fixed")[2,1])-1
exp(summary(rs4_follow$sdr, "fixed")[2,1]+qnorm(0.975) * rs4_follow_se)-1
exp(summary(rs4_follow$sdr, "fixed")[2,1]-qnorm(0.975) * rs4_follow_se)-1

plot(simulateResiduals(rs4_follow))
LOOCV_4_follow <- lapply(1:nrow(rbctconf), function(x){
  m4_follow <- glmmTMB(Incidence~Treatment+Triplet+
                         log(Hist3yr)+offset(log(hdyrsrisk)), 
                       family=genpois, data=rbctconf_follow[-x,])
  return(predict(m4_follow, rbctconf_follow[x,], type = "response"))
})
sq_errors_4_follow <- rep(0, length(LOOCV_4_follow))
abs_errors_4_follow <- rep(0, length(LOOCV_4_follow))

for(i in 1:length(LOOCV_4_follow)){
  sq_errors_4_follow[i] <- (LOOCV_4_follow[[i]]-rbctconf_follow$Incidence[i])^2
  abs_errors_4_follow[i] <- abs(LOOCV_4_follow[[i]]-rbctconf_follow$Incidence[i])
}
cv_mods_follow[which(mod == "rs4_follow"), RMSE:= sqrt(mean(sq_errors_4_follow))]
cv_mods_follow[which(mod == "rs4_follow"), MAE:= mean(abs_errors_4_follow)]
cv_mods_follow





rs4_after<-glmmTMB(Incidence~Treatment+Triplet+
                     log(Hist3yr)+offset(log(hdyrsrisk)), family=genpois,data=rbctconf_after)
summary(rs4_after)
AICc(rs4_after)
check_model(rs4_after)
exp(summary(rs4_after$sdr, "fixed")[2,1])-1
rs4_after_se <- summary(rs4_after$sdr, "fixed")[2,2]
exp(summary(rs4_after$sdr, "fixed")[2,1]+qnorm(0.975) * rs4_after_se)-1
exp(summary(rs4_after$sdr, "fixed")[2,1]-qnorm(0.975) * rs4_after_se)-1

plot(simulateResiduals(rs4_after))

LOOCV_4_after <- lapply(1:nrow(rbctconf), function(x){
  m4_after <- glmmTMB(Incidence~Treatment+Triplet+
                        log(Hist3yr)+offset(log(hdyrsrisk)), 
                      family=genpois, data=rbctconf_after[-x,])
  return(predict(m4_after, rbctconf_after[x,], type = "response"))
})
sq_errors_4_after <- rep(0, length(LOOCV_4_after))
abs_errors_4_after <- rep(0, length(LOOCV_4_after))

for(i in 1:length(LOOCV_4_after)){
  sq_errors_4_after[i] <- (LOOCV_4_after[[i]]-rbctconf_after$Incidence[i])^2
  abs_errors_4_after[i] <- abs(LOOCV_4_after[[i]]-rbctconf_after$Incidence[i])
}
cv_mods_after[which(mod == "rs4_after"), RMSE:= sqrt(mean(sq_errors_4_after))]
cv_mods_after[which(mod == "rs4_after"), MAE:= mean(abs_errors_4_after)]
cv_mods_after




BIC_rs4 <- BIC(rs4)
BIC_rs4
AICc_rs4 <- AICc(rs4)
AICc_rs4

BIC_rs4_follow <- BIC(rs4_follow)
BIC_rs4_follow
AICc_rs4_follow <- AICc(rs4_follow)
AICc_rs4_follow

BIC_rs4_after <- BIC(rs4_after)
BIC_rs4_after
AICc_rs4_after <- AICc(rs4_after)
AICc_rs4_after



#Model 4a----
rs4a<-glmmTMB(Incidence~log(Hist3yr)+offset(log(hdyrsrisk)), family=genpois,data=rbctconf)
summary(rs4a)
AICc(rs4a)
check_model(rs4a)
plot(simulateResiduals(rs4a))
head(rbctconf)


plot(simulateResiduals(rs4a, quantreg = TRUE))

LOOCV_4a <- lapply(1:nrow(rbctconf), function(x){
  m4a <- glmmTMB(Incidence~log(Hist3yr)+offset(log(hdyrsrisk)),
                 family=genpois, data=rbctconf[-x,])
  return(predict(m4a, rbctconf[x,], type = "response"))
})
sq_errors_4a <- rep(0, length(LOOCV_4a))
abs_errors_4a <- rep(0, length(LOOCV_4a))

for(i in 1:length(LOOCV_4a)){
  sq_errors_4a[i] <- (LOOCV_4a[[i]]-rbctconf$Incidence[i])^2
  abs_errors_4a[i] <- abs(LOOCV_4a[[i]]-rbctconf$Incidence[i])
}
cv_mods[which(mod == "rs4a"), RMSE:= sqrt(mean(sq_errors_4a))]
cv_mods[which(mod == "rs4a"), MAE:= mean(abs_errors_4a)]
cv_mods



rs4a_follow<-glmmTMB(Incidence~log(Hist3yr)+offset(log(hdyrsrisk)), 
                     family=genpois,data=rbctconf_follow)
summary(rs4a_follow)
AICc(rs4a_follow)
check_model(rs4a_follow)
plot(simulateResiduals(rs4a_follow))
head(rbctconf)
plot(simulateResiduals(rs4a_follow, quantreg = TRUE))

LOOCV_4a_follow <- lapply(1:nrow(rbctconf_follow), function(x){
  m4a_follow <- glmmTMB(Incidence~log(Hist3yr)+offset(log(hdyrsrisk)),
                        family=genpois, data=rbctconf_follow[-x,])
  return(predict(m4a_follow, rbctconf_follow[x,], type = "response"))
})
sq_errors_4a_follow <- rep(0, length(LOOCV_4a_follow))
abs_errors_4a_follow <- rep(0, length(LOOCV_4a_follow))

for(i in 1:length(LOOCV_4a_follow)){
  sq_errors_4a_follow[i] <- (LOOCV_4a_follow[[i]]-rbctconf_follow$Incidence[i])^2
  abs_errors_4a_follow[i] <- abs(LOOCV_4a_follow[[i]]-rbctconf_follow$Incidence[i])
}
cv_mods_follow[which(mod == "rs4a_follow"), RMSE:= sqrt(mean(sq_errors_4a))]
cv_mods_follow[which(mod == "rs4a_follow"), MAE:= mean(abs_errors_4a)]
cv_mods_follow

cv_mods_follow[order(MAE, decreasing = TRUE)]



rs4a_after<-glmmTMB(Incidence~log(Hist3yr)+offset(log(hdyrsrisk)), 
                    family=genpois,data=rbctconf_after)
plot(simulateResiduals(rs4a_after, quantreg = TRUE))
check_model(rs4a_after)
LOOCV_4a_after <- lapply(1:nrow(rbctconf_after), function(x){
  m4a_after <- glmmTMB(Incidence~log(Hist3yr)+offset(log(hdyrsrisk)),
                       family=genpois, data=rbctconf_after[-x,])
  return(predict(m4a_after, rbctconf_after[x,], type = "response"))
})
sq_errors_4a_after <- rep(0, length(LOOCV_4a_after))
abs_errors_4a_after <- rep(0, length(LOOCV_4a_after))

for(i in 1:length(LOOCV_4a_after)){
  sq_errors_4a_after[i] <- (LOOCV_4a_after[[i]]-rbctconf_after$Incidence[i])^2
  abs_errors_4a_after[i] <- abs(LOOCV_4a_after[[i]]-rbctconf_after$Incidence[i])
}

cv_mods_after[which(mod == "rs4a_after"), RMSE:= sqrt(mean(sq_errors_4a_after))]
cv_mods_after[which(mod == "rs4a_after"), MAE:= mean(abs_errors_4a_after)]


BIC_rs4a <- BIC(rs4a)
BIC_rs4a
AICc_rs4a <- AICc(rs4a)
AICc_rs4a

BIC_rs4a_follow <- BIC(rs4a_follow)
BIC_rs4a_follow
AICc_rs4a_follow <- AICc(rs4a_follow)
AICc_rs4a_follow

BIC_rs4a_after <- BIC(rs4a_after)
BIC_rs4a_after
AICc_rs4a_after <- AICc(rs4a_after)
AICc_rs4a_after


pdf(paste0(file.path(out.dir2, "rs4a_pp_check.pdf")))
check_model(rs4a, check = "pp_check")
dev.off()


pdf(paste0(file.path(out.dir2, "rs4a_after_pp_check.pdf")))
check_model(rs4a_after, check = "pp_check")
dev.off()
par(mfrow = c(1,1))
rs3


#Frequentist GLM Plotting Diagnostics ----
rs2_diagnostics <- check_model(rs2, check = "outliers")
rs2_diag_plots <- plot(rs2_diagnostics, return_list = TRUE)

grid.arrange(rs2_diag_plots+ theme(text = element_text(size = 22),axis.text.x = element_text(size=20),
                                   axis.text.y = element_text(size=20),
                                   plot.title = element_text(size=24, face = "bold", hjust = 0.5),
                                   plot.subtitle = element_text(size=20, hjust = 0.5),
                                   axis.title=element_text(size=20), 
                                   legend.text=element_text(size=18))+ labs(title = "", subtitle = ""),
             nrow = 1)
rs2_diag_plots_grob <- arrangeGrob(rs2_diag_plots+ theme(text = element_text(size = 28),
                                                         axis.text.x = element_text(size=28),
                                                         axis.text.y = element_text(size=28),
                                                         plot.title = element_text(size=24, face = "bold", hjust = 0.5),
                                                         plot.subtitle = element_text(size=20, hjust = 0.5),
                                                         axis.title=element_text(size=26), 
                                                         legend.text=element_text(size=18))+ labs(title = "", subtitle = ""),
                                   nrow = 1)
ggsave(rs2_diag_plots_grob, file = file.path(out.dir2,"rs2_influence_plot.pdf"), h = 14, w = 20)

#rs2 follow influence
rs2_follow_diagnostics <- check_model(rs2_follow, check = "outliers")
rs2_follow_diag_plots <- plot(rs2_follow_diagnostics, return_list = TRUE)
rs2_follow_diag_plots
grid.arrange(rs2_follow_diag_plots+ theme(text = element_text(size = 26),
                                          plot.title = element_text(size=24, face = "bold", hjust = 0.5),
                                          plot.subtitle = element_text(size=20, hjust = 0.5),
                                          axis.title=element_text(size=26), 
                                          legend.text=element_text(size=18))+ labs(title = "", subtitle = ""),
             nrow = 1)
rs2_follow_diag_plots_grob <- arrangeGrob(rs2_follow_diag_plots+ theme(text = element_text(size = 28),
                                                                       axis.text.x = element_text(size=28),
                                                                       axis.text.y = element_text(size=28),
                                                                       plot.title = element_text(size=24, face = "bold", hjust = 0.5),
                                                                       plot.subtitle = element_text(size=20, hjust = 0.5),
                                                                       axis.title=element_text(size=26), 
                                                                       legend.text=element_text(size=18))+ labs(title = "", subtitle = ""),
                                          nrow = 1)
ggsave(rs2_follow_diag_plots_grob, file = file.path(out.dir2,"rs2_follow_influence_plot.pdf"), h = 14, w = 20)

#rs2 after influence
rs2_after_diagnostics <- check_model(rs2_after, check = "outliers")
rs2_after_diag_plots <- plot(rs2_after_diagnostics, return_list = TRUE)

grid.arrange(rs2_after_diag_plots+ theme(text = element_text(size = 20),
                                         axis.text.x = element_text(size=16),
                                         axis.text.y = element_text(size=16),
                                         plot.title = element_text(size=20, face = "bold", hjust = 0.5),
                                         plot.subtitle = element_text(size=16, hjust = 0.5),
                                         axis.title=element_text(size=18), 
                                         legend.text=element_text(size=14)+geom_text(size = 20))+ labs(title = "", subtitle = ""),
             nrow = 1)
rs2_after_diag_plots_grob <- arrangeGrob(rs2_after_diag_plots+ theme(text = element_text(size = 28),
                                                                     axis.text.x = element_text(size=28),
                                                                     axis.text.y = element_text(size=28),
                                                                     plot.title = element_text(size=24, face = "bold", hjust = 0.5),
                                                                     plot.subtitle = element_text(size=20, hjust = 0.5),
                                                                     axis.title=element_text(size=26), 
                                                                     legend.text=element_text(size=18))+ labs(title = "", subtitle = ""),
                                         nrow = 1)
ggsave(rs2_after_diag_plots_grob, file = file.path(out.dir2,"rs2_after_influence_plot.pdf"), h = 14, w = 20)
check_model(rs4a_follow, check = "pp_check")


#rs4a pp check
rs4a_diagnostics <- check_model(rs4a, check = "pp_check")
rs4a_diag_plots <- plot(rs4a_diagnostics, return_list = TRUE)
grid.arrange(rs4a_diag_plots+ theme(axis.text.x = element_text(size=16),
                                    axis.text.y = element_text(size=16),
                                    plot.title = element_text(size=20, face = "bold", hjust = 0.5),
                                    plot.subtitle = element_text(size=16, hjust = 0.5),
                                    axis.title=element_text(size=18), 
                                    legend.text=element_text(size=14)),
             nrow = 1)
rs4a_diag_plots_grob <- arrangeGrob(rs4a_diag_plots+ theme(text = element_text(size = 24),
                                                           axis.text.x = element_text(size=24),
                                                           axis.text.y = element_text(size=24),
                                                           plot.title = element_text(size=24, face = "bold", hjust = 0.5),
                                                           plot.subtitle = element_text(size=20, hjust = 0.5),
                                                           axis.title=element_text(size=20), 
                                                           legend.text=element_text(size=24))+ labs(title = "", subtitle = ""),
                                    nrow = 1)
ggsave(rs4a_diag_plots_grob, file = file.path(out.dir2,"rs4a_pp_check.pdf"), h = 14, w = 20)

#rs4a follow pp check---
rs4a_follow_diagnostics <- check_model(rs4a_follow, check = "pp_check")
rs4a_follow_diag_plots <- plot(rs4a_follow_diagnostics, return_list = TRUE)
grid.arrange(rs4a_follow_diag_plots+theme(text = element_text(size = 24),
                                          axis.text.x = element_text(size=24),
                                          axis.text.y = element_text(size=24),
                                          plot.title = element_text(size=24, face = "bold", hjust = 0.5),
                                          plot.subtitle = element_text(size=20, hjust = 0.5),
                                          axis.title=element_text(size=20), 
                                          legend.text=element_text(size=18))+ labs(title = "", subtitle = ""),
             nrow = 1)
rs4a_follow_diag_plots_grob <- arrangeGrob(rs4a_follow_diag_plots+ theme(text = element_text(size = 24),
                                                                         axis.text.x = element_text(size=24),
                                                                         axis.text.y = element_text(size=24),
                                                                         plot.title = element_text(size=24, face = "bold", hjust = 0.5),
                                                                         plot.subtitle = element_text(size=20, hjust = 0.5),
                                                                         axis.title=element_text(size=20), 
                                                                         legend.text=element_text(size=24))+ labs(title = "", subtitle = ""),
                                           nrow = 1)
ggsave(rs4a_follow_diag_plots_grob, file = file.path(out.dir2,"rs4a_follow_pp_check.pdf"), h = 14, w = 20)

rs4a_after_diagnostics <- check_model(rs4a_after, check = "pp_check")
rs4a_after_diag_plots <- plot(rs4a_after_diagnostics, return_list = TRUE)
grid.arrange(rs4a_after_diag_plots+ theme(axis.text.x = element_text(size=16),
                                          axis.text.y = element_text(size=16),
                                          plot.title = element_text(size=20, face = "bold", hjust = 0.5),
                                          plot.subtitle = element_text(size=16, hjust = 0.5),
                                          axis.title=element_text(size=18), 
                                          legend.text=element_text(size=14)),
             nrow = 1)
rs4a_after_diag_plots_grob <- arrangeGrob(rs4a_after_diag_plots+ theme(text = element_text(size = 24),
                                                                       axis.text.x = element_text(size=24),
                                                                       axis.text.y = element_text(size=24),
                                                                       plot.title = element_text(size=24, face = "bold", hjust = 0.5),
                                                                       plot.subtitle = element_text(size=20, hjust = 0.5),
                                                                       axis.title=element_text(size=20), 
                                                                       legend.text=element_text(size=24))+ labs(title = "", subtitle = ""),
                                          nrow = 1)
ggsave(rs4a_after_diag_plots_grob, file = file.path(out.dir2,"rs4a_after_pp_check.pdf"), h = 14, w = 20)








rs3_diagnostics <- check_model(rs3, check = "pp_check")
rs3_diag_plots <- plot(rs3_diagnostics, return_list = TRUE)
grid.arrange(rs3_diag_plots+ theme(axis.text.x = element_text(size=16),
                                   axis.text.y = element_text(size=16),
                                   plot.title = element_text(size=20, face = "bold", hjust = 0.5),
                                   plot.subtitle = element_text(size=16, hjust = 0.5),
                                   axis.title=element_text(size=18), 
                                   legend.text=element_text(size=14)),
             nrow = 1)
rs3_diag_plots_grob <- arrangeGrob(rs3_diag_plots+ theme(text = element_text(size = 24),
                                                         axis.text.x = element_text(size=24),
                                                         axis.text.y = element_text(size=24),
                                                         plot.title = element_text(size=24, face = "bold", hjust = 0.5),
                                                         plot.subtitle = element_text(size=20, hjust = 0.5),
                                                         axis.title=element_text(size=20), 
                                                         legend.text=element_text(size=24))+ labs(title = "", subtitle = ""),
                                   nrow = 1)
dev.off()
ggsave(rs3_diag_plots_grob, file = file.path(out.dir2,"rs3_pp_check.pdf"), h = 14, w = 20)

#rs3 follow pp check----
rs3_follow_diagnostics <- check_model(rs3_follow, check = "pp_check")
rs3_follow_diag_plots <- plot(rs3_follow_diagnostics, return_list = TRUE)
grid.arrange(rs3_follow_diag_plots+ theme(axis.text.x = element_text(size=16),
                                          axis.text.y = element_text(size=16),
                                          plot.title = element_text(size=20, face = "bold", hjust = 0.5),
                                          plot.subtitle = element_text(size=16, hjust = 0.5),
                                          axis.title=element_text(size=18), 
                                          legend.text=element_text(size=14)),
             nrow = 1)
rs3_follow_diag_plots_grob <- arrangeGrob(rs3_follow_diag_plots+ theme(text = element_text(size = 24),
                                                                       axis.text.x = element_text(size=24),
                                                                       axis.text.y = element_text(size=24),
                                                                       plot.title = element_text(size=24, face = "bold", hjust = 0.5),
                                                                       plot.subtitle = element_text(size=20, hjust = 0.5),
                                                                       axis.title=element_text(size=20), 
                                                                       legend.text=element_text(size=24))+ labs(title = "", subtitle = ""),
                                          nrow = 1)
ggsave(rs3_follow_diag_plots_grob, file = file.path(out.dir2,"rs3_follow_pp_check.pdf"), h = 14, w = 20)


rs3_after_diagnostics <- check_model(rs3_after, check = "pp_check")
rs3_after_diag_plots <- plot(rs3_after_diagnostics, return_list = TRUE)
grid.arrange(rs3_after_diag_plots+ theme(axis.text.x = element_text(size=16),
                                         axis.text.y = element_text(size=16),
                                         plot.title = element_text(size=20, face = "bold", hjust = 0.5),
                                         plot.subtitle = element_text(size=16, hjust = 0.5),
                                         axis.title=element_text(size=18), 
                                         legend.text=element_text(size=14)),
             nrow = 1)
rs3_after_diag_plots_grob <- arrangeGrob(rs3_after_diag_plots+ theme(text = element_text(size = 24),
                                                                     axis.text.x = element_text(size=24),
                                                                     axis.text.y = element_text(size=24),
                                                                     plot.title = element_text(size=24, face = "bold", hjust = 0.5),
                                                                     plot.subtitle = element_text(size=20, hjust = 0.5),
                                                                     axis.title=element_text(size=20), 
                                                                     legend.text=element_text(size=24))+ labs(title = "", subtitle = ""),
                                         nrow = 1)
dev.off()
ggsave(rs3_after_diag_plots_grob, file = file.path(out.dir2,"rs3_after_pp_check.pdf"), h = 14, w = 20)



colSums(influence.measures(rs1)$is.inf)
colSums(influence.measures(rs1_follow)$is.inf)
colSums(influence.measures(rs1_after)$is.inf)
influencePlot(rs1_after)
influencePlot(rs2_after)
check_model(rs1_after)
pdf(paste0(file.path(out.dir2, "rs4a_follow_pp_check.pdf")))
check_model(rs4a_follow, check = "pp_check")
dev.off()

#Model 4c----
rs4c<-glmmTMB(Incidence~Treatment+Triplet+log(Hist3yr)+log(hdyrsrisk), family=genpois,data=rbctconf)
rs4c_se <- summary(rs4c$sdr, "fixed")[2,2]
exp(summary(rs4c$sdr, "fixed")[2,1])-1
exp(summary(rs4c$sdr, "fixed")[2,1]+qnorm(0.975) * rs4c_se)-1
exp(summary(rs4c$sdr, "fixed")[2,1]-qnorm(0.975) * rs4c_se)-1
plot(simulateResiduals(rs4c, quantreg = TRUE))
check_model(rs4c)

LOOCV_4c <- lapply(1:nrow(rbctconf), function(x){
  m4c <- glmmTMB(Incidence ~ Treatment+Triplet+log(Hist3yr)+log(hdyrsrisk),
                 family=genpois, data=rbctconf[-x,])
  return(predict(m4c, rbctconf[x,], type = "response"))
})
sq_errors_4c <- rep(0, length(LOOCV_4c))
abs_errors_4c <- rep(0, length(LOOCV_4c))

for(i in 1:length(LOOCV_4c)){
  sq_errors_4c[i] <- (LOOCV_4c[[i]]-rbctconf$Incidence[i])^2
  abs_errors_4c[i] <- abs(LOOCV_4c[[i]]-rbctconf$Incidence[i])
}
cv_mods[which(mod == "rs4c"), RMSE:= sqrt(mean(sq_errors_4c))]
cv_mods[which(mod == "rs4c"), MAE:= mean(abs_errors_4c)]






rs4c_follow<-glmmTMB(Incidence~Treatment+Triplet+log(Hist3yr)+log(hdyrsrisk), family=genpois,data=rbctconf_follow)
summary(rs4c_follow)

#95% CI's
rs4c_follow_se <- summary(rs4c_follow$sdr, "fixed")[2,2]
exp(summary(rs4c_follow$sdr, "fixed")[2,1])-1
exp(summary(rs4c_follow$sdr, "fixed")[2,1]+qnorm(0.975) * rs4c_follow_se)-1
exp(summary(rs4c_follow$sdr, "fixed")[2,1]-qnorm(0.975) * rs4c_follow_se)-1

plot(simulateResiduals(rs4c_follow, quantreg = TRUE))
check_model(rs4c_follow)

LOOCV_4c_follow <- lapply(1:nrow(rbctconf_follow), function(x){
  m4c_follow <- glmmTMB(Incidence ~ Treatment+Triplet+log(Hist3yr)+log(hdyrsrisk),
                        family=genpois, data=rbctconf_follow[-x,])
  return(predict(m4c_follow, rbctconf_follow[x,], type = "response"))
})
sq_errors_4c_follow <- rep(0, length(LOOCV_4c_follow))
abs_errors_4c_follow <- rep(0, length(LOOCV_4c_follow))

for(i in 1:length(LOOCV_4c_follow)){
  sq_errors_4c_follow[i] <- (LOOCV_4c_follow[[i]]-rbctconf_follow$Incidence[i])^2
  abs_errors_4c_follow[i] <- abs(LOOCV_4c_follow[[i]]-rbctconf_follow$Incidence[i])
}
cv_mods_follow[which(mod == "rs4c_follow"), RMSE:= sqrt(mean(sq_errors_4c_follow))]
cv_mods_follow[which(mod == "rs4c_follow"), MAE:= mean(abs_errors_4c_follow)]


rs4c_after<-glmmTMB(Incidence~Treatment+Triplet+log(Hist3yr)+log(hdyrsrisk), family=genpois,data=rbctconf_after)
summary(rs4c_after)
rs4c_after_se <- summary(rs4c_after$sdr, "fixed")[2,2]

exp(summary(rs4c$sdr, "fixed")[2,1])-1
exp(summary(rs4c_follow$sdr, "fixed")[2,1])-1
exp(summary(rs4c_after$sdr, "fixed")[2,1])-1


exp(summary(rs4c_after$sdr, "fixed")[2,1]+qnorm(0.975) * rs4c_after_se)-1
exp(summary(rs4c_after$sdr, "fixed")[2,1]-qnorm(0.975) * rs4c_after_se)-1
plot(simulateResiduals(rs4c_after, quantreg = TRUE))
check_model(rs4c_after)

LOOCV_4c_after <- lapply(1:nrow(rbctconf_after), function(x){
  m4c_after <- glmmTMB(Incidence ~ Treatment+Triplet+log(Hist3yr)+log(hdyrsrisk),
                       family=genpois, data=rbctconf_after[-x,])
  return(predict(m4c_after, rbctconf_after[x,], type = "response"))
})
sq_errors_4c_after <- rep(0, length(LOOCV_4c_after))
abs_errors_4c_after <- rep(0, length(LOOCV_4c_after))

for(i in 1:length(LOOCV_4c_after)){
  sq_errors_4c_after[i] <- (LOOCV_4c_after[[i]]-rbctconf_after$Incidence[i])^2
  abs_errors_4c_after[i] <- abs(LOOCV_4c_after[[i]]-rbctconf_after$Incidence[i])
}
cv_mods_after[which(mod == "rs4c_after"), RMSE:= sqrt(mean(sq_errors_4c_after))]
cv_mods_after[which(mod == "rs4c_after"), MAE:= mean(abs_errors_4c_after)]



BIC_rs4c <- BIC(rs4c)
BIC_rs4c
AICc_rs4c <- AICc(rs4c)
AICc_rs4c

BIC_rs4c_follow <- BIC(rs4c_follow)
BIC_rs4c_follow
AICc_rs4c_follow <- AICc(rs4c_follow)
AICc_rs4c_follow

BIC_rs4c_after <- BIC(rs4c_after)
BIC_rs4c_after
AICc_rs4c_after <- AICc(rs4c_after)
AICc_rs4c_after


#Model 4d ----
rs4d<-glmmTMB(Incidence~log(Hist3yr)+log(hdyrsrisk), family=genpois,
              data=rbctconf)
summary(rs4d)
check_model(rs4d)
AICc(rs4d)
plot(simulateResiduals(rs4d))

LOOCV_4d <- lapply(1:nrow(rbctconf), function(x){
  m4d <- glmmTMB(Incidence ~ log(Hist3yr)+log(hdyrsrisk),
                 family=genpois, data=rbctconf[-x,])
  return(predict(m4d, rbctconf[x,], type = "response"))
})
sq_errors_4d <- rep(0, length(LOOCV_4d))
abs_errors_4d <- rep(0, length(LOOCV_4d))

for(i in 1:length(LOOCV_4d)){
  sq_errors_4d[i] <- (LOOCV_4d[[i]]-rbctconf$Incidence[i])^2
  abs_errors_4d[i] <- abs(LOOCV_4d[[i]]-rbctconf$Incidence[i])
}
cv_mods[which(mod == "rs4d"), RMSE:= sqrt(mean(sq_errors_4d))]
cv_mods[which(mod == "rs4d"), MAE:= mean(abs_errors_4d)]
cv_mods






rs4d_follow<-glmmTMB(Incidence~log(Hist3yr)+log(hdyrsrisk), family=genpois,
                     data=rbctconf_follow)
summary(rs4d_follow)
exp(summary(rs4d_follow$sdr, "fixed")[2,1])-1
exp(summary(rs4d_follow$sdr, "fixed")[2,1]+qnorm(0.975) * rs4d_follow_se)-1
exp(summary(rs4d_follow$sdr, "fixed")[2,1]-qnorm(0.975) * rs4d_follow_se)-1

check_model(rs4d_follow)
AICc(rs4d_follow)
plot(simulateResiduals(rs4d))

LOOCV_4d_follow <- lapply(1:nrow(rbctconf_follow), function(x){
  m4d_follow <- glmmTMB(Incidence ~ log(Hist3yr)+log(hdyrsrisk),
                        family=genpois, data=rbctconf_follow[-x,])
  return(predict(m4d_follow, rbctconf_follow[x,], type = "response"))
})
sq_errors_4d_follow <- rep(0, length(LOOCV_4d_follow))
abs_errors_4d_follow <- rep(0, length(LOOCV_4d_follow))

for(i in 1:length(LOOCV_4d_follow)){
  sq_errors_4d_follow[i] <- (LOOCV_4d_follow[[i]]-rbctconf_follow$Incidence[i])^2
  abs_errors_4d_follow[i] <- abs(LOOCV_4d_follow[[i]]-rbctconf_follow$Incidence[i])
}
cv_mods_follow[which(mod == "rs4d_follow"), RMSE:= sqrt(mean(sq_errors_4d_follow))]
cv_mods_follow[which(mod == "rs4d_follow"), MAE:= mean(abs_errors_4d_follow)]





rs4d_after<-glmmTMB(Incidence~log(Hist3yr)+log(hdyrsrisk), family=genpois,
                    data=rbctconf_after)

summary(rs4d_after)
check_model(rs4d_after)
AICc(rs4d_after)
plot(simulateResiduals(rs4d))

LOOCV_4d_after <- lapply(1:nrow(rbctconf_after), function(x){
  m4d_after <- glmmTMB(Incidence ~ log(Hist3yr)+log(hdyrsrisk),
                       family=genpois, data=rbctconf_after[-x,])
  return(predict(m4d_after, rbctconf_after[x,], type = "response"))
})
sq_errors_4d_after <- rep(0, length(LOOCV_4d_after))
abs_errors_4d_after <- rep(0, length(LOOCV_4d_after))

for(i in 1:length(LOOCV_4d_after)){
  sq_errors_4d_after[i] <- (LOOCV_4d_after[[i]]-rbctconf_after$Incidence[i])^2
  abs_errors_4d_after[i] <- abs(LOOCV_4d_after[[i]]-rbctconf_after$Incidence[i])
}
cv_mods_after[which(mod == "rs4d_after"), RMSE:= sqrt(mean(sq_errors_4d_after))]
cv_mods_after[which(mod == "rs4d_after"), MAE:= mean(abs_errors_4d_after)]


BIC_rs4d <- BIC(rs4d)
BIC_rs4d
AICc_rs4d <- AICc(rs4d)
AICc_rs4d

BIC_rs4d_follow <- BIC(rs4d_follow)
BIC_rs4d_follow
AICc_rs4d_follow <- AICc(rs4d_follow)
AICc_rs4d_follow

BIC_rs4d_after <- BIC(rs4d_after)
BIC_rs4d_after
AICc_rs4d_after <- AICc(rs4d_after)
AICc_rs4d_after



#Model 5 ----
#Linear Model
rs5<-lm(Incidence/hdyrsrisk~Treatment+Triplet+log(Hist3yr), data=rbctconf)
rs5$coefficients*mean(rbctconf$hdyrsrisk)
summary(rs5)
confint(rs5)
plot(rs5)
rbctconf
ols_test_normality(rs5)
ols_plot_resid_fit(rs5)
ols_plot_resid_hist(rs5)
AICc(rs5)
plot(simulateResiduals(rs5))
LOOCV_5 <- lapply(1:nrow(rbctconf), function(x){
  m5 <- lm(Incidence/hdyrsrisk~Treatment+Triplet+log(Hist3yr), 
           data=rbctconf[-x,])
  return(predict(m5, rbctconf[x,], type = "response"))
})
sq_errors_5 <- rep(0, length(LOOCV_5))
abs_errors_5 <- rep(0, length(LOOCV_5))
for(i in 1:length(LOOCV_5)){
  sq_errors_5[i] <- (LOOCV_5[[i]]*rbctconf$hdyrsrisk[i]-rbctconf$Incidence[i])^2
  abs_errors_5[i] <- abs(LOOCV_5[[i]]*rbctconf$hdyrsrisk[i]-rbctconf$Incidence[i])
}
cv_mods[which(mod == "rs5"), RMSE:= sqrt(mean(sq_errors_5))]
cv_mods[which(mod == "rs5"), MAE:= mean(abs_errors_5)]






rs5_follow<-lm(Incidence/hdyrsrisk~Treatment+Triplet+log(Hist3yr), data=rbctconf_follow)
summary(rs5_follow)
confint(rs5_follow)
AICc(rs5_follow)
ols_test_normality(rs5_follow)
ols_plot_resid_fit(rs5_follow)
ols_plot_resid_hist(rs5_follow)
colSums( influence.measures(rs5_follow)$is.inf )
plot(rs5_follow)
LOOCV_5_follow <- lapply(1:nrow(rbctconf_follow), function(x){
  m5_follow <- lm(Incidence/hdyrsrisk~Treatment+Triplet+log(Hist3yr), 
                  data=rbctconf_follow[-x,])
  return(predict(m5_follow, rbctconf_follow[x,], type = "response"))
})
sq_errors_5_follow <- rep(0, length(LOOCV_5_follow))
abs_errors_5_follow <- rep(0, length(LOOCV_5_follow))
LOOCV_5_follow[[i]]
for(i in 1:length(LOOCV_5_follow)){
  sq_errors_5_follow[i] <- (LOOCV_5_follow[[i]]*rbctconf_follow$hdyrsrisk[i]-rbctconf_follow$Incidence[i])^2
  abs_errors_5_follow[i] <- abs(LOOCV_5_follow[[i]]*rbctconf_follow$hdyrsrisk[i]-rbctconf_follow$Incidence[i])
}
cv_mods_follow[which(mod == "rs5_follow"), RMSE:= sqrt(mean(sq_errors_5_follow))]
cv_mods_follow[which(mod == "rs5_follow"), MAE:= mean(abs_errors_5_follow)]




rs5_after<-lm(Incidence/hdyrsrisk~Treatment+Triplet+log(Hist3yr), data=rbctconf_after)
summary(rs5_after)
confint(rs5_after)
AICc(rs5_after)
ols_test_normality(rs5_after)
ols_plot_resid_fit(rs5_after)
ols_plot_resid_hist(rs5_after)
colSums( influence.measures(rs5_after)$is.inf )
plot(rs5_after)
LOOCV_5_after <- lapply(1:nrow(rbctconf_after), function(x){
  m5_after <- lm(Incidence/hdyrsrisk~Treatment+Triplet+log(Hist3yr), 
                 data=rbctconf_after[-x,])
  return(predict(m5_after, rbctconf_after[x,], type = "response"))
})
sq_errors_5_after <- rep(0, length(LOOCV_5_after))
abs_errors_5_after <- rep(0, length(LOOCV_5_after))
LOOCV_5_after[[i]]
for(i in 1:length(LOOCV_5_after)){
  sq_errors_5_after[i] <- (LOOCV_5_after[[i]]*rbctconf_after$hdyrsrisk[i]-rbctconf_after$Incidence[i])^2
  abs_errors_5_after[i] <- abs(LOOCV_5_after[[i]]*rbctconf_after$hdyrsrisk[i]-rbctconf_after$Incidence[i])
}
cv_mods_after[which(mod == "rs5_after"), RMSE:= sqrt(mean(sq_errors_5_after))]
cv_mods_after[which(mod == "rs5_after"), MAE:= mean(abs_errors_5_after)]


BIC_rs5 <- BIC(rs5)
BIC_rs5
AICc_rs5 <- AICc(rs5)
AICc_rs5

BIC_rs5_follow <- BIC(rs5_follow)
BIC_rs5_follow
AICc_rs5_follow <- AICc(rs5_follow)
AICc_rs5_follow

BIC_rs5_after <- BIC(rs5_after)
BIC_rs5_after
AICc_rs5_after <- AICc(rs5_after)
AICc_rs5_after



#Model 5a ----
rs5a<-lm(Incidence/hdyrsrisk~log(Hist3yr), data=rbctconf)
summary(rs5a)
ols_test_normality(rs5a)
ols_plot_resid_fit(rs5a)
ols_plot_resid_hist(rs5a)
AICc(rs5a)
plot(rs5a)
plot(simulateResiduals(rs5a))
influence.measures(rs5a)$is.inf
colSums( influence.measures(rs5a)$is.inf )

LOOCV_5a <- lapply(1:nrow(rbctconf), function(x){
  m5a <- lm(Incidence/hdyrsrisk~log(Hist3yr), 
            data=rbctconf[-x,])
  return(predict(m5a, rbctconf[x,], type = "response"))
})
sq_errors_5a <- rep(0, length(LOOCV_5a))
abs_errors_5a <- rep(0, length(LOOCV_5a))
LOOCV_5a[[i]]
for(i in 1:length(LOOCV_5a)){
  sq_errors_5a[i] <- (LOOCV_5a[[i]]*rbctconf$hdyrsrisk[i]-rbctconf$Incidence[i])^2
  abs_errors_5a[i] <- abs(LOOCV_5a[[i]]*rbctconf$hdyrsrisk[i]-rbctconf$Incidence[i])
}

cv_mods[which(mod == "rs5a"), RMSE:= sqrt(mean(sq_errors_5a))]
cv_mods[which(mod == "rs5a"), MAE:= mean(abs_errors_5a)]







rs5a_follow<-lm(Incidence/hdyrsrisk~log(Hist3yr), data=rbctconf_follow)
summary(rs5a_follow)

ols_test_normality(rs5a_follow)
ols_plot_resid_fit(rs5a_follow)
ols_plot_resid_hist(rs5a_follow)
AICc(rs5a_follow)
plot(rs5a_follow)
plot(simulateResiduals(rs5a_follow))
colSums( influence.measures(rs5a_follow)$is.inf )
influence.measures(rs5a_follow)
LOOCV_5a_follow <- lapply(1:nrow(rbctconf_follow), function(x){
  m5a_follow <- lm(Incidence/hdyrsrisk~log(Hist3yr), 
                   data=rbctconf_follow[-x,])
  return(predict(m5a_follow, rbctconf_follow[x,], type = "response"))
})
sq_errors_5a_follow <- rep(0, length(LOOCV_5a_follow))
abs_errors_5a_follow <- rep(0, length(LOOCV_5a_follow))
LOOCV_5a_follow[[i]]
for(i in 1:length(LOOCV_5a_follow)){
  sq_errors_5a_follow[i] <- (LOOCV_5a_follow[[i]]*rbctconf_follow$hdyrsrisk[i]-rbctconf_follow$Incidence[i])^2
  abs_errors_5a_follow[i] <- abs(LOOCV_5a_follow[[i]]*rbctconf_follow$hdyrsrisk[i]-rbctconf_follow$Incidence[i])
}
cv_mods_follow[which(mod == "rs5a_follow"), RMSE:= sqrt(mean(sq_errors_5a_follow))]
cv_mods_follow[which(mod == "rs5a_follow"), MAE:= mean(abs_errors_5a_follow)]








rs5a_after<-lm(Incidence/hdyrsrisk~log(Hist3yr), data=rbctconf_after)
summary(rs5a_after)
ols_test_normality(rs5a_after)
ols_plot_resid_fit(rs5a_after)
ols_plot_resid_hist(rs5a_after)
AICc(rs5a_after)
plot(rs5a_after)
plot(simulateResiduals(rs5a_after))
colSums(influence.measures(rs5a_after)$is.inf )
influence.measures(rs5a_after)
rbctconf_after[11, ]
LOOCV_5a_after <- lapply(1:nrow(rbctconf_after), function(x){
  m5a_after <- lm(Incidence/hdyrsrisk~log(Hist3yr), 
                  data=rbctconf_after[-x,])
  return(predict(m5a_after, rbctconf_after[x,], type = "response"))
})
sq_errors_5a_after <- rep(0, length(LOOCV_5a_after))
abs_errors_5a_after <- rep(0, length(LOOCV_5a_after))
LOOCV_5a_after[[i]]
for(i in 1:length(LOOCV_5a_after)){
  sq_errors_5a_after[i] <- (LOOCV_5a_after[[i]]*rbctconf_after$hdyrsrisk[i]-rbctconf_after$Incidence[i])^2
  abs_errors_5a_after[i] <- abs(LOOCV_5a_after[[i]]*rbctconf_after$hdyrsrisk[i]-rbctconf_after$Incidence[i])
}
cv_mods_after[which(mod == "rs5a_after"), RMSE:= sqrt(mean(sq_errors_5a_after))]
cv_mods_after[which(mod == "rs5a_after"), MAE:= mean(abs_errors_5a_after)]



BIC_rs5a <- BIC(rs5a)
BIC_rs5a
AICc_rs5a <- AICc(rs5a)
AICc_rs5a

BIC_rs5a_follow <- BIC(rs5a_follow)
BIC_rs5a_follow
AICc_rs5a_follow <- AICc(rs5a_follow)
AICc_rs5a_follow

BIC_rs5a_after <- BIC(rs5a_after)
BIC_rs5a_after
AICc_rs5a_after <- AICc(rs5a_after)
AICc_rs5a_after



#AKAIKE----
#Initial
cv_mods[which(mod == "rs1"), Akaike:= AICc_rs1]
cv_mods[which(mod == "rs2"), Akaike:= "N/A"]
cv_mods[which(mod == "rs3"), Akaike:= AICc_rs3]
cv_mods[which(mod == "rs3a"), Akaike:= AICc_rs3a]
cv_mods[which(mod == "rs4"), Akaike:= AICc_rs4]
cv_mods[which(mod == "rs4a"), Akaike:= AICc_rs4a]
cv_mods[which(mod == "rs4c"), Akaike:= AICc_rs4c]
cv_mods[which(mod == "rs4d"), Akaike:= AICc_rs4d]
cv_mods[which(mod == "rs5"), Akaike:= AICc_rs5]
cv_mods[which(mod == "rs5a"), Akaike:= AICc_rs5a]

#From follow-up
cv_mods_follow[which(mod == "rs1_follow"), Akaike:= AICc_rs1_follow]
cv_mods_follow[which(mod == "rs2_follow"), Akaike:= "N/A"]
cv_mods_follow[which(mod == "rs3_follow"), Akaike:= AICc_rs3_follow]
cv_mods_follow[which(mod == "rs3a_follow"), Akaike:= AICc_rs3a_follow]
cv_mods_follow[which(mod == "rs4_follow"), Akaike:= AICc_rs4_follow]
cv_mods_follow[which(mod == "rs4a_follow"), Akaike:= AICc_rs4a_follow]
cv_mods_follow[which(mod == "rs4c_follow"), Akaike:= AICc_rs4c_follow]
cv_mods_follow[which(mod == "rs4d_follow"), Akaike:= AICc_rs4d_follow]
cv_mods_follow[which(mod == "rs5_follow"), Akaike:= AICc_rs5_follow]
cv_mods_follow[which(mod == "rs5a_follow"), Akaike:= AICc_rs5a_follow]
cv_mods_follow

#Post-trial
cv_mods_after[which(mod == "rs1_after"), Akaike:= AICc_rs1_after]
cv_mods_after[which(mod == "rs3_after"), Akaike:= AICc_rs3_after]
cv_mods_after[which(mod == "rs3a_after"), Akaike:= AICc_rs3a_after]
cv_mods_after[which(mod == "rs4_after"), Akaike:= AICc_rs4_after]
cv_mods_after[which(mod == "rs4a_after"), Akaike:= AICc_rs4a_after]
cv_mods_after[which(mod == "rs4c_after"), Akaike:= AICc_rs4c_after]
cv_mods_after[which(mod == "rs4d_after"), Akaike:= AICc_rs4d_after]
cv_mods_after[which(mod == "rs5_after"), Akaike:= AICc_rs5_after]
cv_mods_after[which(mod == "rs5a_after"), Akaike:= AICc_rs5a_after]



#BIC----
cv_mods[which(mod == "rs1"), BIC_val:= BIC_rs1]
cv_mods[which(mod == "rs2"), BIC_val:= "N/A"]
cv_mods[which(mod == "rs3"), BIC_val:= BIC_rs3]
cv_mods[which(mod == "rs3a"), BIC_val:= BIC_rs3a]
cv_mods[which(mod == "rs4"), BIC_val:= BIC_rs4]
cv_mods[which(mod == "rs4a"), BIC_val:= BIC_rs4a]
cv_mods[which(mod == "rs4c"), BIC_val:= BIC_rs4c]
cv_mods[which(mod == "rs4d"), BIC_val:= BIC_rs4d]
cv_mods[which(mod == "rs5"), BIC_val:= BIC_rs5]
cv_mods[which(mod == "rs5a"), BIC_val:= BIC_rs5a]

cv_mods_follow[which(mod == "rs1_follow"), BIC_val:= BIC_rs1_follow]
cv_mods_follow[which(mod == "rs2_follow"), BIC_val:= NA]
cv_mods_follow[which(mod == "rs3_follow"), BIC_val:= BIC_rs3_follow]
cv_mods_follow[which(mod == "rs3a_follow"), BIC_val:= BIC_rs3a_follow]
cv_mods_follow[which(mod == "rs4_follow"), BIC_val:= BIC_rs4_follow]
cv_mods_follow[which(mod == "rs4a_follow"), BIC_val:= BIC_rs4a_follow]
cv_mods_follow[which(mod == "rs5_follow"), BIC_val:= BIC_rs5_follow]
cv_mods_follow[which(mod == "rs5a_follow"), BIC_val:= BIC_rs5a_follow]

cv_mods_after[which(mod == "rs1_after"), BIC_val:= BIC_rs1_after]
cv_mods_after[which(mod == "rs2_after"), BIC_val:= NA]
cv_mods_after[which(mod == "rs3_after"), BIC_val:= BIC_rs3_after]
cv_mods_after[which(mod == "rs3a_after"), BIC_val:= BIC_rs3a_after]
cv_mods_after[which(mod == "rs4_after"), BIC_val:= BIC_rs4_after]
cv_mods_after[which(mod == "rs4a_after"), BIC_val:= BIC_rs4a_after]
cv_mods_after[which(mod == "rs4c_after"), BIC_val:= BIC_rs4c_after]
cv_mods_after[which(mod == "rs4d_after"), BIC_val:= BIC_rs4d_after]
cv_mods_after[which(mod == "rs5_after"), BIC_val:= BIC_rs5_after]
cv_mods_after[which(mod == "rs5a_after"), BIC_val:= BIC_rs5a_after]





#Bayesian Analysis ----
#Bayesian Poisson GLM
set.seed(150499)
rs_pois <- rstanarm::stan_glm(Incidence~Treatment+A+B+C+D+E+F+G+H+I+log(Hist3yr),
                    data = rbctconf, family = "poisson",
                    offset = log(Baseline),
                    prior = normal(0, 1),
                    prior_intercept = normal(0, 2),
                    diagnostic_file = file.path(tempdir(), "df.csv"),refresh=0)
loo_rs_pois <- loo(rs_pois, k_threshold = 0.7)
pp_check(rs_pois, plotfun = "dens_overlay")
draws_rs_pois <- as.data.table(rs_pois)
length(which(draws_rs_pois$TreatmentProactive < 0))/nrow(draws_rs_pois)

loo_rs_pois <- loo(rs_pois, k_threshold = 0.7)


rs_pois_after <- stan_glm(Incidence~Treatment+A+B+C+D+E+F+G+H+I+log(Hist3yr)+log(Baseline),
                          data = rbctconf_after, family = "poisson",
                          offset = log(Baseline),
                          prior = normal(0, 1),
                          prior_intercept = normal(0, 2),
                          diagnostic_file = file.path(tempdir(), "df.csv"),refresh=0)
loo_rs_pois_after <- loo(rs_pois_after, k_threshold = 0.7)
pp_check(rs_pois_after, plotfun = "dens_overlay")
draws_rs_pois_after <- as.data.table(rs_pois_after)






rs_pois_follow <- stan_glm(Incidence~Treatment+A+B+C+D+E+F+G+H+I+log(Hist3yr)+log(Baseline),
                           data = rbctconf_follow, family = "poisson",
                           offset = log(Baseline),
                           prior = normal(0, 1),
                           prior_intercept = normal(0, 2),
                           diagnostic_file = file.path(tempdir(), "df.csv"),refresh=0)
loo_rs_pois_follow <- loo(rs_pois_follow, k_threshold = 0.7)
pp_check(rs_pois_follow, plotfun = "dens_overlay")
draws_rs_pois_follow <- as.data.table(rs_pois_follow)





#Posterior Predictive plotting
rs_pois_posterior_predictive <- posterior_predict(rs_pois)
ppc_stat(rbctconf$Incidence, rs_pois_posterior_predictive, stat = "min")

#Posterior vs Prior
draws_rs_pois <- as.data.table(rs_pois)
rs_pois_po_beta_trmnt <- subset(draws_rs_pois, select = "TreatmentProactive")
setnames(rs_pois_po_beta_trmnt, "TreatmentProactive", "posterior_val")
rs_pois_prior_posterior <- data.table()
rs_pois_prior_posterior[, posterior_val:= rs_pois_po_beta_trmnt$posterior_val]

rs_pois_densities <- ggplot(rs_pois_prior_posterior)+
  stat_function(aes(col = "Prior"), fun = dnorm, n = 101, size = 1.5, args = list(mean = 0, sd = 1))+  
  geom_density(aes(posterior_val, col = "Posterior"), size = 2)+
  theme_bw()+
  labs(x = "Treatment Parameter", y = "Density",
       title = element_blank(),
       subtitle = element_blank())+
  theme(text = element_text(size = 28),
        axis.text.x = element_text(size=28),
        axis.text.y = element_text(size=28),
        plot.title = element_text(size=20, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size=16, hjust = 0.5),
        axis.title=element_text(size=28), 
        legend.text=element_text(size=28)+geom_text(size = 28),
        legend.position = "bottom",
        legend.title = element_blank())+
  scale_color_manual(values = plot_val_cols[c(2,1)])+xlim(c(-2.5, 2.5))
rs_pois_densities

ggsave(rs_pois_densities, file = file.path(out.dir2, "rs_pois_densities.pdf"), h = 16, w = 24)





rs1aB_posterior_predictive <- posterior_predict(rs1aB)
ppc_stat(rbctconf$Incidence, rs1aB_posterior_predictive, stat = "min")

#rs1aB (Improved) Prior vs Posterior----
rs1aB_improved_pois_prior_beta_trmnts_dt <- data.table()
rs1aB_improved_prior_posterior <- copy(rs1aB_improved_pois_prior_beta_trmnts_dt)
draws_rs1aB_improved <- as.data.table(rs1aB_improved)
rs1aB_improved_po_beta_trmnt <- subset(draws_rs1aB_improved, select = "TreatmentProactive")
setnames(rs1aB_improved_po_beta_trmnt, "TreatmentProactive", "posterior_val")
rs1aB_improved_prior_posterior[, posterior_val:= rs1aB_improved_po_beta_trmnt$posterior_val]



rs1aB_improved_densities <- ggplot(rs1aB_improved_prior_posterior)+
  stat_function(aes(col = "Prior"), fun = dnorm, n = 101, size = 1.5, args = list(mean = 0, sd = 1))+  
  geom_density(aes(posterior_val, col = "Posterior"), size = 2)+
  theme_bw()+
  labs(x = "Treatment Parameter", y = "Density",
       title = element_blank(),
       subtitle = element_blank())+
  theme(text = element_text(size = 28),
        axis.text.x = element_text(size=28),
        axis.text.y = element_text(size=28),
        plot.title = element_text(size=20, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size=16, hjust = 0.5),
        axis.title=element_text(size=28), 
        legend.text=element_text(size=28)+geom_text(size = 28),
        legend.position = "bottom",
        legend.title = element_blank())+
  scale_color_manual(values = plot_val_cols[c(2,1)])+xlim(c(-2.5, 2.5))
ggsave(rs1aB_improved_densities, file = file.path(out.dir2, "rs1aB_improved_densities.pdf"), h = 16, w = 24)






update_geom_defaults("line", list(size = 2,
                                  width = 2))

update_geom_defaults("pointrange", list(size = 1,
                                        width = 1,
                                        fatten = 1,
                                        linewidth = 1.5))

rs_post_prior <- posterior_vs_prior(rs,
                                    pars = "TreatmentProactive")+ylim(c(-10, 10))+
  labs(title = "Model a.1: Proactive Treatment Parameter")+theme_bw()+
  theme(text = element_text(size = 26),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=26),
        plot.title = element_text(size=30, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size=30, hjust = 0.5),
        axis.title=element_text(size=28), 
        legend.text=element_text(size=14)+geom_text(size = 20),
        legend.position = "none")+
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank())
rs_post_prior$layers[[1]]$aes_params$colour <- c("#2196F3", "red")
rs_post_prior

ggsave(rs_post_prior, file = file.path(out.dir2, "rs_post_prior.pdf"), h = 20, w = 16)

rs_improved_post_prior <- posterior_vs_prior(rs_improved,
                                             pars = "TreatmentProactive")+ylim(c(-10, 10))+
  labs(title = "Model a.2: Proactive Treatment Parameter")+theme_bw()+
  theme(text = element_text(size = 26),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=26),
        plot.title = element_text(size=30, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size=30, hjust = 0.5),
        axis.title=element_text(size=28), 
        legend.text=element_text(size=14)+geom_text(size = 20),
        legend.position = "none")+
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank())
rs_improved_post_prior$layers[[1]]$aes_params$colour <- c("#2196F3", "red")
rs_improved_post_prior
ggsave(rs_improved_post_prior, file = file.path(out.dir2, "rs_improved_post_prior.pdf"), h = 20, w = 16)

rs_post_prior_grob <- arrangeGrob(rs_post_prior, rs_improved_post_prior, nrow = 1)
ggsave(rs_post_prior_grob, file = file.path(out.dir2, "rs_post_prior_grob.pdf"), h = 12, w = 22)



#rs Dispersion---
rs_post_prior_dispersion <- posterior_vs_prior(rs,
                                               pars = "reciprocal_dispersion")+ylim(c(0, 1100))+
  labs(title = "Model a.1: Reciprocal Dispersion Parameter")+theme_bw()+
  theme(text = element_text(size = 26),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=26),
        plot.title = element_text(size=30, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size=30, hjust = 0.5),
        axis.title=element_text(size=28), 
        legend.text=element_text(size=14)+geom_text(size = 20),
        legend.position = "none")+
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank())
rs_post_prior_dispersion$layers[[1]]$aes_params$colour <- c("#2196F3", "red")
ggsave(rs_post_prior_dispersion, file = file.path(out.dir2, "rs_post_prior_dispersion.pdf"), h = 20, w = 16)


rs_improved_post_prior_dispersion <- posterior_vs_prior(rs_improved,
                                                        pars = "reciprocal_dispersion")+ylim(c(0, 1100))+
  labs(title = "Improved Model a.2: Reciprocal Dispersion Parameter")+theme_bw()+
  theme(text = element_text(size = 26),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=26),
        plot.title = element_text(size=30, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size=30, hjust = 0.5),
        axis.title=element_text(size=28), 
        legend.text=element_text(size=14)+geom_text(size = 20),
        legend.position = "none")+
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank())
rs_improved_post_prior_dispersion$layers[[1]]$aes_params$colour <- c("#2196F3", "red")

rs_post_prior_dispersion_grob <- arrangeGrob(rs_post_prior_dispersion, rs_improved_post_prior_dispersion, nrow = 1)
ggsave(rs_post_prior_dispersion_grob, file = file.path(out.dir2, "rs_post_prior_dispersion_grob.pdf"), h = 14, w = 24)


#Model rs ----
rs <- stan_glm.nb(Incidence~Treatment+A+B+C+D+E+F+G+H+I+log(Hist3yr)+log(Baseline),
                  data = rbctconf, prior_intercept=normal(0, 10),
                  diagnostic_file = file.path(tempdir(), "df.csv"),refresh=0
)
rs$prior.info
rs_posterior_predictive <- posterior_predict(rs)
ppc_stat(rbctconf$Incidence, rs_posterior_predictive, stat = "min")
rs_improved_posterior_predictive <- posterior_predict(rs_improved)
ppc_stat(rbctconf$Incidence, rs_improved_posterior_predictive, stat = "min")


ppc_stat(rbctconf$Incidence, rs_posterior_predictive, stat = "max")
ppc_stat(rbctconf$Incidence, rs_improved_posterior_predictive, stat = "max")


rs_pp <- pp_check(rs, plotfun = "dens_overlay")+   theme_bw()+
  theme(text = element_text(size = 30),axis.text.x = element_text(size=30),
        axis.text.y = element_text(size=30),
        plot.title = element_blank(),
        axis.title=element_text(size=30), 
        legend.text=element_text(size=30),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "bottom")+xlim(c(0,450))+
  labs(title = "", x = "Incidence", y = "Density")+
  guides(color = guide_legend(override.aes = list(linewidth = 5)))
rs_pp
ggsave(rs_pp, file = file.path(out.dir2, "rs_pp.pdf"), h = 12, w = 22)

rs_improved_pp <- pp_check(rs_improved, plotfun = "dens_overlay")+   theme_bw()+
  theme(text = element_text(size = 30),axis.text.x = element_text(size=30),
        axis.text.y = element_text(size=30),
        plot.title = element_blank(),
        axis.title=element_text(size=30), 
        legend.text=element_text(size=30),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "bottom")+xlim(c(0,450))+
  labs(title = "", x = "Incidence", y = "Density")+
  guides(color = guide_legend(override.aes = list(linewidth = 5)))
rs_improved_pp
ggsave(rs_improved_pp, file = file.path(out.dir2, "rs_improved_pp.pdf"), h = 12, w = 22)
rs_posterior_predictive <- posterior_predict(rs)
ppc_stat(rbctconf$Incidence, rs_posterior_predictive, stat = "min")
?pcauchy
integrate(function(x) abs(dcauchy(x, 0, 5)), 0, 1)

#Model rs_improved ----
rs_improved <- stan_glm.nb(Incidence~Treatment+A+B+C+D+E+F+G+H+I+log(Hist3yr)+log(Baseline),
                           data = rbctconf, prior_intercept=normal(0, 1, autoscale=TRUE),
                           prior_aux= cauchy(0, 5),
                           prior = normal(0, 1),
                           diagnostic_file = file.path(tempdir(), "df.csv"),refresh=0)

loo_rs_improved <- loo(rs_improved, k_threshold = 0.7)
loo_rs_improved
pp_check(rs, plotfun = "dens_overlay")
pp_check(rs_improved, plotfun = "dens_overlay")
# posterior_predict(rs_improved)
predictive_interval(rs_improved, prob = 0.95)
predictive_interval(rs, prob = 0.95)

#Small values for reciprocal_dispersion correspond to greater dispersion
draws_rs_improved <- as.data.table(rs_improved)
rs_follow <-stan_glm.nb(Incidence~Treatment+A+B+C+D+E+F+G+H+I+log(Hist3yr)+log(Baseline),prior_intercept=normal(0, 10),
                        prior = normal(0, 10),
                        diagnostic_file = file.path(tempdir(), "df.csv"),refresh=0, data = rbctconf_follow)
pp_check(rs_follow, plotfun = "dens_overlay")
draws_rs_follow <- as.data.table(rs_follow)
loo_rs_follow <- loo(rs_follow)
loo_rs_follow

rs_follow_improved <- stan_glm.nb(Incidence~Treatment+A+B+C+D+E+F+G+H+I+log(Hist3yr)+log(Baseline),
                                  data = rbctconf_follow, prior_intercept=normal(0, 5, autoscale=TRUE),
                                  prior_aux= cauchy(0, 5),
                                  prior = normal(0, 1),
                                  diagnostic_file = file.path(tempdir(), "df.csv"),refresh=0)
loo_rs_follow_improved <- loo(rs_follow_improved, k_threshold = 0.7)
loo_rs_follow_improved

rs_after <-stan_glm.nb(Incidence~Treatment+A+B+C+D+E+F+G+H+I+log(Hist3yr)+log(Baseline),prior_intercept=normal(0, 10),
                       prior = normal(0, 10),
                       diagnostic_file = file.path(tempdir(), "df.csv"),refresh=0, data = rbctconf_after)
pp_check(rs_after, plotfun = "dens_overlay")
loo_rs_after <- loo(rs_after, k_threshold = 0.7)
draws_rs_after <- as.data.table(rs_after)

rs_after_improved <- stan_glm.nb(Incidence~Treatment+A+B+C+D+E+F+G+H+I+log(Hist3yr)+log(Baseline),
                                 data = rbctconf_after, 
                                 prior_aux= cauchy(0, 5),
                                 prior = normal(0, 1),
                                 prior_intercept = normal(0, 1, autoscale = TRUE),
                                 diagnostic_file = file.path(tempdir(), "df.csv"),refresh=0)
loo_rs_after_improved <- loo(rs_after_improved, k_threshold = 0.7)
pp_check(rs_after_improved)
draws_rs_after_improved <- as.data.table(rs_after_improved)

###Dispersion: Prior vs Posterior Density----
rs_prior_dispersion_dt <- data.table()
rs_dispersion_prior_posterior <- copy(rs_prior_dispersion_dt)
draws_rs <- as.data.table(rs)
rs_po_dispersion_trmnt <- subset(draws_rs, select = "reciprocal_dispersion")
setnames(rs_po_dispersion_trmnt, "reciprocal_dispersion", "posterior_val")
rs_dispersion_prior_posterior[, posterior_val:= rs_po_dispersion_trmnt$posterior_val]
rs_dispersion_densities <- ggplot(rs_dispersion_prior_posterior)+
  stat_function(aes(col = "Prior"), fun = dexp, n = 101, size = 1.5, args = list(rate = 1))+
  geom_density(aes(posterior_val, col = "Posterior"))+
  theme_bw()+
  labs(x = "Value", y = "Density",
       title = "Model rs: Posterior vs Prior of Reciprocal Dispersion Parameter",
       subtitle = "(on log scale)")+
  labs(x = "Reciprocal Dispersion Parameter", y = "Density",
       title = element_blank(),
       subtitle = element_blank())+
  theme(text = element_text(size = 28),
        axis.text.x = element_text(size=28),
        axis.text.y = element_text(size=28),
        plot.title = element_text(size=20, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size=16, hjust = 0.5),
        axis.title=element_text(size=28), 
        legend.text=element_text(size=28)+geom_text(size = 28),
        legend.position = "bottom",
        legend.title = element_blank())+
  scale_color_manual(values = plot_val_cols[c(2,1)])+
  xlim(c(0, 600))
ggsave(rs_dispersion_densities, file = file.path(out.dir2, "rs_dispersion_densities.pdf"), h = 16, w = 24)
rs_improved_dispersion_densities

rs_improved_prior_dispersion_dt <- data.table()
rs_improved_dispersion_prior_posterior <- copy(rs_improved_prior_dispersion_dt)
draws_rs_improved <- as.data.table(rs_improved)
rs_improved_po_dispersion_trmnt <- subset(draws_rs_improved, select = "reciprocal_dispersion")
setnames(rs_improved_po_dispersion_trmnt, "reciprocal_dispersion", "posterior_val")
rs_improved_dispersion_prior_posterior[, posterior_val:= rs_improved_po_dispersion_trmnt$posterior_val]

rs_improved_dispersion_densities <- 
  ggplot(rs_improved_dispersion_prior_posterior)+
  stat_function(aes(col = "Prior"), fun = dhcauchy, n = 101, args = list(sigma = 5))+
  geom_density(aes(posterior_val, col = "Posterior"))+
  theme_bw()+
  labs(x = "Value", y = "Density",
       title = "Improved Model rs: Posterior vs Prior of Reciprocal Dispersion Parameter",
       subtitle = "(on log scale)") +
  labs(x = "Reciprocal Dispersion Parameter", y = "Density",
       title = element_blank(),
       subtitle = element_blank())+
  theme(text = element_text(size = 28),
        axis.text.x = element_text(size=28),
        axis.text.y = element_text(size=28),
        plot.title = element_text(size=20, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size=16, hjust = 0.5),
        axis.title=element_text(size=28), 
        legend.text=element_text(size=28)+geom_text(size = 28),
        legend.position = "bottom",
        legend.title = element_blank())+
  scale_color_manual(values = plot_val_cols[c(2,1)])+
  xlim(c(0,600))
rs_improved_dispersion_densities
ggsave(rs_improved_dispersion_densities, file = file.path(out.dir2, "rs_improved_dispersion_densities.pdf"), h = 16, w = 24)
####





rs_improved_pois_prior_beta_trmnts_dt <- data.table()
rs_improved_prior_posterior <- copy(rs_improved_pois_prior_beta_trmnts_dt)
draws_rs_improved <- as.data.table(rs_improved)
rs_improved_po_beta_trmnt <- subset(draws_rs_improved, select = "TreatmentProactive")
setnames(rs_improved_po_beta_trmnt, "TreatmentProactive", "posterior_val")
rs_improved_prior_posterior[, posterior_val:= rs_improved_po_beta_trmnt$posterior_val]
update_geom_defaults("density", list(size = 2))
rs_improved_densities <- ggplot(rs_improved_prior_posterior)+
  stat_function(aes(col = "Prior"), fun = dnorm, n = 101, size = 1.5, args = list(mean = 0, sd = 1))+
  scale_x_continuous(limits = c(-2.5, 2.5))+
  geom_density(aes(posterior_val, col = "Posterior"), size = 2)+
  theme_bw()+
  labs(x = "Treatment Parameter", y = "Density",
       title = element_blank(),
       subtitle = element_blank())+
  theme(text = element_text(size = 28),
        axis.text.x = element_text(size=28),
        axis.text.y = element_text(size=28),
        plot.title = element_text(size=20, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size=16, hjust = 0.5),
        axis.title=element_text(size=28), 
        legend.text=element_text(size=28)+geom_text(size = 28),
        legend.position = "bottom",
        legend.title = element_blank())+
  scale_color_manual(values = plot_val_cols[c(2,1)])
rs_improved_densities
ggsave(rs_improved_densities, file = file.path(out.dir2, "rs_improved_densities.pdf"), h = 16, w = 24)







#Bayesian Model 2----
rs1B<-stan_glm.nb(Incidence~Treatment+log(hdyrsrisk)+log(Hist3yr),
                  prior_intercept=normal(0,10), prior=normal(0,10),
                  data=rbctconf,diagnostic_file =file.path(tempdir(),"df.csv"),refresh=0)

pp_check(rs1B, plotfun = "dens_overlay")



loo_rs1B <- loo(rs1B)
draws_rs1B <- as.data.table(rs1B)

rs1B_improved<-stan_glm.nb(Incidence~Treatment+log(hdyrsrisk)+log(Hist3yr),
                           prior_intercept=normal(0,10, autoscale = TRUE), prior=normal(0,10, autoscale = TRUE),
                           prior_aux = cauchy(0, 5, autoscale = TRUE),
                           data=rbctconf,diagnostic_file =file.path(tempdir(),"df.csv"),refresh=0)
pp_check(rs1B, plotfun = "dens_overlay")
pp_check(rs1B_improved, plotfun = "dens_overlay")
draws_rs1B_improved <- as.data.table(rs1B_improved)

rs1B_follow<-stan_glm.nb(Incidence~Treatment+log(hdyrsrisk)+log(Hist3yr),
                         prior_intercept=normal(0,10), prior=normal(0,10),
                         data=rbctconf_follow,diagnostic_file =file.path(tempdir(),"df.csv"),refresh=0)
loo_rs1B_follow <- loo(rs1B_follow)
pp_check(rs1B_follow, plotfun = "dens_overlay")

draws_rs1B_follow <- as.data.table(rs1B_follow)

rs1B_follow_improved<-stan_glm.nb(Incidence~Treatment+log(hdyrsrisk)+log(Hist3yr),
                                  prior_intercept = normal(0, 5, autoscale = TRUE),
                                  prior = normal(0, 1, autoscale = TRUE),
                                  prior_aux= cauchy(0, 5, autoscale = TRUE),
                                  data=rbctconf_follow,diagnostic_file =file.path(tempdir(),"df.csv"),refresh=0)
draws_rs1B_follow_improved <- as.data.table(rs1B_follow_improved)
pp_check(rs1B_follow, plotfun = "dens_overlay")
pp_check(rs1B_follow_improved, plotfun = "dens_overlay")

loo_rs1B_follow_improved <- loo(rs1B_follow_improved, k_threshold = 0.7)

rs1B_after<-stan_glm.nb(Incidence~Treatment+log(hdyrsrisk)+log(Hist3yr),
                        prior_intercept=normal(0,10), prior=normal(0,10),
                        data=rbctconf_after,diagnostic_file =file.path(tempdir(),"df.csv"),refresh=0)
loo_rs1B_after <- loo(rs1B_after)
pp_check(rs1B_after, plotfun = "dens_overlay")

draws_rs1B_after <- as.data.table(rs1B_after)

rs1B_after_improved<-stan_glm.nb(Incidence~Treatment+log(hdyrsrisk)+log(Hist3yr),
                                 prior_intercept = normal(0, 1, autoscale = TRUE),
                                 prior = normal(0, 1, autoscale = TRUE),
                                 prior_aux= cauchy(0, 5, autoscale = TRUE),
                                 data=rbctconf_after,diagnostic_file =file.path(tempdir(),"df.csv"),refresh=0)
draws_rs1B_after_improved <- as.data.table(rs1B_after_improved)
loo_rs1B_after_improved <- loo(rs1B_after_improved, k_threshold = 0.7)

#Model rs1ab----
#Uses herd years as offset
rs1aB<-stan_glm.nb(Incidence~Treatment+log(hdyrsrisk)+log(Hist3yr),
                   offset = log(hdyrsrisk), 
                   prior_intercept=normal(0,10), prior=normal(0,10),
                   data = rbctconf, refresh=0)
pp_check(rs1aB, plotfun = "dens_overlay")
pp_check(rs1B, plotfun = "dens_overlay")

loo_rs1aB <- loo(rs1aB)
draws_rs1aB <- as.data.table(rs1aB)

rs1aB_improved<-stan_glm.nb(Incidence~Treatment+log(hdyrsrisk)+log(Hist3yr),
                            prior_intercept = normal(0, 1, autoscale = TRUE),
                            prior = normal(0, 1, autoscale = TRUE),
                            prior_aux= cauchy(0, 5, autoscale = TRUE),
                            offset = log(hdyrsrisk),
                            data = rbctconf, refresh=0)



rs1aB_follow<-stan_glm.nb(Incidence~Treatment+log(hdyrsrisk)+log(Hist3yr),
                          offset = log(hdyrsrisk),
                          prior_intercept=normal(0,10), prior=normal(0,10),
                          data = rbctconf_follow, refresh=0)

summary(rs1aB_follow, probs=c(0.025,0.5,0.975), digits=3)
loo_rs1aB_follow <- loo(rs1aB_follow)
loo_rs1aB_follow
draws_rs1aB_follow <- as.data.table(rs1aB_follow)

rs1aB_improved_post_plot <- mcmc_areas(draws_rs1aB_improved, pars = "TreatmentProactive", prob = 0.95) + xlab('Log Effect')+
  ggplot2::theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20))+
  theme(text = element_text(size = 22), 
        plot.title = element_text(hjust = 0.5, face = "bold", size = 24))+
  theme(plot.subtitle = element_text(hjust = 0.5, size = 18),
        axis.text.x = element_text(size = 28),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 28)
  )
rs1aB_improved_post_plot
ggsave(rs1aB_improved_post_plot, file = file.path(out.dir2, "rs1aB_improved_post_plot.pdf"), h = 8, w = 16)

rs1aB_follow_improved<-stan_glm.nb(Incidence~Treatment+log(hdyrsrisk)+log(Hist3yr),
                                   prior_intercept = normal(0, 5, autoscale = TRUE),
                                   prior = normal(0, 1, autoscale = TRUE),
                                   prior_aux= cauchy(0, 5, autoscale = TRUE),
                                   offset = log(hdyrsrisk),
                                   data = rbctconf_follow, refresh=0)
pp_check(rs1aB_follow_improved, plotfun = "dens_overlay")
loo_rs1aB_follow_improved <- loo(rs1aB_follow_improved, k_threshold = 0.7)
draws_rs1aB_follow_improved <- as.data.table(rs1aB_follow_improved)


rs1aB_follow_improved_post_plot <- mcmc_areas(draws_rs1aB_follow_improved, pars = "TreatmentProactive", prob = 0.95) + xlab('Log Effect')+
  ggplot2::theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20))+
  theme(text = element_text(size = 22), 
        plot.title = element_text(hjust = 0.5, face = "bold", size = 24))+
  theme(plot.subtitle = element_text(hjust = 0.5, size = 18),
        axis.text.x = element_text(size = 28),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 28)
  )
rs1aB_follow_improved_post_plot
ggsave(rs1aB_follow_improved_post_plot, file = file.path(out.dir2, "rs1aB_follow_improved_post_plot.pdf"), h = 8, w = 16)







rs1aB_after<-stan_glm.nb(Incidence~Treatment+log(hdyrsrisk)+log(Hist3yr),
                         offset = log(hdyrsrisk),
                         prior_intercept=normal(0,10), prior=normal(0,10),
                         data = rbctconf_after, refresh=0)

loo_rs1aB_after <- loo(rs1aB_after)
pp_check(rs1aB_after)
draws_rs1aB_after <- as.data.table(rs1aB_after)

rs1aB_improved_post_plot <- mcmc_areas(draws_rs1aB_improved, pars = "TreatmentProactive", prob = 0.95) + xlab('Log Effect')+
  ggplot2::theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13))+
  theme(plot.subtitle = element_text(hjust = 0.5, size = 12))+
  labs(title = "Posterior Distribution of Treatment Parameter", subtitle = "with medians and 95% credible intervals")
rs1aB_improved_post_plot


rs1aB_after_improved<-stan_glm.nb(Incidence~Treatment+log(hdyrsrisk)+log(Hist3yr),
                                  prior_intercept = normal(0, 1, autoscale = TRUE),
                                  prior = normal(0, 1, autoscale = TRUE),
                                  prior_aux= cauchy(0, 5, autoscale = TRUE),
                                  offset = log(hdyrsrisk),
                                  data = rbctconf_after, refresh=0)
pp_check(rs1aB_after, plotfun = "dens_overlay")
pp_check(rs1aB_after_improved, plotfun = "dens_overlay")
loo_rs1aB_after_improved <- loo(rs1aB_after_improved, k_threshold = 0.7)
draws_rs1aB_after_improved <- as.data.table(rs1aB_after_improved)


rs1aB_after_improved_post_plot <- mcmc_areas(draws_rs1aB_after_improved, pars = "TreatmentProactive", prob = 0.95) + xlab('Log Effect')+
  ggplot2::theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))+
  theme(plot.subtitle = element_text(hjust = 0.5, size = 15),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        text = element_text(size = 15))+
  labs(title = "Posterior Distribution of Treatment Parameter", subtitle = "with medians and 95% credible intervals")
rs1aB_after_improved_post_plot
ggsave(rs1aB_after_improved_post_plot, file = file.path(out.dir2, "rs1aB_after_improved_post_plot.pdf"), h = 8, w = 16)





rs2B<-stan_glm.nb(Incidence~log(hdyrsrisk)+log(Hist3yr),
                  prior_intercept=normal(0,10), prior=normal(0,10),
                  data=rbctconf,diagnostic_file = file.path(tempdir(), "df.csv"),refresh=0)
loo_rs2B <- loo(rs2B, k_threshold = 0.7)

pp_check(rs2B, plotfun = "dens_overlay")
pp_check(rs1B, plotfun = "dens_overlay")

rs2B_improved<-stan_glm.nb(Incidence~log(hdyrsrisk)+log(Hist3yr),
                           prior_intercept=normal(0,1, autoscale = TRUE), prior=normal(0,1, autoscale = TRUE),
                           prior_aux = cauchy(0, 5, autoscale = TRUE),
                           data=rbctconf,diagnostic_file = file.path(tempdir(), "df.csv"),refresh=0)
loo_rs2B_improved <- loo(rs2B_improved, k_threshold = 0.7)

pp_check(rs2B_improved, plotfun = "dens_overlay")
pp_check(rs1B, plotfun = "dens_overlay")


rs2B_follow<-stan_glm.nb(Incidence~log(hdyrsrisk)+log(Hist3yr),
                         prior_intercept=normal(0,10), prior=normal(0,10),
                         data=rbctconf_follow,diagnostic_file = file.path(tempdir(), "df.csv"),refresh=0)
loo_rs2B_follow <- loo(rs2B_follow)
pp_check(rs2B_follow, plotfun = "dens_overlay")



rs2B_follow_improved<-stan_glm.nb(Incidence~log(hdyrsrisk)+log(Hist3yr),
                                  prior_intercept=normal(0,1, autoscale = TRUE), prior=normal(0,1, autoscale = TRUE),
                                  prior_aux = cauchy(0, 5, autoscale = TRUE),
                                  data=rbctconf_follow,diagnostic_file = file.path(tempdir(), "df.csv"),refresh=0)
loo_rs2B_follow_improved <- loo(rs2B_follow_improved, k_threshold = 0.7)
loo_rs2B_follow_improved
pp_check(rs2B_follow_improved, plotfun = "dens_overlay")

rs2B_after<-stan_glm.nb(Incidence~log(hdyrsrisk)+log(Hist3yr),
                        prior_intercept=normal(0,10), prior=normal(0,10),
                        data=rbctconf_after,diagnostic_file = file.path(tempdir(), "df.csv"),refresh=0)
loo_rs2B_after <- loo(rs2B_after)
loo_rs2B_after
pp_check(rs2B_after, plotfun = "dens_overlay")



rs2B_after_improved<-stan_glm.nb(Incidence~log(hdyrsrisk)+log(Hist3yr),
                                 prior_intercept=normal(0,10, autoscale = TRUE),
                                 prior_aux = cauchy(0, 5, autoscale = TRUE),
                                 prior = normal(0, 10),
                                 data=rbctconf_after,diagnostic_file = file.path(tempdir(), "df.csv"),refresh=0)
loo_rs2B_after_improved <- loo(rs2B_after_improved, k_threshold = 0.7)
loo_rs2B_after_improved
pp_check(rs2B_after_improved, plotfun = "dens_overlay")






#Posterior Median Estimates of exponentiated (log-scale) treatment parameter -----
quantile(exp(draws_rs_pois$TreatmentProactive), probs = c(0.025, 0.5, 0.975))-1
quantile(exp(draws_rs_pois_after$TreatmentProactive), probs = c(0.025, 0.5, 0.975))-1

quantile(exp(draws_rs$TreatmentProactive), probs = c(0.025, 0.5, 0.975))-1
quantile(exp(draws_rs_follow_improved$TreatmentProactive), probs = c(0.025, 0.5, 0.975))-1
quantile(exp(draws_rs_after$TreatmentProactive), probs = c(0.025, 0.5, 0.975))-1

quantile(exp(draws_rs_improved$TreatmentProactive), probs = c(0.025, 0.5, 0.975))-1
quantile(exp(draws_rs_follow$TreatmentProactive), probs = c(0.025, 0.5, 0.975))-1
quantile(exp(draws_rs_after_improved$TreatmentProactive), probs = c(0.025, 0.5, 0.975))-1


quantile(exp(draws_rs1B$TreatmentProactive), probs = c(0.025, 0.5, 0.975))-1
quantile(exp(draws_rs1B_follow$TreatmentProactive), probs = c(0.025, 0.5, 0.975))-1
quantile(exp(draws_rs1B_after$TreatmentProactive), probs = c(0.025, 0.5, 0.975))-1

quantile(exp(draws_rs1B_improved$TreatmentProactive), probs = c(0.025, 0.5, 0.975))-1
quantile(exp(draws_rs1B_follow_improved$TreatmentProactive), probs = c(0.025, 0.5, 0.975))-1
quantile(exp(draws_rs1B_after_improved$TreatmentProactive), probs = c(0.025, 0.5, 0.975))-1


quantile(exp(draws_rs1aB$TreatmentProactive), probs = c(0.025, 0.5, 0.975))-1
quantile(exp(draws_rs1aB_follow$TreatmentProactive), probs = c(0.025, 0.5, 0.975))-1
quantile(exp(draws_rs1aB_after$TreatmentProactive), probs = c(0.025, 0.5, 0.975))-1

quantile(exp(draws_rs1aB_improved$TreatmentProactive), probs = c(0.025, 0.5, 0.975))-1
quantile(exp(draws_rs1aB_follow_improved$TreatmentProactive), probs = c(0.025, 0.5, 0.975))-1
quantile(exp(draws_rs1aB_after_improved$TreatmentProactive), probs = c(0.025, 0.5, 0.975))-1


#Bayesian Plotting -----
rs1aB_improved_post_plot <- mcmc_areas(draws_rs1aB_improved, pars = "TreatmentProactive", prob = 0.95) + xlab('Log Effect')+
  ggplot2::theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20))+
  theme(plot.subtitle = element_text(hjust = 0.5, size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_blank(),
        text = element_text(size = 18))+
  labs(title = "Improved Model rs1aB: Posterior Distribution of Treatment Parameter", subtitle = "with medians and 95% credible intervals")
rs1aB_improved_post_plot
ggsave(rs1aB_improved_post_plot, file = file.path(out.dir2, "rs_1aBimproved_follow_post_plot.pdf"), h = 8, w = 16)





#Model rs1 and Poisson Effects Together ----
tmp <- copy(draws_rs_improved)
dens_nb <- density(tmp$TreatmentProactive)
nb_df <- data.frame(x=dens_nb$x, y=dens_nb$y)
probs <- c(0.025, 0.975)
quantiles_nb <- quantile(tmp$TreatmentProactive, probs)
dens_nb$quant <- factor(findInterval(nb_df$x,quantiles_nb))
df <- data.table(x = dens_nb$x,
                 y = dens_nb$y,
                 quant = dens_nb$quant)

tmp <- copy(draws_rs_pois)
dens_pois <- density(tmp$TreatmentProactive)
pois_df <- data.frame(x=dens_pois$x, y=dens_pois$y)
probs <- c(0.025, 0.975)
quantiles_pois <- quantile(tmp$TreatmentProactive, probs)
quantiles_pois
dens_pois$quant <- factor(findInterval(pois_df$x,quantiles_pois), 
                          labels = c(4, 5, 6))

df_pois <- data.table(x = dens_pois$x,
                      y = dens_pois$y,
                      quant = dens_pois$quant)
df <- rbind(df, df_pois)
df[, variable:= c(rep("nb", nrow(df)/2),
                  rep("pois", nrow(df)/2))]
pois_1ab_initial_plot <- ggplot(
  data = df
)+
  geom_line(
    aes(
      x=x,
      y=y,
      color = variable
    )
  )+  geom_ribbon(
    aes(
      x=x, ymin=0, ymax=y, fill = quant,
      color = variable
    ),
    alpha = 0.25, show.legend=FALSE
  )+geom_vline(aes(xintercept = 0), linetype = "dashed",
               linewidth = 0.4)+
  scale_fill_manual(values = c(nb_quant_cols,
                               pois_quant_cols),
                    guide = "none")+theme_bw()+
  scale_color_manual(values = c(nb_quant_cols[2],pois_quant_cols[2]),
                     name = "Model", labels = c("Negative Binomial",
                                                "Poisson"))+
  theme(legend.position = "bottom")+
  scale_x_continuous(breaks = seq(-0.8, 0.8, by = 0.1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20))+
  theme(text = element_text(size = 22), 
        plot.title = element_text(hjust = 0.5, face = "bold", size = 24))+
  theme(plot.subtitle = element_text(hjust = 0.5, size = 18),
        axis.text.x = element_text(size = 28),  
        legend.title=element_text(size=30), 
        legend.text=element_text(size=30),
        axis.text.y = element_text(size = 28),
        axis.title = element_text(size = 28))+
  guides(
    fill = guide_legend(show = FALSE) 
  )+ xlab('Log Effect')+ylab("Density")
pois_1ab_initial_plot <- pois_1ab_initial_plot+ guides(color = guide_legend(override.aes = list(linewidth = 4)))
pois_1ab_initial_plot
ggsave(pois_1ab_initial_plot, file = file.path(out.dir2, "pois_1ab_initial_plot.pdf"), h = 12, w = 24)  
 


#Model 1aB and Poisson Effects Together (follow-up) ----
tmp <- copy(draws_rs1aB_follow_improved)
dens_nb <- density(tmp$TreatmentProactive)
nb_df <- data.frame(x=dens_nb$x, y=dens_nb$y)
probs <- c(0.025, 0.975)
quantiles_nb <- quantile(tmp$TreatmentProactive, probs)
quantiles_nb
factor(findInterval(nb_df$x,quantiles_nb))
dens_nb$quant <- factor(findInterval(nb_df$x,quantiles_nb))
df <- data.table(x = dens_nb$x,
                 y = dens_nb$y,
                 quant = dens_nb$quant)

tmp <- copy(draws_rs_pois_follow)
dens_pois <- density(tmp$TreatmentProactive)
pois_df <- data.frame(x=dens_pois$x, y=dens_pois$y)
probs <- c(0.025, 0.975)
quantiles_pois <- quantile(tmp$TreatmentProactive, probs)
quantiles_pois
factor(findInterval(pois_df$x,quantiles_pois))

dens_pois$quant <- factor(findInterval(pois_df$x,quantiles_pois), 
                          labels = c(4, 5, 6))
df_pois <- data.table(x = dens_pois$x,
                      y = dens_pois$y,
                      quant = dens_pois$quant)

df <- rbind(df, df_pois)
df[, variable:= c(rep("nb", nrow(df)/2),
                  rep("pois", nrow(df)/2))]
pois_1ab_follow_plot <- ggplot(
  data = df
)+
  geom_line(
    aes(
      x=x,
      y=y,
      color = variable
    )
  )+  geom_ribbon(
    aes(
      x=x, ymin=0, ymax=y, fill = quant,
      color = variable
    ),
    alpha = 0.25, show.legend=FALSE
  )+geom_vline(aes(xintercept = 0), linetype = "dashed",
               linewidth = 0.4)+
  scale_fill_manual(values = c(nb_quant_cols,
                               pois_quant_cols),
                    guide = "none")+theme_bw()+
  scale_color_manual(values = c(nb_quant_cols[2],pois_quant_cols[2]),
                     name = "Model", labels = c("Negative Binomial",
                                                "Poisson"))+
  theme(legend.position = "bottom")+
  scale_x_continuous(breaks = seq(-1.4, 1.4, by = 0.2),
                     labels = scales::number_format(accuracy = 0.01))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20))+
  theme(text = element_text(size = 22), 
        plot.title = element_text(hjust = 0.5, face = "bold", size = 24))+
  theme(plot.subtitle = element_text(hjust = 0.5, size = 18),
        axis.text.x = element_text(size = 28),  
        legend.title=element_text(size=30), 
        legend.text=element_text(size=30),
        axis.text.y = element_text(size = 28),
        axis.title = element_text(size = 28))+
  guides(
    fill = guide_legend(show = FALSE) 
  )+ xlab('Log Effect')+ylab("Density")
pois_1ab_follow_plot <- pois_1ab_follow_plot+ guides(color = guide_legend(override.aes = list(linewidth = 4)))

ggsave(pois_1ab_follow_plot, file = file.path(out.dir2, "pois_1ab_follow_plot.pdf"), h = 12, w = 24) 



#Model 1aB and Poisson Effects Together (post-trial) ----
tmp <- copy(draws_rs1aB_after_improved)
dens_nb <- density(tmp$TreatmentProactive)
nb_df <- data.frame(x=dens_nb$x, y=dens_nb$y)
probs <- c(0.025, 0.975)
quantiles_nb <- quantile(tmp$TreatmentProactive, probs)
quantiles_nb
factor(findInterval(nb_df$x,quantiles_nb))
dens_nb$quant <- factor(findInterval(nb_df$x,quantiles_nb))
df <- data.table(x = dens_nb$x,
                 y = dens_nb$y,
                 quant = dens_nb$quant)

tmp <- copy(draws_rs_pois_after)
dens_pois <- density(tmp$TreatmentProactive)
pois_df <- data.frame(x=dens_pois$x, y=dens_pois$y)
probs <- c(0.025, 0.975)
quantiles_pois <- quantile(tmp$TreatmentProactive, probs)
quantiles_pois
factor(findInterval(pois_df$x,quantiles_pois))

dens_pois$quant <- factor(findInterval(pois_df$x,quantiles_pois), 
                          labels = c(4, 5, 6))
df_pois <- data.table(x = dens_pois$x,
                      y = dens_pois$y,
                      quant = dens_pois$quant)

df <- rbind(df, df_pois)
df[, variable:= c(rep("nb", nrow(df)/2),
                  rep("pois", nrow(df)/2))]
pois_rs_after_plot <- ggplot(
  data = df
)+
  geom_line(
    aes(
      x=x,
      y=y,
      color = variable
    )
  )+  geom_ribbon(
    aes(
      x=x, ymin=0, ymax=y, fill = quant,
      color = variable
    ),
    alpha = 0.25, show.legend=FALSE
  )+geom_vline(aes(xintercept = 0), linetype = "dashed",
               linewidth = 0.4)+
  scale_fill_manual(values = c(nb_quant_cols,
                               pois_quant_cols),
                    guide = "none")+theme_bw()+
  scale_color_manual(values = c(nb_quant_cols[2],pois_quant_cols[2]),
                     name = "Model", labels = c("Negative Binomial",
                                                "Poisson"))+
  theme(legend.position = "bottom")+
  scale_x_continuous(breaks = seq(-1.4, 1.4, by = 0.2),
                     labels = scales::number_format(accuracy = 0.01))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20))+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20))+
  theme(text = element_text(size = 22), 
        plot.title = element_text(hjust = 0.5, face = "bold", size = 24))+
  theme(plot.subtitle = element_text(hjust = 0.5, size = 18),
        axis.text.x = element_text(size = 28),  
        legend.title=element_text(size=30), 
        legend.text=element_text(size=30),
        axis.text.y = element_text(size = 28),
        axis.title = element_text(size = 28))+
  guides(
    fill = guide_legend(show = FALSE) 
  )+ xlab('Log Effect')+ylab("Density")
pois_rs_after_plot <- pois_rs_after_plot + guides(color = guide_legend(override.aes = list(linewidth = 4)))
ggsave(pois_rs_after_plot, file = file.path(out.dir2, "pois_rs_after_plot.pdf"), h = 12, w = 24) 


rs2B_pp <- pp_check(rs2B, plotfun = "dens_overlay")+ 
  theme(text = element_text(size = 20),axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        plot.title = element_text(size=22, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size=20, hjust = 0.5),
        axis.title=element_text(size=20), 
        legend.text=element_text(size=18))+
  theme(legend.position = "bottom")+xlim(c(0,260))+
  labs(title = "Posterior Predictive Check of Model rs2B", x = "Incidence", y = "Density")
rs2B_pp

rs1aB_pp <- pp_check(rs1aB, plotfun = "dens_overlay")+ 
  theme(text = element_text(size = 20),axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        plot.title = element_text(size=22, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size=20, hjust = 0.5),
        axis.title=element_text(size=20), 
        legend.text=element_text(size=18))+
  theme(legend.position = "bottom")+
  xlim(c(0,260))+
  labs(title = "Posterior Predictive Check of Model rs1aB", x = "Incidence", y = "Density")
rs1aB_pp
ggsave(rs1aB_pp, file = file.path(out.dir2, "rs1aB_pp.pdf"), h = 14, w = 8)

rs1aB_improved_pp <- pp_check(rs1aB_improved, plotfun = "dens_overlay")+ 
  theme(text = element_text(size = 20),axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        plot.title = element_text(size=22, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size=20, hjust = 0.5),
        axis.title=element_text(size=20), 
        legend.text=element_text(size=18))+
  theme(legend.position = "bottom")+xlim(c(0,260))+
  labs(title = "Posterior Predictive Check of Improved Model rs1aB", x = "Incidence", y = "Density")

ggsave(rs1aB_improved_pp, file = file.path(out.dir2, "rs1aB_improved_pp.pdf"), h = 14, w = 8)

rs1aB_improved_posterior_predictive <- posterior_predict(rs1aB_improved)
ppc_max_rs1aB <- ppc_stat(rbctconf$Incidence, rs1aB_posterior_predictive, stat = "max")+
  xlim(c(0, 400))+ 
  theme(text = element_text(size = 28),axis.text.x = element_text(size=28),
        axis.text.y = element_text(size=28),
        plot.title = element_text(size=28, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size=20, hjust = 0.5),
        axis.title=element_text(size=28), 
        legend.text=element_text(size=28))+
  theme(legend.position = "bottom")+
  labs(title = "Model rs1aB Posterior Predictive Check: Test Statistic (Maximum)", x = "Maximum Incidence", y = "Count")

ppc_max_rs1aB_improved <-ppc_stat(rbctconf$Incidence, rs1aB_improved_posterior_predictive, stat = "max")+
  xlim(c(0, 400))+ 
  theme(text = element_text(size = 28),axis.text.x = element_text(size=28),
        axis.text.y = element_text(size=28),
        plot.title = element_text(size=28, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size=20, hjust = 0.5),
        axis.title=element_text(size=28), 
        legend.text=element_text(size=28))+
  theme(legend.position = "bottom")+
  labs(title = "Improved Model rs1aB Posterior Predictive Check: Test Statistic (Maximum)", x = "Maximum Incidence", y = "Count")

rs_pois_posterior_predictive <- posterior_predict(rs_pois)
ppc_max_rs_pois <-ppc_stat(rbctconf$Incidence, rs_pois_posterior_predictive, stat = "max")+
  xlim(c(0, 400))+ 
  theme(text = element_text(size = 28),axis.text.x = element_text(size=28),
        axis.text.y = element_text(size=28),
        plot.title = element_text(size=28, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size=20, hjust = 0.5),
        axis.title=element_text(size=28), 
        legend.text=element_text(size=28))+
  theme(legend.position = "bottom")+
  labs(title = "Poisson GLM Posterior Predictive Check: Test Statistic (Maximum)", x = "Maximum Incidence", y = "Count")
ppc_max_rs_pois
grid.arrange(ppc_max_rs1aB, ppc_max_rs1aB_improved, ppc_max_rs_pois, nrow = 3)
col_max_ppcs_initial <- arrangeGrob(ppc_max_rs1aB, ppc_max_rs1aB_improved, ppc_max_rs_pois, nrow = 3)
ggsave(col_max_ppcs_initial, file = file.path(out.dir2, "col_max_ppcs_initial.pdf"), h = 24, w = 17)
dim(rs_pois_posterior_predictive)
rs_pois_pp <- pp_check(rs_pois, plotfun = "dens_overlay")+ 
  theme(text = element_text(size = 20),axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        plot.title = element_text(size=22, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size=20, hjust = 0.5),
        axis.title=element_text(size=20), 
        legend.text=element_text(size=18))+xlim(c(0,260))+
  theme(legend.position = "bottom")+
  labs(title = "Posterior Predictive Check of Poisson Model", x = "Incidence", y = "Density")
rs_pois_pp
ggsave(rs_pois_pp, file = file.path(out.dir2, "rs_pois_pp.pdf"), h = 14, w = 8)

grid.arrange(rs1aB_pp, rs1aB_improved_pp, rs2B_pp, rs_pois_pp, nrow = 2)
grid_pps_initial_period <- arrangeGrob(rs1aB_pp, rs1aB_improved_pp, rs2B_pp, rs_pois_pp, nrow = 2)
ggsave(grid_pps_initial_period, file = file.path(out.dir2, "grid_pps_initial_period.pdf"), h = 14, w = 17)



#TOTAL BREAKDOWNS (FREQUENTIST ANALYSIS) -----

#Model 7----
rbct_allbreakdowns$Triplet<-relevel(factor(rbct_allbreakdowns$Triplet), ref="J")
rbct_allbreakdowns$Treatment<-relevel(factor(rbct_allbreakdowns$Treatment), ref="Survey-only")

rbct_allbreakdowns_follow$Triplet<-relevel(factor(rbct_allbreakdowns_follow$Triplet), ref="J")
rbct_allbreakdowns_follow$Treatment<-relevel(factor(rbct_allbreakdowns_follow$Treatment), ref="Survey-only")


rbct_allbreakdowns_after$Triplet<-relevel(factor(rbct_allbreakdowns_after$Triplet), ref="J")
rbct_allbreakdowns_after$Treatment<-relevel(factor(rbct_allbreakdowns_after$Treatment), ref="Survey-only")

rs7<-glm(round(Incidence/0.85)~Treatment+Triplet+log(Hist3yr)+log(Baseline),quasipoisson, data=rbct_allbreakdowns)
summary(rs7)

rs7_se <- sqrt(diag(vcov(rs7)))
exp(rs7$coefficients[2])-1

exp(rs7$coefficients[2]-qnorm(0.975) * rs7_se[2])-1
exp(rs7$coefficients[2]+qnorm(0.975) * rs7_se[2])-1

rs7_follow<-glm(Incidence~Treatment+Triplet+log(Hist3yr)+log(Baseline),quasipoisson, data=rbct_allbreakdowns_follow)
summary(rs7_follow)
exp(rs7_follow$coefficients[2])-1

rs7_follow_se <- sqrt(diag(vcov(rs7_follow)))
exp(rs7_follow$coefficients[2]-qnorm(0.975) * rs7_follow_se[2])-1
exp(rs7_follow$coefficients[2]+qnorm(0.975) * rs7_follow_se[2])-1


rs7_after <-glm(Incidence~Treatment+Triplet+log(Hist3yr)+log(Baseline),
                quasipoisson, data=rbct_allbreakdowns_after)
summary(rs7_after)
exp(rs7_after$coefficients[2])-1

rs7_after_se <- sqrt(diag(vcov(rs7_after)))
exp(rs7_after$coefficients[2]-qnorm(0.975) * rs7_after_se[2])-1
exp(rs7_after$coefficients[2]+qnorm(0.975) * rs7_after_se[2])-1

colSums( influence.measures(rs7)$is.inf )
influence.measures(rs7_after)$is.inf
colSums( influence.measures(rs7_follow)$is.inf )
colSums( influence.measures(rs7_after)$is.inf )
rbctconf_after[c(13,14)]


#Model 7a
rs7a<-glmmTMB(round(Incidence/0.85)~Treatment+Triplet+log(Hist3yr)+log(hrdyrsrisk),family=genpois, data= rbct_allbreakdowns)
summary(rs7a)
AICc(rs7)
BIC(rs7)
AICc(rs7a)
BIC(rs7a)
AICc(rs7b)
BIC(rs7b)

check_model(rs7a)
plot(simulateResiduals(rs7a))
summary(rs7a$sdr, "fixed")
rs7a_se <- summary(rs7a$sdr, "fixed")[2,2]
exp(summary(rs7a$sdr, "fixed")[2,1])-1
exp(summary(rs7a$sdr, "fixed")[2,1]+qnorm(0.975) * rs7a_se)-1
exp(summary(rs7a$sdr, "fixed")[2,1]-qnorm(0.975) * rs7a_se)-1

rs7a_follow2<-glmmTMB((Incidence)~Treatment+Triplet+log(Hist3yr)+log(hdyrsrisk),family=genpois, data= rbct_allbreakdowns_follow)
summary(rs7a_follow2) #Without division by 0.85
rs7a_follow<-glmmTMB(round(Incidence/0.85)~Treatment+Triplet+log(Hist3yr)+log(hdyrsrisk),family=genpois, data= rbct_allbreakdowns_follow)
summary(rs7a_follow)
check_model(rs7a_follow)
AICc(rs7a_follow)
plot(simulateResiduals(rs7a_follow), quantreg = TRUE)
summary(rs7a_follow$sdr, "fixed")[2,1]
rs7a_follow_se <- summary(rs7a_follow$sdr, "fixed")[2,2]
exp(summary(rs7a_follow$sdr, "fixed")[2,1])-1

exp(summary(rs7a_follow$sdr, "fixed")[2,1]+qnorm(0.975) * rs7a_follow_se)-1
exp(summary(rs7a_follow$sdr, "fixed")[2,1]-qnorm(0.975) * rs7a_follow_se)-1


mod_vector_all <- c("rs7", "rs7a", "rs7b", "rs8")
all_cv_mods <- data.table(mod = mod_vector_all, 
                                         RMSE = rep(0, length(mod_vector_all)),
                                         MAE = rep(0, length(mod_vector_all)), 
                                         Rsquared = rep(0, length(mod_vector_all)))

mod_vector_all_follow <- c("rs7_follow", "rs7a_follow", "rs7b_follow",
                           "rs8_follow")
all_cv_mods_follow <- data.table(mod = mod_vector_all_follow, 
                                                RMSE = rep(0, length(mod_vector_all_follow)),
                                                MAE = rep(0, length(mod_vector_all_follow)), 
                                                Rsquared = rep(0, length(mod_vector_all_follow)))

mod_vector_all_after <- c("rs7_after", "rs7a_after", "rs7b_after",
                          "rs8_after")
all_cv_mods_after <- data.table(mod = mod_vector_all_after, 
                                               RMSE = rep(0, length(mod_vector_all_after)),
                                               MAE = rep(0, length(mod_vector_all_after)), 
                                               Rsquared = rep(0, length(mod_vector_all_after)))

rs7_cv <- train(round(Incidence/0.85)~Treatment+Triplet+log(Hist3yr)+log(Baseline), family = quasipoisson, data=rbct_allbreakdowns, method = "glm",
                trControl = train.control)
all_cv_mods[which(mod == "rs7"), RMSE:= rs7_cv$results$RMSE]
all_cv_mods[which(mod == "rs7"), Rsquared:= rs7_cv$results$Rsquared]
all_cv_mods[which(mod == "rs7"), MAE:= rs7_cv$results$MAE]

rs7_follow_cv <- train(Incidence~Treatment+Triplet+log(Hist3yr)+log(Baseline), family = quasipoisson, data=rbct_allbreakdowns_follow, method = "glm",
                       trControl = train.control)
all_cv_mods_follow[which(mod == "rs7_follow"), RMSE:= rs7_follow_cv$results$RMSE]
all_cv_mods_follow[which(mod == "rs7_follow"), Rsquared:= rs7_follow_cv$results$Rsquared]
all_cv_mods_follow[which(mod == "rs7_follow"), MAE:= rs7_follow_cv$results$MAE]

rs7_after_cv <- train(Incidence~Treatment+Triplet+log(Hist3yr)+log(Baseline), family = quasipoisson, data=rbct_allbreakdowns_after, method = "glm",
                      trControl = train.control)
all_cv_mods_after[which(mod == "rs7_after"), RMSE:= rs7_after_cv$results$RMSE]
all_cv_mods_after[which(mod == "rs7_after"), Rsquared:= rs7_after_cv$results$Rsquared]
all_cv_mods_after[which(mod == "rs7_after"), MAE:= rs7_after_cv$results$MAE]


LOOCV_7a <- lapply(1:nrow(rbct_allbreakdowns), function(x){
  m7a <- glmmTMB(round(Incidence/0.85)~Treatment+Triplet+log(Hist3yr)+log(hrdyrsrisk),family=genpois, 
                 data=rbct_allbreakdowns[-x,])
  return(predict(m7a, rbct_allbreakdowns[x,], type = "response"))
})
sq_errors_7a <- rep(0, length(LOOCV_7a))
abs_errors_7a <- rep(0, length(LOOCV_7a))

for(i in 1:length(LOOCV_7a)){
  sq_errors_7a[i] <- (LOOCV_7a[[i]]-rbct_allbreakdowns$Incidence[i]/0.85)^2
  abs_errors_7a[i] <- abs(LOOCV_7a[[i]]-rbct_allbreakdowns$Incidence[i]/0.85)
}
all_cv_mods[which(mod == "rs7a"), RMSE:= sqrt(mean(sq_errors_7a))]
all_cv_mods[which(mod == "rs7a"), Rsquared:= r2_func(rbct_allbreakdowns$Incidence/0.85, unlist(LOOCV_7a))]
all_cv_mods[which(mod == "rs7a"), MAE:= mean(abs_errors_7a)]
all_cv_mods[which(mod == "rs7a"), Akaike:= AICc(rs7a)]
all_cv_mods[which(mod == "rs7a"), BIC_val:= BIC(rs7a)]



LOOCV_7a_follow <- lapply(1:nrow(rbct_allbreakdowns_follow), function(x){
  m7a_follow <- glmmTMB(round(Incidence/0.85)~Treatment+Triplet+log(Hist3yr)+log(hdyrsrisk),family=genpois, 
                        data=rbct_allbreakdowns_follow[-x,])
  return(predict(m7a_follow, rbct_allbreakdowns_follow[x,], type = "response"))
})
sq_errors_7a_follow <- rep(0, length(LOOCV_7a_follow))
abs_errors_7a_follow <- rep(0, length(LOOCV_7a_follow))

for(i in 1:length(LOOCV_7a_follow)){
  sq_errors_7a_follow[i] <- (LOOCV_7a_follow[[i]]-rbct_allbreakdowns_follow$Incidence[i]/0.85)^2
  abs_errors_7a_follow[i] <- abs(LOOCV_7a_follow[[i]]-rbct_allbreakdowns_follow$Incidence[i]/0.85)
}
all_cv_mods_follow[which(mod == "rs7a_follow"), RMSE:= sqrt(mean(sq_errors_7a_follow))]
all_cv_mods_follow[which(mod == "rs7a_follow"), Rsquared:= r2_func(rbct_allbreakdowns_follow$Incidence/0.85, unlist(LOOCV_7a_follow))]
all_cv_mods_follow[which(mod == "rs7a_follow"), MAE:= mean(abs_errors_7a_follow)]
all_cv_mods_follow[which(mod == "rs7a_follow"), Akaike:= AICc(rs7a_follow)]
all_cv_mods_follow[which(mod == "rs7a_follow"), BIC_val:= BIC(rs7a_follow)]
all_cv_mods_follow


rs7a_after<-glmmTMB(round(Incidence/0.85)~Treatment+Triplet+log(Hist3yr)+log(hdyrsrisk),family=genpois, data= rbct_allbreakdowns_after)
rs7a_after2<-glmmTMB((Incidence)~Treatment+Triplet+log(Hist3yr)+log(hdyrsrisk),family=genpois, data= rbct_allbreakdowns_after)

summary(rs7a_after)
AICc(rs7a_after)
plot(simulateResiduals(rs7a_after))
check_model(rs7a_after2)
summary(rs7a_after$sdr, "fixed")[2,1]
rs7a_after_se <- summary(rs7a_after$sdr, "fixed")[2,2]
exp(summary(rs7a_after$sdr, "fixed")[2,1])-1

exp(summary(rs7a_after$sdr, "fixed")[2,1]+qnorm(0.975) * rs7a_after_se)-1
exp(summary(rs7a_after$sdr, "fixed")[2,1]-qnorm(0.975) * rs7a_after_se)-1


LOOCV_7a_after <- lapply(1:nrow(rbct_allbreakdowns_after), function(x){
  m7a_after <- glmmTMB(round(Incidence/0.85)~Treatment+Triplet+log(Hist3yr)+log(hdyrsrisk),family=genpois, 
                       data=rbct_allbreakdowns_after[-x,])
  return(predict(m7a_after, rbct_allbreakdowns_after[x,], type = "response"))
})
sq_errors_7a_after <- rep(0, length(LOOCV_7a_after))
abs_errors_7a_after <- rep(0, length(LOOCV_7a_after))

for(i in 1:length(LOOCV_7a_after)){
  sq_errors_7a_after[i] <- (LOOCV_7a_after[[i]]-rbct_allbreakdowns_after$Incidence[i]/0.85)^2
  abs_errors_7a_after[i] <- abs(LOOCV_7a_after[[i]]-rbct_allbreakdowns_after$Incidence[i]/0.85)
}
all_cv_mods_after[which(mod == "rs7a_after"), RMSE:= sqrt(mean(sq_errors_7a_after))]
all_cv_mods_after[which(mod == "rs7a_after"), Rsquared:= r2_func(rbct_allbreakdowns_after$Incidence/0.85, unlist(LOOCV_7a_after))]
all_cv_mods_after[which(mod == "rs7a_after"), MAE:= mean(abs_errors_7a_after)]
all_cv_mods_after[which(mod == "rs7a_after"), Akaike:= AICc(rs7a_after)]
all_cv_mods_after[which(mod == "rs7a_after"), BIC_val:= BIC(rs7a_after)]






#Model 7b-----
rs7b<-glmmTMB(round(Incidence/0.85)~log(Hist3yr)+log(hrdyrsrisk),
              family=genpois, data=rbct_allbreakdowns)
rs7b_follow<-glmmTMB(round(Incidence/0.85)~log(Hist3yr)+log(hdyrsrisk),
                     family=genpois, data=rbct_allbreakdowns_follow) #Without 0.85 division
rs7b_follow2<-glmmTMB((Incidence)~log(Hist3yr)+log(hdyrsrisk),
                      family=genpois, data=rbct_allbreakdowns_follow)

rs7b_after<-glmmTMB(round(Incidence/0.85)~log(Hist3yr)+log(hdyrsrisk),
                    family=genpois, data=rbct_allbreakdowns_after)
rs7b_after2<-glmmTMB((Incidence)~log(Hist3yr)+log(hdyrsrisk),
                     family=genpois, data=rbct_allbreakdowns_after)

summary(rs7b)
plot(simulateResiduals(rs7b_follow))
check_model(rs7b_follow)

check_model(rs7b_after)

check_model(rs7b_after)
check_model(rs7b2)

plot(simulateResiduals(rs7b_after))

LOOCV_7b <- lapply(1:nrow(rbct_allbreakdowns), function(x){
  m7b <- glmmTMB(round(Incidence/0.85)~log(Hist3yr)+log(hrdyrsrisk),family=genpois, 
                 data=rbct_allbreakdowns[-x,])
  return(predict(m7b, rbct_allbreakdowns[x,], type = "response"))
})
sq_errors_7b <- rep(0, length(LOOCV_7b))
abs_errors_7b <- rep(0, length(LOOCV_7b))

for(i in 1:length(LOOCV_7b)){
  sq_errors_7b[i] <- (LOOCV_7b[[i]]-rbct_allbreakdowns$Incidence[i]/0.85)^2
  abs_errors_7b[i] <- abs(LOOCV_7b[[i]]-rbct_allbreakdowns$Incidence[i]/0.85)
}
all_cv_mods[which(mod == "rs7b"), RMSE:= sqrt(mean(sq_errors_7b))]
all_cv_mods[which(mod == "rs7b"), Rsquared:= r2_func(round(rbct_allbreakdowns$Incidence/0.85), unlist(LOOCV_7b))]
all_cv_mods[which(mod == "rs7b"), MAE:= mean(abs_errors_7b)]
all_cv_mods[which(mod == "rs7b"), Akaike:= AICc(rs7b)]
all_cv_mods[which(mod == "rs7b"), BIC_val:= BIC(rs7b)]

LOOCV_7b_follow <- lapply(1:nrow(rbct_allbreakdowns_follow), function(x){
  m7b_follow <- glmmTMB(round(Incidence/0.85)~log(Hist3yr)+log(hdyrsrisk),family=genpois, 
                        data=rbct_allbreakdowns_follow[-x,])
  return(predict(m7b_follow, rbct_allbreakdowns_follow[x,], type = "response"))
})
sq_errors_7b_follow <- rep(0, length(LOOCV_7b_follow))
abs_errors_7b_follow <- rep(0, length(LOOCV_7b_follow))

for(i in 1:length(LOOCV_7b_follow)){
  sq_errors_7b_follow[i] <- (LOOCV_7b_follow[[i]]- round(rbct_allbreakdowns_follow$Incidence[i]/0.85))^2
  abs_errors_7b_follow[i] <- abs(LOOCV_7b_follow[[i]]-round(rbct_allbreakdowns_follow$Incidence[i]/0.85))
}
all_cv_mods_follow[which(mod == "rs7b_follow"), RMSE:= sqrt(mean(sq_errors_7b_follow))]
all_cv_mods_follow[which(mod == "rs7b_follow"), Rsquared:= r2_func(round(rbct_allbreakdowns_follow$Incidence/0.85), unlist(LOOCV_7b_follow))]
all_cv_mods_follow[which(mod == "rs7b_follow"), MAE:= mean(abs_errors_7b_follow)]
all_cv_mods_follow[which(mod == "rs7b_follow"), Akaike:= AICc(rs7b_follow)]
all_cv_mods_follow[which(mod == "rs7b_follow"), BIC_val:= BIC(rs7b_follow)]



LOOCV_7b_after <- lapply(1:nrow(rbct_allbreakdowns_after), function(x){
  m7b_after <- glmmTMB(round(Incidence/0.85)~log(Hist3yr)+log(hdyrsrisk),family=genpois, 
                       data=rbct_allbreakdowns_after[-x,])
  return(predict(m7b_after, rbct_allbreakdowns_after[x,], type = "response"))
})
sq_errors_7b_after <- rep(0, length(LOOCV_7b_after))
abs_errors_7b_after <- rep(0, length(LOOCV_7b_after))

for(i in 1:length(LOOCV_7b_after)){
  sq_errors_7b_after[i] <- (LOOCV_7b_after[[i]]-rbct_allbreakdowns_after$Incidence[i]/0.85)^2
  abs_errors_7b_after[i] <- abs(LOOCV_7b_after[[i]]-rbct_allbreakdowns_after$Incidence[i]/0.85)
}
all_cv_mods_after[which(mod == "rs7b_after"), RMSE:= sqrt(mean(sq_errors_7b_after))]
all_cv_mods_after[which(mod == "rs7b_after"), Rsquared:= r2_func(rbct_allbreakdowns_after$Incidence/0.85, unlist(LOOCV_7b_after))]
all_cv_mods_after[which(mod == "rs7b_after"), MAE:= mean(abs_errors_7b_after)]
all_cv_mods_after[which(mod == "rs7b_after"), Akaike:= AICc(rs7b_after)]
all_cv_mods_after[which(mod == "rs7b_after"), BIC_val:= BIC(rs7b_after)]

#Model 8 -----
rs8 <- glm(round(Incidence/0.85)~Treatment+Triplet+log(Hist3yr)+log(Baseline), family = poisson, data = rbct_allbreakdowns)
rs8_follow <- glm(round(Incidence/0.85)~Treatment+Triplet+log(Hist3yr)+log(Baseline), family = poisson, data = rbct_allbreakdowns_follow)
rs8_after <- glm(round(Incidence/0.85)~Treatment+Triplet+log(Hist3yr)+log(Baseline), family = poisson, data = rbct_allbreakdowns_after)

summary(rs8_follow)
summary(rs8_after)

rs8_cv <- train(round(Incidence/0.85)~Treatment+Triplet+log(Hist3yr)+log(Baseline), family = poisson, data=rbct_allbreakdowns, method = "glm",
                trControl = train.control)
all_cv_mods[which(mod == "rs8"), RMSE:= rs8_cv$results$RMSE]
all_cv_mods[which(mod == "rs8"), Rsquared:= rs8_cv$results$Rsquared]
all_cv_mods[which(mod == "rs8"), MAE:= rs8_cv$results$MAE]
all_cv_mods[which(mod == "rs8"), Akaike:= AICc(rs8)]
all_cv_mods[which(mod == "rs8"), BIC_val:= BIC(rs8)]



rs8_follow_cv <- train(Incidence~Treatment+Triplet+log(Hist3yr)+log(Baseline), family = poisson, data=rbct_allbreakdowns_follow, method = "glm",
                       trControl = train.control)

all_cv_mods_follow[which(mod == "rs8_follow"), RMSE:= rs8_follow_cv$results$RMSE]
all_cv_mods_follow[which(mod == "rs8_follow"), Rsquared:= rs8_follow_cv$results$Rsquared]
all_cv_mods_follow[which(mod == "rs8_follow"), MAE:= rs8_follow_cv$results$MAE]
all_cv_mods_follow[which(mod == "rs8_follow"), Akaike:= AICc(rs8_follow)]
all_cv_mods_follow[which(mod == "rs8_follow"), BIC_val:= BIC(rs8_follow)]
all_cv_mods_follow



rs8_after_cv <- train(round(Incidence/0.85)~Treatment+Triplet+log(Hist3yr)+log(Baseline), family = poisson, data=rbct_allbreakdowns_after, method = "glm",
                      trControl = train.control)

all_cv_mods_after[which(mod == "rs8_after"), RMSE:= rs8_after_cv$results$RMSE]
all_cv_mods_after[which(mod == "rs8_after"), Rsquared:= rs8_after_cv$results$Rsquared]
all_cv_mods_after[which(mod == "rs8_after"), MAE:= rs8_after_cv$results$MAE]
all_cv_mods_after[which(mod == "rs8_after"), Akaike:= AICc(rs8_after)]
all_cv_mods_after[which(mod == "rs8_after"), BIC_val:= BIC(rs8_after)]
all_cv_mods_after



