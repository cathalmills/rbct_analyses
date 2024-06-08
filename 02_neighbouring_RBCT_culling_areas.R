# install.packages("DHARMa")

library(glmmTMB)
library(MuMIn)
library(rstanarm)
library(bridgesampling)
library(DHARMa)
require(AER)
require(vcd)
require(performance)
require(cowplot)
require(rstan)
require(boot)
require(olsrr)
require(bc)
require(lme4)
require(caret)
require(bayesplot)
require(gridExtra)
require(vctrs)
require(purrr)
require(loo)


#Directory Management etc ----
data.dir2 <- ""
out.dir2 <- ""

plot_val_cols <- c("#2196F3", "red")

#Load in data ----
rbctconf_outside <- data.table(read_xlsx(file.path(data.dir2,"outside_areas_rbct.xlsx"), sheet = 1))

rbctconf_outside_follow <- data.table(read_xlsx(file.path(data.dir2,"outside_areas_rbct.xlsx"), sheet = 2))

rbctconf_outside_after <- data.table(read_xlsx(file.path(data.dir2,"outside_areas_rbct.xlsx"), sheet = 3))


rbctconf_outside$Triplet<-relevel(factor(rbctconf_outside$Triplet), ref="J")
rbctconf_outside$Treatment<-relevel(factor(rbctconf_outside$Treatment), ref="Survey-only")

rbctconf_outside_follow$Triplet<-relevel(factor(rbctconf_outside_follow$Triplet), ref="J")
rbctconf_outside_follow$Treatment<-relevel(factor(rbctconf_outside_follow$Treatment), ref="Survey-only")


rbctconf_outside_after$Triplet<-relevel(factor(rbctconf_outside_after$Triplet), ref="J")
rbctconf_outside_after$Treatment<-relevel(factor(rbctconf_outside_after$Treatment), ref="Survey-only")


#FREQUENTIST ANALYSIS -------

#Model from 2006 paper----
#Fit to initial period

rs1_outside<-glm(Incidence~Treatment+Triplet+log(Hist3yr)+log(Baseline), family = poisson, data=rbctconf_outside)
summary(rs1_outside)
plot(simulateResiduals(rs1_outside))
simulateResiduals(rs1_outside)
testResiduals(rs1_outside)

#Fit to period from follow-up cull onwards
rs1_outside_follow<-glm(Incidence~Treatment+Triplet+log(Hist3yr)+log(Baseline), family = poisson, data=rbctconf_outside_follow)
summary(rs1_outside_follow)

#Fit to post-trial period
rs1_outside_after<-glm(Incidence~Treatment+Triplet+log(Hist3yr)+log(Baseline), family = poisson, data=rbctconf_outside_after)
summary(rs1_outside_after)


rs1_outside_se <- sqrt(diag(vcov(rs1_outside)))*sqrt(rs1_outside$deviance/rs1_outside$df.residual)
rs1_outside_follow_se <- sqrt(diag(vcov(rs1_outside_follow)))*sqrt(rs1_outside_follow$deviance/rs1_outside_follow$df.residual)
rs1_outside_after_se <- sqrt(diag(vcov(rs1_outside_after))) * sqrt(rs1_outside_after$deviance/rs1_outside_after$df.residual)


exp(rs1_outside$coefficients[2])-1
exp(rs1_outside$coefficients[2]-qnorm(0.975) * rs1_outside_se[2])-1
exp(rs1_outside$coefficients[2]+qnorm(0.975) * rs1_outside_se[2])-1

exp(rs1_outside_follow$coefficients[2])-1
exp(rs1_outside_follow$coefficients[2]-qnorm(0.975) * rs1_outside_follow_se[2])-1
exp(rs1_outside_follow$coefficients[2]+qnorm(0.975) * rs1_outside_follow_se[2])-1

exp(rs1_outside_after$coefficients[2])-1
exp(rs1_outside_after$coefficients[2]-qnorm(0.975) * rs1_outside_after_se[2])-1
exp(rs1_outside_after$coefficients[2]+qnorm(0.975) * rs1_outside_after_se[2])-1

mod_vector <- c("rs1", "rs2", "rs3", "rs3a", "rs4", "rs4a", "rs4c",
                "rs5", "rs5a")
# Set up LOOCV
train.control <- trainControl(method = "LOOCV")

#Data.table will store model names and LOOCV accuracy metrics
cv_mods_outside <- data.table(mod = mod_vector, 
                                             RMSE = rep(0, length(mod_vector)),
                                             MAE = rep(0, length(mod_vector)))
rs1_outside_cv <- train(Incidence~Treatment+Triplet+log(Hist3yr)+log(Baseline), family = poisson, data=rbctconf_outside, method = "glm",
                        trControl = train.control)
cv_mods_outside[which(mod == "rs1_outside"), RMSE:= rs1_outside_cv$results$RMSE]
cv_mods_outside[which(mod == "rs1_outside"), MAE:= rs1_outside_cv$results$MAE]

mod_vector_follow <- c("rs1_outside_follow", "rs2_outside_follow", "rs3_follow", 
                       "rs3a_outside_follow", "rs4_follow", "rs4a_follow", "rs4c_follow",
                       "rs5_follow", "rs5a_follow")
cv_mods_outside_follow <- data.table(mod = mod_vector_follow, 
                                                    RMSE = rep(0, length(mod_vector_follow)),
                                                    MAE = rep(0, length(mod_vector_follow)))

rs1_outside_follow_cv <- train(Incidence~Treatment+Triplet+log(Hist3yr)+log(Baseline), family = poisson, data=rbctconf_outside_follow, method = "glm",
                               trControl = train.control)
cv_mods_outside_follow
cv_mods_outside_follow[which(mod == "rs1_outside_follow"), RMSE:= rs1_outside_follow_cv$results$RMSE]
cv_mods_outside_follow[which(mod == "rs1_outside_follow"), MAE:= rs1_outside_follow_cv$results$MAE]
cv_mods_outside_follow
mod_vector_after <- c("rs1_outside_after", "rs2_outside_after", "rs3_after", 
                      "rs3a_outside_after", "rs4_after",
                      "rs4a_after", "rs4c_after",
                      "rs5_after", "rs5a_after")
cv_mods_outside_after <- data.table(mod = mod_vector_after, 
                                                   RMSE = rep(0, length(mod_vector_after)),
                                                   MAE = rep(0, length(mod_vector_after)))

rs1_outside_after_cv <- train(Incidence~Treatment+Triplet+log(Hist3yr)+log(Baseline), family = poisson, data=rbctconf_outside_after, method = "glm",
                              trControl = train.control)
cv_mods_outside_after[which(mod == "rs1_outside_after"), RMSE:= rs1_outside_after_cv$results$RMSE]
cv_mods_outside_after[which(mod == "rs1_outside_after"), MAE:= rs1_outside_after_cv$results$MAE]
cv_mods_outside_after

performance::check_model(rs1_outside_after)

plot(simulateResiduals(rs1_outside))
plot(simulateResiduals(rs1_outside_after))
plot(simulateResiduals(rs1_outside_after))



BIC_rs1_outside <- BIC(rs1_outside)
BIC_rs1_outside
AICc_rs1_outside <- AICc(rs1_outside)
AICc_rs1_outside

BIC_rs1_outside_follow <- BIC(rs1_outside_follow)
BIC_rs1_outside_follow
AICc_rs1_outside_follow <- AICc(rs1_outside_follow)
AICc_rs1_outside_follow

BIC_rs1_outside_after <- BIC(rs1_outside_after)
BIC_rs1_outside_after
AICc_rs1_outside_after <- AICc(rs1_outside_after)
AICc_rs1_outside_after

#Model 2-----
rs2_outside<-glm(Incidence~Treatment+Triplet+log(Hist3yr)+offset(log(Baseline)), family = quasipoisson, data=rbctconf_outside)
check_model(rs2_outside)

rs2_outside_follow<-glm(Incidence~Treatment+Triplet+log(Hist3yr)+offset(log(Baseline)), family = quasipoisson, data=rbctconf_outside_follow)

rs2_outside_after<-glm(Incidence~Treatment+Triplet+log(Hist3yr)+offset(log(Baseline)), family = quasipoisson, data=rbctconf_outside_after)
summary(rs2_outside_after)
check_model(rs2_outside_after)


rs2_outside_cv <- train(Incidence~Treatment+Triplet+log(Hist3yr)+offset(log(Baseline)), family = quasipoisson, data=rbctconf_outside, method = "glm",
                        trControl = train.control)
cv_mods_outside[which(mod == "rs2_outside"), RMSE:= rs2_outside_cv$results$RMSE]
cv_mods_outside[which(mod == "rs2_outside"), MAE:= rs2_outside_cv$results$MAE]
cv_mods_outside

rs2_outside_follow_cv <- train(Incidence~Treatment+Triplet+log(Hist3yr)+offset(log(Baseline)), family = quasipoisson, data=rbctconf_outside_follow, method = "glm",
                               trControl = train.control)
cv_mods_outside_follow[which(mod == "rs2_outside_follow"), RMSE:= rs2_outside_follow_cv$results$RMSE]
cv_mods_outside_follow[which(mod == "rs2_outside_follow"), MAE:= rs2_outside_follow_cv$results$MAE]
cv_mods_outside_follow

rs2_outside_after_cv <- train(Incidence~Treatment+Triplet+log(Hist3yr)+offset(log(Baseline)), family = quasipoisson, data=rbctconf_outside_after, method = "glm",
                              trControl = train.control)
cv_mods_outside_after[which(mod == "rs2_outside_after"), RMSE:= rs2_outside_after_cv$results$RMSE]
cv_mods_outside_after[which(mod == "rs2_outside_after"), MAE:= rs2_outside_after_cv$results$MAE]
cv_mods_outside_after


rs2_outside_se <- sqrt(diag(vcov(rs2_outside)))
exp(rs2_outside$coefficients[2])-1
exp(rs2_outside$coefficients[2]-qnorm(0.975) * rs2_outside_se[2])-1
exp(rs2_outside$coefficients[2]+qnorm(0.975) * rs2_outside_se[2])-1

rs2_outside_follow_se <- sqrt(diag(vcov(rs2_outside_follow)))
exp(rs2_outside_follow$coefficients[2])-1
exp(rs2_outside_follow$coefficients[2]-qnorm(0.975) * rs2_outside_follow_se[2])-1
exp(rs2_outside_follow$coefficients[2]+qnorm(0.975) * rs2_outside_follow_se[2])-1

rs2_outside_after_se <- sqrt(diag(vcov(rs2_outside_after)))
exp(rs2_outside_after$coefficients[2])-1
exp(rs2_outside_after$coefficients[2]-qnorm(0.975) * rs2_outside_after_se[2])-1
exp(rs2_outside_after$coefficients[2]+qnorm(0.975) * rs2_outside_after_se[2])-1



#Model 3----
rs3_outside <- glmmTMB(Incidence~Treatment+Triplet+log(Hist3yr)+log(Baseline),family=genpois, data=rbctconf_outside)
summary(rs3_outside)
check_model(rs3_outside) 

rs3_outside_follow<-glmmTMB(Incidence~Treatment+Triplet+log(Hist3yr)+log(Baseline),family=genpois, data=rbctconf_outside_follow)
summary(rs3_outside_follow)
check_model(rs3_outside_follow)
plot(simulateResiduals(rs3_outside_follow))

rs3_outside_after<-glmmTMB(Incidence~Treatment+Triplet+log(Hist3yr)+log(Baseline),family=genpois, data=rbctconf_outside_after)
DHARMa::testDispersion(rs3_outside_after)

plot(simulateResiduals(rs3_outside_after))
summary(rs3_outside_after_new)
check_model(rs3_outside_after)

rs3_outside_se <- summary(rs3_outside$sdr, "fixed")[2,2]
exp(summary(rs3_outside$sdr, "fixed")[2,1])-1
exp(summary(rs3_outside$sdr, "fixed")[2,1]+qnorm(0.975) * rs3_outside_se)-1
exp(summary(rs3_outside$sdr, "fixed")[2,1]-qnorm(0.975) * rs3_outside_se)-1

rs3_outside_follow_se <- summary(rs3_outside_follow$sdr, "fixed")[2,2]
exp(summary(rs3_outside_follow$sdr, "fixed")[2,1])-1
exp(summary(rs3_outside_follow$sdr, "fixed")[2,1]-qnorm(0.975) * rs3_outside_follow_se)-1
exp(summary(rs3_outside_follow$sdr, "fixed")[2,1]+qnorm(0.975) * rs3_outside_follow_se)-1

rs3_outside_after_se <- summary(rs3_outside_after$sdr, "fixed")[2,2]
exp(summary(rs3_outside_after$sdr, "fixed")[2,1])-1
exp(summary(rs3_outside_after$sdr, "fixed")[2,1]-qnorm(0.975) * rs3_outside_after_se)-1
exp(summary(rs3_outside_after$sdr, "fixed")[2,1]+qnorm(0.975) * rs3_outside_after_se)-1




BIC_rs3_outside <- BIC(rs3_outside)
BIC_rs3_outside
AICc_rs3_outside <- AICc(rs3_outside)
AICc_rs3_outside

BIC_rs3_outside_follow <- BIC(rs3_outside_follow)
BIC_rs3_outside_follow
AICc_rs3_outside_follow <- AICc(rs3_outside_follow)
AICc_rs3_outside_follow

BIC_rs3_outside_after <- BIC(rs3_outside_after)
BIC_rs3_outside_after
AICc_rs3_outside_after <- AICc(rs3_outside_after)
AICc_rs3_outside_after


LOOCV <- lapply(1:nrow(rbctconf_outside), function(x){
  m1 <- glmmTMB(Incidence~Treatment+Triplet+log(Hist3yr)+log(Baseline),
                family=genpois, data=rbctconf_outside[-x,])
  return(predict(m1, rbctconf_outside[x,], type = "response"))
})
sq_errors <- rep(0, length(LOOCV))
abs_errors <- rep(0, length(LOOCV))

for(i in 1:length(LOOCV)){
  sq_errors[i] <- (LOOCV[[i]]-rbctconf_outside$Incidence[i])^2
  abs_errors[i] <- abs(LOOCV[[i]]-rbctconf_outside$Incidence[i])
}
cv_mods_outside[which(mod == "rs3"), mod:= "rs3_outside"]
cv_mods_outside[which(mod == "rs3_outside"), RMSE:= sqrt(mean(sq_errors))]
cv_mods_outside[which(mod == "rs3_outside"), MAE:= mean(abs_errors)]
cv_mods_outside



LOOCV_3_follow <- lapply(1:nrow(rbctconf_outside_follow), function(x){
  m1_3_follow <- glmmTMB(Incidence~Treatment+Triplet+log(Hist3yr)+log(Baseline),
                         family=genpois, data=rbctconf_outside_follow[-x,])
  return(predict(m1_3_follow, rbctconf_outside_follow[x,], type = "response"))
})
sq_errors_3_follow <- rep(0, length(LOOCV_3_follow))
abs_errors_3_follow <- rep(0, length(LOOCV_3_follow))
# rbctconf_outside_follow <- subset(rbctconf_outside_follow, 
#                           select = colnames(rbctconf_outside_follow)[1:26])
for(i in 1:length(LOOCV_3_follow)){
  sq_errors_3_follow[i] <- (LOOCV_3_follow[[i]]-rbctconf_outside_follow$Incidence[i])^2
  abs_errors_3_follow[i] <- abs(LOOCV_3_follow[[i]]-rbctconf_outside_follow$Incidence[i])
}
cv_mods_outside_follow[which(mod == "rs3_follow"), mod:= "rs3_outside_follow"]
cv_mods_outside_follow[which(mod == "rs3_outside_follow"), RMSE:= sqrt(mean(sq_errors_3_follow))]
cv_mods_outside_follow[which(mod == "rs3_outside_follow"), MAE:= mean(abs_errors_3_follow)]
cv_mods_outside_follow

LOOCV_3_after <- lapply(1:nrow(rbctconf_outside_after), function(x){
  m1_3_after <- glmmTMB(Incidence~Treatment+Triplet+log(Hist3yr)+ log(Baseline),
                        family=genpois, data=rbctconf_outside_after[-x,])
  return(predict(m1_3_after, rbctconf_outside_after[x,], type = "response"))
})
sq_errors_3_after <- rep(0, length(LOOCV_3_after))
abs_errors_3_after <- rep(0, length(LOOCV_3_after))
for(i in 1:length(LOOCV_3_after)){
  sq_errors_3_after[i] <- (LOOCV_3_after[[i]]-rbctconf_outside_after$Incidence[i])^2
  abs_errors_3_after[i] <- abs(LOOCV_3_after[[i]]-rbctconf_outside_after$Incidence[i])
}
cv_mods_outside_after[which(mod == "rs3_after"), mod:= "rs3_outside_after"]
cv_mods_outside_after[which(mod == "rs3_outside_after"), RMSE:= sqrt(mean(sq_errors_3_after))]
cv_mods_outside_after[which(mod == "rs3_outside_after"), MAE:= mean(abs_errors_3_after)]
cv_mods_outside_after


#Model 3a ----
rs3a_outside<- glmmTMB(Incidence~log(Hist3yr)+log(Baseline),family=genpois, data=rbctconf_outside)
summary(rs3a_outside)
plot(simulateResiduals(rs3a_outside))
check_model(rs3a_outside) 
#Posterior Predictive check indicates some potential misfit

rs3a_outside_follow<-glmmTMB(Incidence~log(Hist3yr)+log(Baseline),family=genpois, data=rbctconf_outside_follow)
summary(rs3a_outside_follow)
plot(simulateResiduals(rs3a_outside_follow))
check_model(rs3a_outside_follow)

rs3a_outside_after<-glmmTMB(Incidence~log(Hist3yr)+log(Baseline),family=genpois, data=rbctconf_outside_after)

summary(rs3a_outside_after)
plot(simulateResiduals(rs3a_outside_after))
check_model(rs3a_outside_after)


BIC_rs3a_outside <- BIC(rs3a_outside)
BIC_rs3a_outside
AICc_rs3a_outside <- AICc(rs3a_outside)
AICc_rs3a_outside

BIC_rs3a_outside_follow <- BIC(rs3a_outside_follow)
BIC_rs3a_outside_follow
AICc_rs3a_outside_follow <- AICc(rs3a_outside_follow)
AICc_rs3a_outside_follow

BIC_rs3a_outside_after <- BIC(rs3a_outside_after)
BIC_rs3a_outside_after
AICc_rs3a_outside_after <- AICc(rs3a_outside_after)
AICc_rs3a_outside_after


LOOCV_3a <- lapply(1:nrow(rbctconf_outside), function(x){
  m3a <- glmmTMB(Incidence~log(Hist3yr)+log(Baseline),
                 family=genpois, data=rbctconf_outside[-x,])
  return(predict(m3a, rbctconf_outside[x,], type = "response"))
})
sq_errors_3a <- rep(0, length(LOOCV_3a))
abs_errors_3a <- rep(0, length(LOOCV_3a))

for(i in 1:length(LOOCV_3a)){
  sq_errors_3a[i] <- (LOOCV_3a[[i]]-rbctconf_outside$Incidence[i])^2
  abs_errors_3a[i] <- abs(LOOCV_3a[[i]]-rbctconf_outside$Incidence[i])
}
cv_mods_outside[which(mod == "rs3a_outside"), RMSE:= sqrt(mean(sq_errors_3a))]
cv_mods_outside[which(mod == "rs3a_outside"), MAE:= mean(abs_errors_3a)]
cv_mods_outside


LOOCV_3a_follow <- lapply(1:nrow(rbctconf_outside_follow), function(x){
  m3a_follow <- glmmTMB(Incidence~log(Hist3yr)+log(Baseline),
                        family=genpois, data=rbctconf_outside_follow[-x,])
  return(predict(m3a_follow, rbctconf_outside_follow[x,], type = "response"))
})
sq_errors_3a_follow <- rep(0, length(LOOCV_3a_follow))
abs_errors_3a_follow <- rep(0, length(LOOCV_3a_follow))

for(i in 1:length(LOOCV_3a_follow)){
  sq_errors_3a_follow[i] <- (LOOCV_3a_follow[[i]]-rbctconf_outside_follow$Incidence[i])^2
  abs_errors_3a_follow[i] <- abs(LOOCV_3a_follow[[i]]-rbctconf_outside_follow$Incidence[i])
}
cv_mods_outside_follow[which(mod == "rs3a_outside_follow"), RMSE:= sqrt(mean(sq_errors_3a_follow))]
cv_mods_outside_follow[which(mod == "rs3a_outside_follow"), MAE:= mean(abs_errors_3a_follow)]
cv_mods_outside_follow






LOOCV_3a_after <- lapply(1:nrow(rbctconf_outside_after), function(x){
  m3a_after <- glmmTMB(Incidence~log(Hist3yr)+log(Baseline),
                       family=genpois, data=rbctconf_outside_after[-x,])
  return(predict(m3a_after, rbctconf_outside_after[x,], type = "response"))
})
sq_errors_3a_after <- rep(0, length(LOOCV_3a_after))
abs_errors_3a_after <- rep(0, length(LOOCV_3a_after))

for(i in 1:length(LOOCV_3a_after)){
  sq_errors_3a_after[i] <- (LOOCV_3a_after[[i]]-rbctconf_outside_after$Incidence[i])^2
  abs_errors_3a_after[i] <- abs(LOOCV_3a_after[[i]]-rbctconf_outside_after$Incidence[i])
}
cv_mods_outside_after[which(mod == "rs3a_outside_after"), RMSE:= sqrt(mean(sq_errors_3a_after))]
cv_mods_outside_after[which(mod == "rs3a_outside_after"), MAE:= mean(abs_errors_3a_after)]
cv_mods_outside_after


rs3c_after<-glmmTMB(Incidence~Treatment + log(Hist3yr)+offset(log(Baseline))+Triplet,family=genpois, data=rbctconf_outside_after)

#Model 4 ----
rs4_outside<-glmmTMB(Incidence~Treatment+Triplet+
                       log(Hist3yr)+offset(log(hdyrsrisk)), family=genpois,data=rbctconf_outside)
summary(rs4_outside)
AICc(rs4_outside)
check_model(rs4_outside)
rs4_outside_se <- summary(rs4_outside$sdr, "fixed")[2,2]
summary(rs4_outside$sdr, "fixed")[2,1]
exp(summary(rs4_outside$sdr, "fixed")[2,1])-1

exp(summary(rs4_outside$sdr, "fixed")[2,1]+qnorm(0.975) * rs4_outside_se)-1
exp(summary(rs4_outside$sdr, "fixed")[2,1]-qnorm(0.975) * rs4_outside_se)-1

LOOCV_4 <- lapply(1:nrow(rbctconf_outside), function(x){
  m4 <- glmmTMB(Incidence~Treatment+Triplet+
                  log(Hist3yr)+offset(log(hdyrsrisk)), 
                family=genpois, data=rbctconf_outside[-x,])
  return(predict(m4, rbctconf_outside[x,], type = "response"))
})
sq_errors_4 <- rep(0, length(LOOCV_4))
abs_errors_4 <- rep(0, length(LOOCV_4))

for(i in 1:length(LOOCV_4)){
  sq_errors_4[i] <- (LOOCV_4[[i]]-rbctconf_outside$Incidence[i])^2
  abs_errors_4[i] <- abs(LOOCV_4[[i]]-rbctconf_outside$Incidence[i])
}
cv_mods_outside[which(mod == "rs4_outside"), RMSE:= sqrt(mean(sq_errors_4))]
cv_mods_outside[which(mod == "rs4_outside"), MAE:= mean(abs_errors_4)]


rs4_outside_follow<-glmmTMB(Incidence~Treatment+Triplet+
                              log(Hist3yr)+offset(log(hdyrsrisk)), family=genpois,data=rbctconf_outside_follow)
summary(rs4_outside_follow)

rs4_outside_follow_se <- summary(rs4_outside_follow$sdr, "fixed")[2,2]
exp(summary(rs4_outside_follow$sdr, "fixed")[2,1])-1
exp(summary(rs4_outside_follow$sdr, "fixed")[2,1]+qnorm(0.975) * rs4_outside_follow_se)-1
exp(summary(rs4_outside_follow$sdr, "fixed")[2,1]-qnorm(0.975) * rs4_outside_follow_se)-1

plot(simulateResiduals(rs4_outside_follow))
LOOCV_4_follow <- lapply(1:nrow(rbctconf_outside), function(x){
  m4_follow <- glmmTMB(Incidence~Treatment+Triplet+
                         log(Hist3yr)+offset(log(hdyrsrisk)), 
                       family=genpois, data=rbctconf_outside_follow[-x,])
  return(predict(m4_follow, rbctconf_outside_follow[x,], type = "response"))
})
sq_errors_4_follow <- rep(0, length(LOOCV_4_follow))
abs_errors_4_follow <- rep(0, length(LOOCV_4_follow))

for(i in 1:length(LOOCV_4_follow)){
  sq_errors_4_follow[i] <- (LOOCV_4_follow[[i]]-rbctconf_outside_follow$Incidence[i])^2
  abs_errors_4_follow[i] <- abs(LOOCV_4_follow[[i]]-rbctconf_outside_follow$Incidence[i])
}
# tmp <- data.table(mod = "rs4_outside_follow", RMSE = 0, MAE = 0, Rsquared = 0)
# cv_mods_outside_follow <- rbind(cv_mods_outside_follow, tmp)
cv_mods_outside_follow[which(mod == "rs4_follow"), mod:= "rs4_outside_follow"]
cv_mods_outside_follow[which(mod == "rs4_outside_follow"), RMSE:= sqrt(mean(sq_errors_4_follow))]
cv_mods_outside_follow[which(mod == "rs4_outside_follow"), MAE:= mean(abs_errors_4_follow)]






check_model(rs4_outside_follow)


rs4_outside_after<-glmmTMB(Incidence~Treatment+Triplet+
                             log(Hist3yr)+offset(log(hdyrsrisk)), family=genpois,data=rbctconf_outside_after)
summary(rs4_outside_after)
AICc(rs4_outside_after)
check_model(rs4_outside_after)
exp(summary(rs4_outside_after$sdr, "fixed")[2,1])-1
rs4_outside_after_se <- summary(rs4_outside_after$sdr, "fixed")[2,2]
exp(summary(rs4_outside_after$sdr, "fixed")[2,1]+qnorm(0.975) * rs4_outside_after_se)-1
exp(summary(rs4_outside_after$sdr, "fixed")[2,1]-qnorm(0.975) * rs4_outside_after_se)-1

plot(simulateResiduals(rs4_outside_after))

LOOCV_4_after <- lapply(1:nrow(rbctconf_outside), function(x){
  m4_after <- glmmTMB(Incidence~Treatment+Triplet+
                        log(Hist3yr)+offset(log(hdyrsrisk)), 
                      family=genpois, data=rbctconf_outside_after[-x,])
  return(predict(m4_after, rbctconf_outside_after[x,], type = "response"))
})
sq_errors_4_after <- rep(0, length(LOOCV_4_after))
abs_errors_4_after <- rep(0, length(LOOCV_4_after))

for(i in 1:length(LOOCV_4_after)){
  sq_errors_4_after[i] <- (LOOCV_4_after[[i]]-rbctconf_outside_after$Incidence[i])^2
  abs_errors_4_after[i] <- abs(LOOCV_4_after[[i]]-rbctconf_outside_after$Incidence[i])
}
cv_mods_outside_after[which(mod == "rs4_outside_after"), RMSE:= sqrt(mean(sq_errors_4_after))]
cv_mods_outside_after[which(mod == "rs4_outside_after"), MAE:= mean(abs_errors_4_after)]


BIC_rs4_outside <- BIC(rs4_outside)
BIC_rs4_outside
AICc_rs4_outside <- AICc(rs4_outside)
AICc_rs4_outside

BIC_rs4_outside_follow <- BIC(rs4_outside_follow)
BIC_rs4_outside_follow
AICc_rs4_outside_follow <- AICc(rs4_outside_follow)
AICc_rs4_outside_follow

BIC_rs4_outside_after <- BIC(rs4_outside_after)
BIC_rs4_outside_after
AICc_rs4_outside_after <- AICc(rs4_outside_after)
AICc_rs4_outside_after



#Model 4a----
rs4a_outside<-glmmTMB(Incidence~log(Hist3yr)+offset(log(hdyrsrisk)), family=genpois,data=rbctconf_outside)
summary(rs4a_outside)
AICc(rs4a_outside)
check_model(rs4a_outside)
plot(simulateResiduals(rs4a_outside))
head(rbctconf_outside)


plot(simulateResiduals(rs4a_outside, quantreg = TRUE))

LOOCV_4a <- lapply(1:nrow(rbctconf_outside), function(x){
  m4a <- glmmTMB(Incidence~log(Hist3yr)+offset(log(hdyrsrisk)),
                 family=genpois, data=rbctconf_outside[-x,])
  return(predict(m4a, rbctconf_outside[x,], type = "response"))
})
sq_errors_4a <- rep(0, length(LOOCV_4a))
abs_errors_4a <- rep(0, length(LOOCV_4a))

for(i in 1:length(LOOCV_4a)){
  sq_errors_4a[i] <- (LOOCV_4a[[i]]-rbctconf_outside$Incidence[i])^2
  abs_errors_4a[i] <- abs(LOOCV_4a[[i]]-rbctconf_outside$Incidence[i])
}
cv_mods_outside[which(mod == "rs4a_outside"), RMSE:= sqrt(mean(sq_errors_4a))]
cv_mods_outside[which(mod == "rs4a_outside"), MAE:= mean(abs_errors_4a)]
cv_mods_outside



rs4a_outside_follow<-glmmTMB(Incidence~log(Hist3yr)+offset(log(hdyrsrisk)), 
                             family=genpois,data=rbctconf_outside_follow)
summary(rs4a_outside_follow)
check_model(rs4a_outside_follow)
plot(simulateResiduals(rs4a_outside_follow))
plot(simulateResiduals(rs4a_outside_follow, quantreg = TRUE))

LOOCV_4a_follow <- lapply(1:nrow(rbctconf_outside_follow), function(x){
  m4a_follow <- glmmTMB(Incidence~log(Hist3yr)+offset(log(hdyrsrisk)),
                        family=genpois, data=rbctconf_outside_follow[-x,])
  return(predict(m4a_follow, rbctconf_outside_follow[x,], type = "response"))
})
sq_errors_4a_follow <- rep(0, length(LOOCV_4a_follow))
abs_errors_4a_follow <- rep(0, length(LOOCV_4a_follow))

for(i in 1:length(LOOCV_4a_follow)){
  sq_errors_4a_follow[i] <- (LOOCV_4a_follow[[i]]-rbctconf_outside_follow$Incidence[i])^2
  abs_errors_4a_follow[i] <- abs(LOOCV_4a_follow[[i]]-rbctconf_outside_follow$Incidence[i])
}
tmp <- data.table(mod = "rs4a_outside_follow", RMSE = 0, MAE = 0, Rsquared = 0)
cv_mods_outside_follow <- rbind(cv_mods_outside_follow, tmp)
cv_mods_outside_follow[which(mod == "rs4a_outside_follow"), RMSE:= sqrt(mean(sq_errors_4a))]
cv_mods_outside_follow[which(mod == "rs4a_outside_follow"), MAE:= mean(abs_errors_4a)]
cv_mods_outside_follow





rs4a_outside_after<-glmmTMB(Incidence~log(Hist3yr)+offset(log(hdyrsrisk)), 
                            family=genpois,data=rbctconf_outside_after)
plot(simulateResiduals(rs4a_outside_after, quantreg = TRUE))
check_model(rs4a_outside_after)
LOOCV_4a_after <- lapply(1:nrow(rbctconf_outside_after), function(x){
  m4a_after <- glmmTMB(Incidence~log(Hist3yr)+offset(log(hdyrsrisk)),
                       family=genpois, data=rbctconf_outside_after[-x,])
  return(predict(m4a_after, rbctconf_outside_after[x,], type = "response"))
})
sq_errors_4a_after <- rep(0, length(LOOCV_4a_after))
abs_errors_4a_after <- rep(0, length(LOOCV_4a_after))

for(i in 1:length(LOOCV_4a_after)){
  sq_errors_4a_after[i] <- (LOOCV_4a_after[[i]]-rbctconf_outside_after$Incidence[i])^2
  abs_errors_4a_after[i] <- abs(LOOCV_4a_after[[i]]-rbctconf_outside_after$Incidence[i])
}
tmp <- data.table(mod = "rs4a_outside_after", RMSE = 0, MAE = 0, Rsquared = 0)
cv_mods_outside_after <- rbind(cv_mods_outside_after, tmp)
cv_mods_outside_after[which(mod == "rs4a_outside_after"), RMSE:= sqrt(mean(sq_errors_4a_after))]
cv_mods_outside_after[which(mod == "rs4a_outside_after"), MAE:= mean(abs_errors_4a_after)]



BIC_rs4a_outside <- BIC(rs4a_outside)
BIC_rs4a_outside
AICc_rs4a_outside <- AICc(rs4a_outside)
AICc_rs4a_outside

BIC_rs4a_outside_follow <- BIC(rs4a_outside_follow)
BIC_rs4a_outside_follow
AICc_rs4a_outside_follow <- AICc(rs4a_outside_follow)
AICc_rs4a_outside_follow

BIC_rs4a_outside_after <- BIC(rs4a_outside_after)
BIC_rs4a_outside_after
AICc_rs4a_outside_after <- AICc(rs4a_outside_after)









#Model 4c----
rs4c_outside<-glmmTMB(Incidence~Treatment+Triplet+log(Hist3yr)+log(hdyrsrisk), family=genpois,data=rbctconf_outside)

rs4c_outside_se <- summary(rs4c_outside$sdr, "fixed")[2,2]
exp(summary(rs4c_outside$sdr, "fixed")[2,1])-1
exp(summary(rs4c_outside$sdr, "fixed")[2,1]+qnorm(0.975) * rs4c_outside_se)-1
exp(summary(rs4c_outside$sdr, "fixed")[2,1]-qnorm(0.975) * rs4c_outside_se)-1
plot(simulateResiduals(rs4c_outside, quantreg = TRUE))
check_model(rs4c_outside)

LOOCV_4c <- lapply(1:nrow(rbctconf_outside), function(x){
  m4c <- glmmTMB(Incidence ~ Treatment+Triplet+log(Hist3yr)+log(hdyrsrisk),
                 family=genpois, data=rbctconf_outside[-x,])
  return(predict(m4c, rbctconf_outside[x,], type = "response"))
})
sq_errors_4c <- rep(0, length(LOOCV_4c))
abs_errors_4c <- rep(0, length(LOOCV_4c))

for(i in 1:length(LOOCV_4c)){
  sq_errors_4c[i] <- (LOOCV_4c[[i]]-rbctconf_outside$Incidence[i])^2
  abs_errors_4c[i] <- abs(LOOCV_4c[[i]]-rbctconf_outside$Incidence[i])
}
cv_mods_outside[which(mod == "rs4c_outside"), RMSE:= sqrt(mean(sq_errors_4c))]
cv_mods_outside[which(mod == "rs4c_outside"), MAE:= mean(abs_errors_4c)]
cv_mods_outside







rs4c_outside_follow<-glmmTMB(Incidence~Treatment+Triplet+log(Hist3yr)+log(hdyrsrisk), family=genpois,data=rbctconf_outside_follow)
summary(rs4c_outside_follow)
AICc(rs4c_outside_follow)
rs4c_outside_follow_se <- summary(rs4c_outside_follow$sdr, "fixed")[2,2]
exp(summary(rs4c_outside_follow$sdr, "fixed")[2,1])-1
exp(summary(rs4c_outside_follow$sdr, "fixed")[2,1]+qnorm(0.975) * rs4c_outside_follow_se)-1
exp(summary(rs4c_outside_follow$sdr, "fixed")[2,1]-qnorm(0.975) * rs4c_outside_follow_se)-1
plot(simulateResiduals(rs4c_outside_follow, quantreg = TRUE))
check_model(rs4c_outside_follow)

LOOCV_4c_follow <- lapply(1:nrow(rbctconf_outside_follow), function(x){
  m4c_follow <- glmmTMB(Incidence ~ Treatment+Triplet+log(Hist3yr)+log(hdyrsrisk),
                        family=genpois, data=rbctconf_outside_follow[-x,])
  return(predict(m4c_follow, rbctconf_outside_follow[x,], type = "response"))
})
sq_errors_4c_follow <- rep(0, length(LOOCV_4c_follow))
abs_errors_4c_follow <- rep(0, length(LOOCV_4c_follow))

for(i in 1:length(LOOCV_4c_follow)){
  sq_errors_4c_follow[i] <- (LOOCV_4c_follow[[i]]-rbctconf_outside_follow$Incidence[i])^2
  abs_errors_4c_follow[i] <- abs(LOOCV_4c_follow[[i]]-rbctconf_outside_follow$Incidence[i])
}
cv_mods_outside_follow[which(mod == "rs4c_outside_follow"), RMSE:= sqrt(mean(sq_errors_4c_follow))]
cv_mods_outside_follow[which(mod == "rs4c_outside_follow"), MAE:= mean(abs_errors_4c_follow)]
cv_mods_outside_follow




rs4c_outside_after<-glmmTMB(Incidence~Treatment+Triplet+log(Hist3yr)+log(hdyrsrisk), family=genpois,data=rbctconf_outside_after)
summary(rs4c_outside_after)
AICc(rs4c_outside_after)
rs4c_outside_after_se <- summary(rs4c_outside_after$sdr, "fixed")[2,2]

exp(summary(rs4c_outside$sdr, "fixed")[2,1])-1
exp(summary(rs4c_outside_after$sdr, "fixed")[2,1])-1
exp(summary(rs4c_outside_after$sdr, "fixed")[2,1])-1


exp(summary(rs4c_outside_after$sdr, "fixed")[2,1]+qnorm(0.975) * rs4c_outside_after_se)-1
exp(summary(rs4c_outside_after$sdr, "fixed")[2,1]-qnorm(0.975) * rs4c_outside_after_se)-1
plot(simulateResiduals(rs4c_outside_after, quantreg = TRUE))
check_model(rs4c_outside_after)

LOOCV_4c_after <- lapply(1:nrow(rbctconf_outside_after), function(x){
  m4c_after <- glmmTMB(Incidence ~ Treatment+Triplet+log(Hist3yr)+log(hdyrsrisk),
                       family=genpois, data=rbctconf_outside_after[-x,])
  return(predict(m4c_after, rbctconf_outside_after[x,], type = "response"))
})
sq_errors_4c_after <- rep(0, length(LOOCV_4c_after))
abs_errors_4c_after <- rep(0, length(LOOCV_4c_after))

for(i in 1:length(LOOCV_4c_after)){
  sq_errors_4c_after[i] <- (LOOCV_4c_after[[i]]-rbctconf_outside_after$Incidence[i])^2
  abs_errors_4c_after[i] <- abs(LOOCV_4c_after[[i]]-rbctconf_outside_after$Incidence[i])
}
cv_mods_outside_after[which(mod == "rs4c_outside_after"), RMSE:= sqrt(mean(sq_errors_4c_after))]
cv_mods_outside_after[which(mod == "rs4c_outside_after"), MAE:= mean(abs_errors_4c_after)]
cv_mods_outside_after



BIC_rs4c_outside <- BIC(rs4c_outside)
BIC_rs4c_outside
AICc_rs4c_outside <- AICc(rs4c_outside)
AICc_rs4c_outside

BIC_rs4c_outside_after <- BIC(rs4c_outside_after)
BIC_rs4c_outside_after
AICc_rs4c_outside_after <- AICc(rs4c_outside_after)
AICc_rs4c_outside_after

BIC_rs4c_outside_after <- BIC(rs4c_outside_after)
BIC_rs4c_outside_after
AICc_rs4c_outside_after <- AICc(rs4c_outside_after)
AICc_rs4c_outside_after


#Model 4d ----
rs4d_outside<-glmmTMB(Incidence~log(Hist3yr)+log(hdyrsrisk), family=genpois,
                      data=rbctconf_outside)
summary(rs4d_outside)
check_model(rs4d_outside)
AICc(rs4d_outside)
plot(simulateResiduals(rs4d_outside))

LOOCV_4d <- lapply(1:nrow(rbctconf_outside), function(x){
  m4d <- glmmTMB(Incidence ~ log(Hist3yr)+log(hdyrsrisk),
                 family=genpois, data=rbctconf_outside[-x,])
  return(predict(m4d, rbctconf_outside[x,], type = "response"))
})
sq_errors_4d <- rep(0, length(LOOCV_4d))
abs_errors_4d <- rep(0, length(LOOCV_4d))

for(i in 1:length(LOOCV_4d)){
  sq_errors_4d[i] <- (LOOCV_4d[[i]]-rbctconf_outside$Incidence[i])^2
  abs_errors_4d[i] <- abs(LOOCV_4d[[i]]-rbctconf_outside$Incidence[i])
}
cv_mods_outside[which(mod == "rs4d_outside"), RMSE:= sqrt(mean(sq_errors_4d))]
cv_mods_outside[which(mod == "rs4d_outside"), MAE:= mean(abs_errors_4d)]
cv_mods_outside






rs4d_outside_follow<-glmmTMB(Incidence~log(Hist3yr)+log(hdyrsrisk), family=genpois,
                             data=rbctconf_outside_follow)
summary(rs4d_outside_follow)

check_model(rs4d_outside_follow)
AICc(rs4d_outside_follow)
plot(simulateResiduals(rs4d_outside))

LOOCV_4d_follow <- lapply(1:nrow(rbctconf_outside_follow), function(x){
  m4d_follow <- glmmTMB(Incidence ~ log(Hist3yr)+log(hdyrsrisk),
                        family=genpois, data=rbctconf_outside_follow[-x,])
  return(predict(m4d_follow, rbctconf_outside_follow[x,], type = "response"))
})
sq_errors_4d_follow <- rep(0, length(LOOCV_4d_follow))
abs_errors_4d_follow <- rep(0, length(LOOCV_4d_follow))

for(i in 1:length(LOOCV_4d_follow)){
  sq_errors_4d_follow[i] <- (LOOCV_4d_follow[[i]]-rbctconf_outside_follow$Incidence[i])^2
  abs_errors_4d_follow[i] <- abs(LOOCV_4d_follow[[i]]-rbctconf_outside_follow$Incidence[i])
}
tmp <- data.table(mod = "rs4d_outside_follow", RMSE = 0, MAE = 0, Rsquared = 0)
cv_mods_outside_follow <- rbind(cv_mods_outside_follow, tmp)
cv_mods_outside_follow[which(mod == "rs4d_outside_follow"), RMSE:= sqrt(mean(sq_errors_4d_follow))]
cv_mods_outside_follow[which(mod == "rs4d_outside_follow"), MAE:= mean(abs_errors_4d_follow)]
cv_mods_outside_follow





rs4d_outside_after<-glmmTMB(Incidence~log(Hist3yr)+log(hdyrsrisk), family=genpois,
                            data=rbctconf_outside_after)

summary(rs4d_outside_after)
check_model(rs4d_outside_after)
AICc(rs4d_outside_after)
plot(simulateResiduals(rs4d_outside))

LOOCV_4d_after <- lapply(1:nrow(rbctconf_outside_after), function(x){
  m4d_after <- glmmTMB(Incidence ~ log(Hist3yr)+log(hdyrsrisk),
                       family=genpois, data=rbctconf_outside_after[-x,])
  return(predict(m4d_after, rbctconf_outside_after[x,], type = "response"))
})
sq_errors_4d_after <- rep(0, length(LOOCV_4d_after))
abs_errors_4d_after <- rep(0, length(LOOCV_4d_after))

for(i in 1:length(LOOCV_4d_after)){
  sq_errors_4d_after[i] <- (LOOCV_4d_after[[i]]-rbctconf_outside_after$Incidence[i])^2
  abs_errors_4d_after[i] <- abs(LOOCV_4d_after[[i]]-rbctconf_outside_after$Incidence[i])
}
cv_mods_outside_after[which(mod == "rs4d_outside_after"), RMSE:= sqrt(mean(sq_errors_4d_after))]
cv_mods_outside_after[which(mod == "rs4d_outside_after"), MAE:= mean(abs_errors_4d_after)]
cv_mods_outside_after


BIC_rs4d_outside <- BIC(rs4d_outside)
BIC_rs4d_outside
AICc_rs4d_outside <- AICc(rs4d_outside)
AICc_rs4d_outside

BIC_rs4d_outside_after <- BIC(rs4d_outside_after)
BIC_rs4d_outside_after
AICc_rs4d_outside_after <- AICc(rs4d_outside_after)
AICc_rs4d_outside_after

BIC_rs4d_outside_after <- BIC(rs4d_outside_after)
BIC_rs4d_outside_after
AICc_rs4d_outside_after <- AICc(rs4d_outside_after)
AICc_rs4d_outside_after



#Model 5 ----
rs5_outside<-lm(Incidence/hdyrsrisk~Treatment+Triplet+log(Hist3yr), data=rbctconf_outside)
rs5_outside$coefficients*mean(rbctconf_outside$hdyrsrisk)
summary(rs5_outside)
# confint(rs5_outside)*mean(rbctconf_outside$hdyrsrisk)
confint(rs5_outside)
plot(rs5_outside)
rbctconf_outside
ols_test_normality(rs5_outside)
ols_plot_resid_fit(rs5_outside)
ols_plot_resid_hist(rs5_outside)
?shapiro.test
AICc(rs5_outside)
plot(simulateResiduals(rs5_outside))

LOOCV_5 <- lapply(1:nrow(rbctconf_outside), function(x){
  m5 <- lm(Incidence/hdyrsrisk~Treatment+Triplet+log(Hist3yr), 
           data=rbctconf_outside[-x,])
  return(predict(m5, rbctconf_outside[x,], type = "response"))
})
sq_errors_5 <- rep(0, length(LOOCV_5))
abs_errors_5 <- rep(0, length(LOOCV_5))
LOOCV_5[[i]]
for(i in 1:length(LOOCV_5)){
  sq_errors_5[i] <- (LOOCV_5[[i]]*rbctconf_outside$hdyrsrisk[i]-rbctconf_outside$Incidence[i])^2
  abs_errors_5[i] <- abs(LOOCV_5[[i]]*rbctconf_outside$hdyrsrisk[i]-rbctconf_outside$Incidence[i])
}
cv_mods_outside[which(mod == "rs5_outside"), RMSE:= sqrt(mean(sq_errors_5))]
cv_mods_outside[which(mod == "rs5_outside"), MAE:= mean(abs_errors_5)]
cv_mods_outside






rs5_outside_follow<-lm(Incidence/hdyrsrisk~Treatment+Triplet+log(Hist3yr), data=rbctconf_outside_follow)
summary(rs5_outside_follow)
confint(rs5_outside_follow)
AICc(rs5_outside_follow)
ols_test_normality(rs5_outside_follow)
ols_plot_resid_fit(rs5_outside_follow)
ols_plot_resid_hist(rs5_outside_follow)
colSums( influence.measures(rs5_outside_follow)$is.inf )
plot(rs5_outside_follow)
cv_mods_outside_follow

LOOCV_5_follow <- lapply(1:nrow(rbctconf_outside_follow), function(x){
  m5_follow <- lm(Incidence/hdyrsrisk~Treatment+Triplet+log(Hist3yr), 
                  data=rbctconf_outside_follow[-x,])
  return(predict(m5_follow, rbctconf_outside_follow[x,], type = "response"))
})
sq_errors_5_follow <- rep(0, length(LOOCV_5_follow))
abs_errors_5_follow <- rep(0, length(LOOCV_5_follow))
LOOCV_5_follow[[i]]
for(i in 1:length(LOOCV_5_follow)){
  sq_errors_5_follow[i] <- (LOOCV_5_follow[[i]]*rbctconf_outside_follow$hdyrsrisk[i]-rbctconf_outside_follow$Incidence[i])^2
  abs_errors_5_follow[i] <- abs(LOOCV_5_follow[[i]]*rbctconf_outside_follow$hdyrsrisk[i]-rbctconf_outside_follow$Incidence[i])
}
cv_mods_outside_follow[which(mod == "rs5_outside_follow"), RMSE:= sqrt(mean(sq_errors_5_follow))]
cv_mods_outside_follow[which(mod == "rs5_outside_follow"), MAE:= mean(abs_errors_5_follow)]




rs5_outside_after<-lm(Incidence/hdyrsrisk~Treatment+Triplet+log(Hist3yr), data=rbctconf_outside_after)
summary(rs5_outside_after)
confint(rs5_outside_after)
AICc(rs5_outside_after)
ols_test_normality(rs5_outside_after)
ols_plot_resid_fit(rs5_outside_after)
ols_plot_resid_hist(rs5_outside_after)
colSums( influence.measures(rs5_outside_after)$is.inf )

LOOCV_5_after <- lapply(1:nrow(rbctconf_outside_after), function(x){
  m5_after <- lm(Incidence/hdyrsrisk~Treatment+Triplet+log(Hist3yr), 
                 data=rbctconf_outside_after[-x,])
  return(predict(m5_after, rbctconf_outside_after[x,], type = "response"))
})
sq_errors_5_after <- rep(0, length(LOOCV_5_after))
abs_errors_5_after <- rep(0, length(LOOCV_5_after))
LOOCV_5_after[[i]]
for(i in 1:length(LOOCV_5_after)){
  sq_errors_5_after[i] <- (LOOCV_5_after[[i]]*rbctconf_outside_after$hdyrsrisk[i]-rbctconf_outside_after$Incidence[i])^2
  abs_errors_5_after[i] <- abs(LOOCV_5_after[[i]]*rbctconf_outside_after$hdyrsrisk[i]-rbctconf_outside_after$Incidence[i])
}
cv_mods_outside_after[which(mod == "rs5_outside_after"), RMSE:= sqrt(mean(sq_errors_5_after))]
cv_mods_outside_after[which(mod == "rs5_outside_after"), MAE:= mean(abs_errors_5_after)]
cv_mods_outside_after




BIC_rs5_outside <- BIC(rs5_outside)
BIC_rs5_outside
AICc_rs5_outside <- AICc(rs5_outside)
AICc_rs5_outside

BIC_rs5_outside_follow <- BIC(rs5_outside_follow)
BIC_rs5_outside_follow
AICc_rs5_outside_follow <- AICc(rs5_outside_follow)
AICc_rs5_outside_follow

BIC_rs5_outside_after <- BIC(rs5_outside_after)
BIC_rs5_outside_after
AICc_rs5_outside_after <- AICc(rs5_outside_after)
AICc_rs5_outside_after

#Model 5a ----
rs5a_outside<-lm(Incidence/hdyrsrisk~log(Hist3yr), data=rbctconf_outside)
summary(rs5a_outside)

ols_test_normality(rs5a_outside)
ols_plot_resid_fit(rs5a_outside)
ols_plot_resid_hist(rs5a_outside)
AICc(rs5a_outside)
plot(rs5a_outside)
plot(simulateResiduals(rs5a_outside))
influence.measures(rs5a_outside)$is.inf
colSums( influence.measures(rs5a_outside)$is.inf )
rbctconf_outside[11, ]
colSums( influence.measures(rs1_outside)$is.inf )

LOOCV_5a <- lapply(1:nrow(rbctconf_outside), function(x){
  m5a <- lm(Incidence/hdyrsrisk~log(Hist3yr), 
            data=rbctconf_outside[-x,])
  return(predict(m5a, rbctconf_outside[x,], type = "response"))
})
sq_errors_5a <- rep(0, length(LOOCV_5a))
abs_errors_5a <- rep(0, length(LOOCV_5a))
LOOCV_5a[[i]]
for(i in 1:length(LOOCV_5a)){
  sq_errors_5a[i] <- (LOOCV_5a[[i]]*rbctconf_outside$hdyrsrisk[i]-rbctconf_outside$Incidence[i])^2
  abs_errors_5a[i] <- abs(LOOCV_5a[[i]]*rbctconf_outside$hdyrsrisk[i]-rbctconf_outside$Incidence[i])
}
cv_mods_outside[which(mod == "rs5a_outside"), RMSE:= sqrt(mean(sq_errors_5a))]
cv_mods_outside[which(mod == "rs5a_outside"), MAE:= mean(abs_errors_5a)]
cv_mods_outside







rs5a_outside_after<-lm(Incidence/hdyrsrisk~log(Hist3yr), data=rbctconf_outside_after)
summary(rs5a_outside_after)

ols_test_normality(rs5a_outside_after)
ols_plot_resid_fit(rs5a_outside_after)
ols_plot_resid_hist(rs5a_outside_after)
AICc(rs5a_outside_after)
plot(rs5a_outside_after)
plot(simulateResiduals(rs5a_outside_after))
colSums( influence.measures(rs5a_outside_after)$is.inf )
influence.measures(rs5a_outside_after)

LOOCV_5a_after <- lapply(1:nrow(rbctconf_outside_after), function(x){
  m5a_after <- lm(Incidence/hdyrsrisk~log(Hist3yr), 
                  data=rbctconf_outside_after[-x,])
  return(predict(m5a_after, rbctconf_outside_after[x,], type = "response"))
})
sq_errors_5a_after <- rep(0, length(LOOCV_5a_after))
abs_errors_5a_after <- rep(0, length(LOOCV_5a_after))
LOOCV_5a_after[[i]]
for(i in 1:length(LOOCV_5a_after)){
  sq_errors_5a_after[i] <- (LOOCV_5a_after[[i]]*rbctconf_outside_after$hdyrsrisk[i]-rbctconf_outside_after$Incidence[i])^2
  abs_errors_5a_after[i] <- abs(LOOCV_5a_after[[i]]*rbctconf_outside_after$hdyrsrisk[i]-rbctconf_outside_after$Incidence[i])
}
cv_mods_outside_after[which(mod == "rs5a_outside_after"), RMSE:= sqrt(mean(sq_errors_5a_after))]
cv_mods_outside_after[which(mod == "rs5a_outside_after"), MAE:= mean(abs_errors_5a_after)]
cv_mods_outside_after








rs5a_outside_after<-lm(Incidence/hdyrsrisk~log(Hist3yr), data=rbctconf_outside_after)
summary(rs5a_outside_after)
ols_test_normality(rs5a_outside_after)
ols_plot_resid_fit(rs5a_outside_after)
ols_plot_resid_hist(rs5a_outside_after)
AICc(rs5a_outside_after)
plot(rs5a_outside_after)
plot(simulateResiduals(rs5a_outside_after))
colSums(influence.measures(rs5a_outside_after)$is.inf )
influence.measures(rs5a_outside_after)
rbctconf_outside_after[11, ]

LOOCV_5a_after <- lapply(1:nrow(rbctconf_outside_after), function(x){
  m5a_after <- lm(Incidence/hdyrsrisk~log(Hist3yr), 
                  data=rbctconf_outside_after[-x,])
  return(predict(m5a_after, rbctconf_outside_after[x,], type = "response"))
})
sq_errors_5a_after <- rep(0, length(LOOCV_5a_after))
abs_errors_5a_after <- rep(0, length(LOOCV_5a_after))
LOOCV_5a_after[[i]]
for(i in 1:length(LOOCV_5a_after)){
  sq_errors_5a_after[i] <- (LOOCV_5a_after[[i]]*rbctconf_outside_after$hdyrsrisk[i]-rbctconf_outside_after$Incidence[i])^2
  abs_errors_5a_after[i] <- abs(LOOCV_5a_after[[i]]*rbctconf_outside_after$hdyrsrisk[i]-rbctconf_outside_after$Incidence[i])
}
cv_mods_outside_after[which(mod == "rs5a_outside_after"), RMSE:= sqrt(mean(sq_errors_5a_after))]
cv_mods_outside_after[which(mod == "rs5a_outside_after"), MAE:= mean(abs_errors_5a_after)]
cv_mods_outside_after[!duplicated(cv_mods_outside_after)]
cv_mods_outside_after <- cv_mods_outside_after[!duplicated(cv_mods_outside_after)]
cv_mods_outside_after

loo_rs




BIC_rs5a_outside <- BIC(rs5a_outside)
BIC_rs5a_outside
AICc_rs5a_outside <- AICc(rs5a_outside)
AICc_rs5a_outside

BIC_rs5a_outside_follow <- BIC(rs5a_outside_follow)
BIC_rs5a_outside_follow
AICc_rs5a_outside_follow <- AICc(rs5a_outside_follow)
AICc_rs5a_outside_follow

BIC_rs5a_outside_after <- BIC(rs5a_outside_after)
BIC_rs5a_outside_after
AICc_rs5a_outside_after <- AICc(rs5a_outside_after)
AICc_rs5a_outside_after

#Model 6b----
rs6b_outside<-glmmTMB(Incidence~offset(log(hdyrsrisk)),family=genpois, data=rbctconf_outside)
summary(rs6b_outside)
AICc(rs6b_outside)
check_model(rs6b_outside)
plot(simulateResiduals(rs6b_outside))

LOOCV_6b <- lapply(1:nrow(rbctconf_outside), function(x){
  m6b <- glmmTMB(Incidence~offset(log(hdyrsrisk)), 
                 family=genpois, data=rbctconf_outside[-x,])
  return(predict(m6b, rbctconf_outside[x,], type = "response"))
})
sq_errors_6b <- rep(0, length(LOOCV_6b))
abs_errors_6b <- rep(0, length(LOOCV_6b))

for(i in 1:length(LOOCV_6b)){
  sq_errors_6b[i] <- (LOOCV_6b[[i]]-rbctconf_outside$Incidence[i])^2
  abs_errors_6b[i] <- abs(LOOCV_6b[[i]]-rbctconf_outside$Incidence[i])
}

cv_mods_outside[which(mod == "rs6b_outside"), RMSE:= sqrt(mean(sq_errors_6b))]
cv_mods_outside[which(mod == "rs6b_outside"), MAE:= mean(abs_errors_6b)]
cv_mods_outside



plot(simulateResiduals(rs6b_outside, n = 2500, quantreg = TRUE))



rs6b_outside_follow<-glmmTMB(Incidence~offset(log(hdyrsrisk)), family=genpois,data=rbctconf_outside_follow)
summary(rs6b_outside_follow)
AICc(rs6b_outside_follow)
rs6b_outside_follow_se <- summary(rs6b_outside_follow$sdr, "fixed")[2,2]
check_model(rs6b_outside_follow)
plot(simulateResiduals(rs6b_outside_follow))

LOOCV_6b_follow <- lapply(1:nrow(rbctconf_outside), function(x){
  m6b_follow <- glmmTMB(Incidence~offset(log(hdyrsrisk)), 
                        family=genpois, data=rbctconf_outside_follow[-x,])
  return(predict(m6b_follow, rbctconf_outside_follow[x,], type = "response"))
})
sq_errors_6b_follow <- rep(0, length(LOOCV_6b_follow))
abs_errors_6b_follow <- rep(0, length(LOOCV_6b_follow))

for(i in 1:length(LOOCV_6b_follow)){
  sq_errors_6b_follow[i] <- (LOOCV_6b_follow[[i]]-rbctconf_outside_follow$Incidence[i])^2
  abs_errors_6b_follow[i] <- abs(LOOCV_6b_follow[[i]]-rbctconf_outside_follow$Incidence[i])
}
cv_mods_outside_follow[which(mod == "rs6b_outside_follow"), RMSE:= sqrt(mean(sq_errors_6b_follow))]
cv_mods_outside_follow[which(mod == "rs6b_outside_follow"), MAE:= mean(abs_errors_6b_follow)]
cv_mods_outside_follow



rs6b_outside_after<-glmmTMB(Incidence~offset(log(hdyrsrisk)),family=genpois, data=rbctconf_outside_after)
summary(rs6b_outside_after)
AICc(rs6b_outside_after)
check_model(rs6b_outside_after)

plot(simulateResiduals(rs6b_outside_after))

LOOCV_6b_after <- lapply(1:nrow(rbctconf_outside), function(x){
  m6b_after <- glmmTMB(Incidence~offset(log(hdyrsrisk)),family=genpois, data=rbctconf_outside_after[-x,])
  return(predict(m6b_after, rbctconf_outside_after[x,], type = "response"))
})
sq_errors_6b_after <- rep(0, length(LOOCV_6b_after))
abs_errors_6b_after <- rep(0, length(LOOCV_6b_after))

for(i in 1:length(LOOCV_6b_after)){
  sq_errors_6b_after[i] <- (LOOCV_6b_after[[i]]-rbctconf_outside_after$Incidence[i])^2
  abs_errors_6b_after[i] <- abs(LOOCV_6b_after[[i]]-rbctconf_outside_after$Incidence[i])
}
cv_mods_outside_after[which(mod == "rs6b_outside_after"), RMSE:= sqrt(mean(sq_errors_6b_after))]
cv_mods_outside_after[which(mod == "rs6b_outside_after"), MAE:= mean(abs_errors_6b_after)]
cv_mods_outside_after




BIC_rs6b_outside <- BIC(rs6b_outside)
BIC_rs6b_outside
AICc_rs6b_outside <- AICc(rs6b_outside)
AICc_rs6b_outside

BIC_rs6b_outside_follow <- BIC(rs6b_outside_follow)
BIC_rs6b_outside_follow
AICc_rs6b_outside_follow <- AICc(rs6b_outside_follow)
AICc_rs6b_outside_follow

BIC_rs6b_outside_after <- BIC(rs6b_outside_after)
BIC_rs6b_outside_after
AICc_rs6b_outside_after <- AICc(rs6b_outside_after)
AICc_rs6b_outside_after



#DIAGNOSTICS-----

#rs2_outside Influence----
rs2_outside_diagnostics <- check_model(rs2_outside, check = "outliers")
rs2_outside_diag_plots <- plot(rs2_outside_diagnostics, return_list = TRUE)
update_geom_defaults("text", list(size = 20))
?update_geom_defaults
update_geom_defaults("size", )
rs2_outside_diag_plots
rs2_outside_diag_plots
rs4d
grid.arrange(rs2_outside_diag_plots$OUTLIERS+ theme(text = element_text(size = 20),
                                                    axis.text.x = element_text(size=16),
                                                    axis.text.y = element_text(size=16),
                                                    plot.title = element_text(size=20, face = "bold", hjust = 0.5),
                                                    plot.subtitle = element_text(size=16, hjust = 0.5),
                                                    axis.title=element_text(size=18), 
                                                    legend.text=element_text(size=14)+geom_text(size = 20)))
dev.off()
rs2_outside_diag_plots_grob <- arrangeGrob(rs2_outside_diag_plots$OUTLIERS+ theme(text = element_text(size = 20),axis.text.x = element_text(size=20),
                                                                                  axis.text.y = element_text(size=20),
                                                                                  plot.title = element_text(size=24, face = "bold", hjust = 0.5),
                                                                                  plot.subtitle = element_text(size=20, hjust = 0.5),
                                                                                  axis.title=element_text(size=20), 
                                                                                  legend.text=element_text(size=18)),
                                           nrow = 1)
ggsave(rs2_outside_diag_plots_grob, file = file.path(out.dir2,"rs2_outside_influence_plot.pdf"), h = 14, w = 20)

#rs2_outside follow influence----
rs2_outside_follow_diagnostics <- check_model(rs2_outside_follow, check = "outliers")
rs2_outside_follow_diag_plots <- plot(rs2_outside_follow_diagnostics, return_list = TRUE)
rs2_outside_follow_diag_plots
grid.arrange(rs2_outside_follow_diag_plots$OUTLIERS+ theme(text = element_text(size = 20),
                                                           axis.text.x = element_text(size=16),
                                                           axis.text.y = element_text(size=16),
                                                           plot.title = element_text(size=20, face = "bold", hjust = 0.5),
                                                           plot.subtitle = element_text(size=16, hjust = 0.5),
                                                           axis.title=element_text(size=18), 
                                                           legend.text=element_text(size=14)+geom_text(size = 20)),
             nrow = 1)
rs2_outside_follow_diag_plots_grob <- arrangeGrob(rs2_outside_follow_diag_plots$OUTLIERS+ theme(axis.text.x = element_text(size=20),
                                                                                                axis.text.y = element_text(size=20),
                                                                                                plot.title = element_text(size=24, face = "bold", hjust = 0.5),
                                                                                                plot.subtitle = element_text(size=20, hjust = 0.5),
                                                                                                axis.title=element_text(size=20), 
                                                                                                legend.text=element_text(size=18)),
                                                  nrow = 1)
dev.off()
ggsave(rs2_outside_follow_diag_plots_grob, file = file.path(out.dir2,"rs2_outside_follow_influence_plot.pdf"), h = 14, w = 20)

#rs2_outside after influence-----
rs2_outside_after_diagnostics <- check_model(rs2_outside_after, check = "outliers")
rs2_outside_after_diag_plots <- plot(rs2_outside_after_diagnostics, return_list = TRUE)
rs2_outside_after_diag_plots
grid.arrange(rs2_outside_after_diag_plots+ theme(text = element_text(size = 20),
                                                 axis.text.x = element_text(size=16),
                                                 axis.text.y = element_text(size=16),
                                                 plot.title = element_text(size=20, face = "bold", hjust = 0.5),
                                                 plot.subtitle = element_text(size=16, hjust = 0.5),
                                                 axis.title=element_text(size=18), 
                                                 legend.text=element_text(size=14)+geom_text(size = 20)),
             nrow = 1)
rs2_outside_after_diag_plots_grob <- arrangeGrob(rs2_outside_after_diag_plots+ theme(axis.text.x = element_text(size=20),
                                                                                     axis.text.y = element_text(size=20),
                                                                                     plot.title = element_text(size=24, face = "bold", hjust = 0.5),
                                                                                     plot.subtitle = element_text(size=20, hjust = 0.5),
                                                                                     axis.title=element_text(size=20), 
                                                                                     legend.text=element_text(size=18)),
                                                 nrow = 1)
dev.off()
ggsave(rs2_outside_after_diag_plots_grob, file = file.path(out.dir2,"rs2_outside_after_influence_plot.pdf"), h = 14, w = 20)


#rs3 DHARMA ----
dharma_rs3_diagnostics <- (simulateResiduals(rs3))
dharma_rs3_diag_plots <- plot(dharma_rs3_diagnostics, return_list = TRUE)
dharma_rs3_diag_plots
grid.arrange(dharma_rs3_diag_plots+ theme(text = element_text(size = 20),
                                          axis.text.x = element_text(size=16),
                                          axis.text.y = element_text(size=16),
                                          plot.title = element_text(size=20, face = "bold", hjust = 0.5),
                                          plot.subtitle = element_text(size=16, hjust = 0.5),
                                          axis.title=element_text(size=18), 
                                          legend.text=element_text(size=14)+geom_text(size = 20)))
dev.off()




#rs4a_outside pp check----
rs4a_outside_diagnostics <- check_model(rs4a_outside, check = "pp_check")
rs4a_outside_diag_plots <- plot(rs4a_outside_diagnostics, return_list = TRUE)
grid.arrange(rs4a_outside_diag_plots+ theme(axis.text.x = element_text(size=16),
                                            axis.text.y = element_text(size=16),
                                            plot.title = element_text(size=20, face = "bold", hjust = 0.5),
                                            plot.subtitle = element_text(size=16, hjust = 0.5),
                                            axis.title=element_text(size=18), 
                                            legend.text=element_text(size=14)),
             nrow = 1)
rs4a_outside_diag_plots_grob <- arrangeGrob(rs4a_outside_diag_plots+ theme(axis.text.x = element_text(size=20),
                                                                           axis.text.y = element_text(size=20),
                                                                           plot.title = element_text(size=24, face = "bold", hjust = 0.5),
                                                                           plot.subtitle = element_text(size=20, hjust = 0.5),
                                                                           axis.title=element_text(size=20), 
                                                                           legend.text=element_text(size=18)),
                                            nrow = 1)
dev.off()
ggsave(rs4a_outside_diag_plots_grob, file = file.path(out.dir2,"rs4a_outside_pp_check.pdf"), h = 14, w = 20)

#rs4a_outside follow pp check----
rs4a_outside_follow_diagnostics <- check_model(rs4a_outside_follow, check = "pp_check")
rs4a_outside_follow_diag_plots <- plot(rs4a_outside_follow_diagnostics, return_list = TRUE)
grid.arrange(rs4a_outside_follow_diag_plots+ theme(axis.text.x = element_text(size=16),
                                                   axis.text.y = element_text(size=16),
                                                   plot.title = element_text(size=20, face = "bold", hjust = 0.5),
                                                   plot.subtitle = element_text(size=16, hjust = 0.5),
                                                   axis.title=element_text(size=18), 
                                                   legend.text=element_text(size=14)),
             nrow = 1)
rs4a_outside_follow_diag_plots_grob <- arrangeGrob(rs4a_outside_follow_diag_plots+ theme(axis.text.x = element_text(size=20),
                                                                                         axis.text.y = element_text(size=20),
                                                                                         plot.title = element_text(size=24, face = "bold", hjust = 0.5),
                                                                                         plot.subtitle = element_text(size=20, hjust = 0.5),
                                                                                         axis.title=element_text(size=20), 
                                                                                         legend.text=element_text(size=18)),
                                                   nrow = 1)
dev.off()
ggsave(rs4a_outside_follow_diag_plots_grob, file = file.path(out.dir2,"rs4a_outside_follow_pp_check.pdf"), h = 14, w = 20)

#rs4a_outside after pp check----

rs4a_outside_after_diagnostics <- check_model(rs4a_outside_after, check = "pp_check")
rs4a_outside_after_diag_plots <- plot(rs4a_outside_after_diagnostics, return_list = TRUE)
grid.arrange(rs4a_outside_after_diag_plots+ theme(axis.text.x = element_text(size=16),
                                                  axis.text.y = element_text(size=16),
                                                  plot.title = element_text(size=20, face = "bold", hjust = 0.5),
                                                  plot.subtitle = element_text(size=16, hjust = 0.5),
                                                  axis.title=element_text(size=18), 
                                                  legend.text=element_text(size=14)),
             nrow = 1)
rs4a_outside_after_diag_plots_grob <- arrangeGrob(rs4a_outside_after_diag_plots+ theme(axis.text.x = element_text(size=20),
                                                                                       axis.text.y = element_text(size=20),
                                                                                       plot.title = element_text(size=24, face = "bold", hjust = 0.5),
                                                                                       plot.subtitle = element_text(size=20, hjust = 0.5),
                                                                                       axis.title=element_text(size=20), 
                                                                                       legend.text=element_text(size=18)),
                                                  nrow = 1)
dev.off()
ggsave(rs4a_outside_after_diag_plots_grob, file = file.path(out.dir2,"rs4a_outside_after_pp_check.pdf"), h = 14, w = 20)








####
#rs6b pp check ----

rs6b_outside_diagnostics <- check_model(rs6b_outside, check = "pp_check")
rs6b_outside_diag_plots <- plot(rs6b_outside_diagnostics, return_list = TRUE)
grid.arrange(rs6b_outside_diag_plots+ theme(axis.text.x = element_text(size=16),
                                            axis.text.y = element_text(size=16),
                                            plot.title = element_text(size=20, face = "bold", hjust = 0.5),
                                            plot.subtitle = element_text(size=16, hjust = 0.5),
                                            axis.title=element_text(size=18), 
                                            legend.text=element_text(size=14)),
             nrow = 1)
rs6b_outside_diag_plots_grob <- arrangeGrob(rs6b_outside_diag_plots+ theme(axis.text.x = element_text(size=20),
                                                                           axis.text.y = element_text(size=20),
                                                                           plot.title = element_text(size=24, face = "bold", hjust = 0.5),
                                                                           plot.subtitle = element_text(size=20, hjust = 0.5),
                                                                           axis.title=element_text(size=20), 
                                                                           legend.text=element_text(size=18)),
                                            nrow = 1)
dev.off()
ggsave(rs6b_outside_diag_plots_grob, file = file.path(out.dir2,"rs6b_outside_pp_check.pdf"), h = 14, w = 20)






#rs6b follow pp check ----

rs6b_outside_follow_diagnostics <- check_model(rs6b_outside_follow, check = "pp_check")
rs6b_outside_follow_diag_plots <- plot(rs6b_outside_follow_diagnostics, return_list = TRUE)
grid.arrange(rs6b_outside_follow_diag_plots+ theme(axis.text.x = element_text(size=16),
                                                   axis.text.y = element_text(size=16),
                                                   plot.title = element_text(size=20, face = "bold", hjust = 0.5),
                                                   plot.subtitle = element_text(size=16, hjust = 0.5),
                                                   axis.title=element_text(size=18), 
                                                   legend.text=element_text(size=14)),
             nrow = 1)
rs6b_outside_follow_diag_plots_grob <- arrangeGrob(rs6b_outside_follow_diag_plots+ theme(axis.text.x = element_text(size=20),
                                                                                         axis.text.y = element_text(size=20),
                                                                                         plot.title = element_text(size=24, face = "bold", hjust = 0.5),
                                                                                         plot.subtitle = element_text(size=20, hjust = 0.5),
                                                                                         axis.title=element_text(size=20), 
                                                                                         legend.text=element_text(size=18)),
                                                   nrow = 1)
dev.off()
ggsave(rs6b_outside_follow_diag_plots_grob, file = file.path(out.dir2,"rs6b_outside_follow_pp_check.pdf"), h = 14, w = 20)



#rs6b post pp check ----

rs6b_outside_after_diagnostics <- check_model(rs6b_outside_after, check = "pp_check")
rs6b_outside_after_diag_plots <- plot(rs6b_outside_after_diagnostics, return_list = TRUE)
grid.arrange(rs6b_outside_after_diag_plots+ theme(axis.text.x = element_text(size=16),
                                                  axis.text.y = element_text(size=16),
                                                  plot.title = element_text(size=20, face = "bold", hjust = 0.5),
                                                  plot.subtitle = element_text(size=16, hjust = 0.5),
                                                  axis.title=element_text(size=18), 
                                                  legend.text=element_text(size=14)),
             nrow = 1)
rs6b_outside_after_diag_plots_grob <- arrangeGrob(rs6b_outside_after_diag_plots+ theme(axis.text.x = element_text(size=20),
                                                                                       axis.text.y = element_text(size=20),
                                                                                       plot.title = element_text(size=24, face = "bold", hjust = 0.5),
                                                                                       plot.subtitle = element_text(size=20, hjust = 0.5),
                                                                                       axis.title=element_text(size=20), 
                                                                                       legend.text=element_text(size=18)),
                                                  nrow = 1)
dev.off()
ggsave(rs6b_outside_after_diag_plots_grob, file = file.path(out.dir2,"rs6b_outside_after_pp_check.pdf"), h = 14, w = 20)



#rs3 pp check----
rs3_outside_diagnostics <- check_model(rs3_outside, check = "pp_check")
rs3_outside_diag_plots <- plot(rs3_outside_diagnostics, return_list = TRUE)
grid.arrange(rs3_outside_diag_plots+ theme(axis.text.x = element_text(size=16),
                                           axis.text.y = element_text(size=16),
                                           plot.title = element_text(size=20, face = "bold", hjust = 0.5),
                                           plot.subtitle = element_text(size=16, hjust = 0.5),
                                           axis.title=element_text(size=18), 
                                           legend.text=element_text(size=14)),
             nrow = 1)
rs3_outside_diag_plots_grob <- arrangeGrob(rs3_outside_diag_plots+ theme(axis.text.x = element_text(size=20),
                                                                         axis.text.y = element_text(size=20),
                                                                         plot.title = element_text(size=24, face = "bold", hjust = 0.5),
                                                                         plot.subtitle = element_text(size=20, hjust = 0.5),
                                                                         axis.title=element_text(size=20), 
                                                                         legend.text=element_text(size=18)),
                                           nrow = 1)
dev.off()
ggsave(rs3_outside_diag_plots_grob, file = file.path(out.dir2,"rs3_outside_pp_check.pdf"), h = 14, w = 20)

#rs3_outside after pp check----
rs3_outside_after_diagnostics <- check_model(rs3_outside_after, check = "pp_check")
rs3_outside_after_diag_plots <- plot(rs3_outside_after_diagnostics, return_list = TRUE)
grid.arrange(rs3_outside_after_diag_plots+ theme(axis.text.x = element_text(size=16),
                                                 axis.text.y = element_text(size=16),
                                                 plot.title = element_text(size=20, face = "bold", hjust = 0.5),
                                                 plot.subtitle = element_text(size=16, hjust = 0.5),
                                                 axis.title=element_text(size=18), 
                                                 legend.text=element_text(size=14)),
             nrow = 1)
rs3_outside_after_diag_plots_grob <- arrangeGrob(rs3_outside_after_diag_plots+ theme(axis.text.x = element_text(size=20),
                                                                                     axis.text.y = element_text(size=20),
                                                                                     plot.title = element_text(size=24, face = "bold", hjust = 0.5),
                                                                                     plot.subtitle = element_text(size=20, hjust = 0.5),
                                                                                     axis.title=element_text(size=20), 
                                                                                     legend.text=element_text(size=18)),
                                                 nrow = 1)
dev.off()
ggsave(rs3_outside_after_diag_plots_grob, file = file.path(out.dir2,"rs3_outside_after_pp_check.pdf"), h = 14, w = 20)

#rs3_outside after pp check----

rs3_outside_follow_diagnostics <- check_model(rs3_outside_follow, check = "pp_check")
rs3_outside_follow_diag_plots <- plot(rs3_outside_follow_diagnostics, return_list = TRUE)
grid.arrange(rs3_outside_follow_diag_plots+ theme(axis.text.x = element_text(size=16),
                                                  axis.text.y = element_text(size=16),
                                                  plot.title = element_text(size=20, face = "bold", hjust = 0.5),
                                                  plot.subtitle = element_text(size=16, hjust = 0.5),
                                                  axis.title=element_text(size=18), 
                                                  legend.text=element_text(size=14)),
             nrow = 1)
rs3_outside_follow_diag_plots_grob <- arrangeGrob(rs3_outside_follow_diag_plots+ theme(axis.text.x = element_text(size=20),
                                                                                       axis.text.y = element_text(size=20),
                                                                                       plot.title = element_text(size=24, face = "bold", hjust = 0.5),
                                                                                       plot.subtitle = element_text(size=20, hjust = 0.5),
                                                                                       axis.title=element_text(size=20), 
                                                                                       legend.text=element_text(size=18)),
                                                  nrow = 1)
dev.off()
ggsave(rs3_outside_follow_diag_plots_grob, file = file.path(out.dir2,"rs3_outside_follow_pp_check.pdf"), h = 14, w = 20)



colSums(influence.measures(rs1_outside)$is.inf)
colSums(influence.measures(rs1_outside_after)$is.inf)
colSums(influence.measures(rs1_outside_after)$is.inf)
influencePlot(rs1_outside_after)
influencePlot(rs2_outside_after)
check_model(rs1_outside_after)
pdf(paste0(file.path(out.dir2, "rs4a_after_pp_check.pdf")))
check_model(rs4a_after, check = "pp_check")
dev.off()


#AKAIKE----
cv_mods_outside[which(mod == "rs1_outside"), Akaike:= AICc_rs1_outside]
# cv_mods_outside[which(mod == "rs2_outside"), Akaike:= "N/A"]
cv_mods_outside[which(mod == "rs3_outside"), Akaike:= AICc_rs3]
cv_mods_outside[which(mod == "rs3a_outside"), Akaike:= AICc_rs3a_outside]
cv_mods_outside[which(mod == "rs4_outside"), Akaike:= AICc_rs4]
cv_mods_outside[which(mod == "rs4a_outside"), Akaike:= AICc_rs4a_outside]
cv_mods_outside[which(mod == "rs4c_outside"), Akaike:= AICc_rs4c_outside]
cv_mods_outside[which(mod == "rs4d_outside"), Akaike:= AICc_rs4d_outside]
cv_mods_outside[which(mod == "rs5_outside"), Akaike:= AICc_rs5]
cv_mods_outside[which(mod == "rs5a_outside"), Akaike:= AICc_rs5a_outside]
cv_mods_outside[which(mod == "rs6b_outside"), Akaike:= AICc_rs6b_outside]



cv_mods_outside_follow[which(mod == "rs1_outside_follow"), Akaike:= AICc_rs1_outside_follow]
cv_mods_outside_follow[which(mod == "rs3_outside_follow"), Akaike:= AICc_rs3_follow]
cv_mods_outside_follow[which(mod == "rs3a_outside_follow"), Akaike:= AICc_rs3a_outside_follow]
cv_mods_outside_follow[which(mod == "rs4_outside_follow"), Akaike:= AICc_rs4_follow]
cv_mods_outside_follow[which(mod == "rs4a_outside_follow"), Akaike:= AICc_rs4a_outside_follow]
cv_mods_outside_follow[which(mod == "rs4c_outside_follow"), Akaike:= AICc_rs4c_outside_follow]
cv_mods_outside_follow[which(mod == "rs4d_outside_follow"), Akaike:= AICc_rs4d_outside_follow]
cv_mods_outside_follow[which(mod == "rs5_outside_follow"), Akaike:= AICc_rs5_follow]
cv_mods_outside_follow[which(mod == "rs5a_outside_follow"), Akaike:= AICc_rs5a_outside_follow]
cv_mods_outside_follow[which(mod == "rs6b_outside_follow"), Akaike:= AICc_rs6b_outside_follow]

cv_mods_outside_after[which(mod == "rs1_outside_after"), Akaike:= AICc_rs1_outside_after]
cv_mods_outside_after[which(mod == "rs2_outside_after"), Akaike:= "N/A"]
cv_mods_outside_after[which(mod == "rs3_outside_after"), Akaike:= AICc_rs3_outside_after]
cv_mods_outside_after[which(mod == "rs3a_outside_after"), Akaike:= AICc_rs3a_outside_after]
cv_mods_outside_after[which(mod == "rs4_outside_after"), Akaike:= AICc_rs4_outside_after]
cv_mods_outside_after[which(mod == "rs4a_outside_after"), Akaike:= AICc_rs4a_outside_after]
cv_mods_outside_after[which(mod == "rs4c_outside_after"), Akaike:= AICc_rs4c_outside_after]
cv_mods_outside_after[which(mod == "rs4d_outside_after"), Akaike:= AICc_rs4d_outside_after]
cv_mods_outside_after[which(mod == "rs5_outside_after"), Akaike:= AICc_rs5_outside_after]
cv_mods_outside_after[which(mod == "rs5a_outside_after"), Akaike:= AICc_rs5a_outside_after]
cv_mods_outside_after[which(mod == "rs6b_outside_after"), Akaike:= AICc_rs6b_outside_after]


#BIC----
cv_mods_outside[which(mod == "rs1_outside"), BIC_val:= BIC_rs1_outside]
cv_mods_outside[which(mod == "rs2_outside"), BIC_val:= "N/A"]
cv_mods_outside[which(mod == "rs3_outside"), BIC_val:= BIC_rs3_outside]
cv_mods_outside[which(mod == "rs3a_outside"), BIC_val:= BIC_rs3a_outside]
cv_mods_outside[which(mod == "rs4_outside"), BIC_val:= BIC_rs4_outside]
cv_mods_outside[which(mod == "rs4a_outside"), BIC_val:= BIC_rs4a_outside]
cv_mods_outside[which(mod == "rs4a_outside"), BIC_val:= BIC_rs4a_outside]
cv_mods_outside[which(mod == "rs4c_outside"), BIC_val:= BIC_rs4c_outside]
cv_mods_outside[which(mod == "rs4d_outside"), BIC_val:= BIC_rs4d_outside]
cv_mods_outside[which(mod == "rs5_outside"), BIC_val:= BIC_rs5_outside]
cv_mods_outside[which(mod == "rs5a_outside"), BIC_val:= BIC_rs5a_outside]
cv_mods_outside[which(mod == "rs6b_outside"), BIC_val:= BIC_rs6b_outside]

cv_mods_outside_follow[which(mod == "rs1_outside_follow"), BIC_val:= BIC_rs1_outside_follow]
cv_mods_outside_follow[which(mod == "rs2_outside_follow"), BIC_val:= NA]
cv_mods_outside_follow[which(mod == "rs3_outside_follow"), BIC_val:= BIC_rs3_outside_follow]
cv_mods_outside_follow[which(mod == "rs3a_outside_follow"), BIC_val:= BIC_rs3a_outside_follow]
cv_mods_outside_follow[which(mod == "rs4_outside_follow"), BIC_val:= BIC_rs4_outside_follow]
cv_mods_outside_follow[which(mod == "rs4a_outside_follow"), BIC_val:= BIC_rs4a_outside_follow]
cv_mods_outside_follow[which(mod == "rs4c_outside_follow"), BIC_val:= BIC_rs4c_outside_follow]
cv_mods_outside_follow[which(mod == "rs4d_outside_follow"), BIC_val:= BIC_rs4d_outside_follow]
cv_mods_outside_follow[which(mod == "rs5_outside_follow"), BIC_val:= BIC_rs5_outside_follow]
cv_mods_outside_follow[which(mod == "rs5a_outside_follow"), BIC_val:= BIC_rs5a_outside_follow]
cv_mods_outside_follow[which(mod == "rs6b_outside_follow"), BIC_val:= BIC_rs6b_outside_follow]

cv_mods_outside_after[which(mod == "rs1_outside_after"), BIC_val:= BIC_rs1_outside_after]
cv_mods_outside_after[which(mod == "rs2_outside_after"), BIC_val:= NA]
cv_mods_outside_after[which(mod == "rs3_outside_after"), BIC_val:= BIC_rs3_outside_after]
cv_mods_outside_after[which(mod == "rs3a_outside_after"), BIC_val:= BIC_rs3a_outside_after]
cv_mods_outside_after[which(mod == "rs4_outside_after"), BIC_val:= BIC_rs4_outside_after]
cv_mods_outside_after[which(mod == "rs4a_outside_after"), BIC_val:= BIC_rs4a_outside_after]
cv_mods_outside_after[which(mod == "rs4c_outside_after"), BIC_val:= BIC_rs4c_outside_after]
cv_mods_outside_after[which(mod == "rs4d_outside_after"), BIC_val:= BIC_rs4d_outside_after]
cv_mods_outside_after[which(mod == "rs5_outside_after"), BIC_val:= BIC_rs5_outside_after]
cv_mods_outside_after[which(mod == "rs5a_outside_after"), BIC_val:= BIC_rs5a_outside_after]
cv_mods_outside_after[which(mod == "rs6b_outside_after"), BIC_val:= BIC_rs6b_outside_after]





#Bayesian Analysis ----
rs_pois_outside_follow <- stan_glm(Incidence~Treatment+A+B+C+D+E+F+G+H+I+log(Hist3yr)+log(Baseline),
                                   data = rbctconf_outside_follow, family = "poisson",
                                   offset = log(Baseline),
                                   prior = normal(0, 1),
                                   prior_intercept = normal(0, 2),
                                   diagnostic_file = file.path(tempdir(), "df.csv"),refresh=0)
loo_rs_pois_outside_follow <- loo(rs_pois_outside_follow, k_threshold = 0.7)
loo_rs_pois_outside_follow
pp_check(rs_pois_outside_follow, plotfun = "dens_overlay")



length(which(exp(draws_rs_pois_outside_follow$TreatmentProactive) > 1))/nrow(draws_rs_pois_outside_follow)

draws_rs_pois_outside_follow <- as.data.table(rs_pois_outside_follow)

loo_rs_pois_outside_follow <- loo(rs_pois_outside_follow, k_threshold = 0.7)
loo_rs_pois_outside_follow





set.seed(150499)
rs_pois_outside <- stan_glm(Incidence~Treatment+A+B+C+D+E+F+G+H+I+log(Hist3yr)+log(Baseline),
                            data = rbctconf_outside, family = "poisson",
                            offset = log(Baseline),
                            prior = normal(0, 1),
                            prior_intercept = normal(0, 2),
                            diagnostic_file = file.path(tempdir(), "df.csv"),refresh=0)

loo_rs_pois_outside <- loo(rs_pois_outside, k_threshold = 0.7)
pp_check(rs_pois_outside, plotfun = "dens_overlay")





summary(rs_pois_outside, probs=c(0.025,0.5,0.975), digits=3)
summary(rs_pois_outside, probs=c(0.05,0.5,0.95), digits=3)
draws_rs_pois_outside <- as.data.table(rs_pois_outside)
loo_rs_pois_outside <- loo(rs_pois_outside, k_threshold = 0.7)
rs_pois_outside_posterior_predictive <- posterior_predict(rs_pois_outside)
ppc_stat(rbctconf_outside$Incidence, rs_pois_outside_posterior_predictive, stat = "max")



rs_pois_outside_after <- stan_glm(Incidence~Treatment+A+B+C+D+E+F+G+H+I+log(Hist3yr)+log(Baseline),
                                  data = rbctconf_outside_after, family = "poisson",
                                  offset = log(Baseline),
                                  prior = normal(0, 1),
                                  prior_intercept = normal(0, 1),
                                  diagnostic_file = file.path(tempdir(), "df.csv"),refresh=0)

loo_rs_pois_outside_after <- loo(rs_pois_outside_after, k_threshold = 0.7)
pp_check(rs_pois_outside_after, plotfun = "dens_overlay")
draws_rs_pois_outside_after <- as.data.table(rs_pois_outside_after)


###Dispersion: Prior vs Posterior Density----
rs_prior_dispersion_dt <- data.table()
rs_outside_dispersion_prior_posterior <- copy(rs_prior_dispersion_dt)
draws_rs <- as.data.table(rs_outside)
rs_po_dispersion_trmnt <- subset(draws_rs, select = "reciprocal_dispersion")
setnames(rs_po_dispersion_trmnt, "reciprocal_dispersion", "posterior_val")
rs_outside_dispersion_prior_posterior[, posterior_val:= rs_po_dispersion_trmnt$posterior_val]

rs_outside_dispersion_densities <- ggplot(rs_outside_dispersion_prior_posterior)+
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
  xlim(c(0,750))
rs_outside_dispersion_densities
ggsave(rs_outside_dispersion_densities, file = file.path(out.dir2, "rs_outside_dispersion_densities.pdf"), h = 16, w = 18)

rs_outside_improved_prior_dispersion_dt <- data.table()
rs_outside_improved_dispersion_prior_posterior <- copy(rs_outside_improved_prior_dispersion_dt)
draws_rs_outside_improved <- as.data.table(rs_outside_improved)
rs_outside_improved_po_dispersion_trmnt <- subset(draws_rs_outside_improved, select = "reciprocal_dispersion")
setnames(rs_outside_improved_po_dispersion_trmnt, "reciprocal_dispersion", "posterior_val")
rs_outside_improved_dispersion_prior_posterior[, posterior_val:= rs_outside_improved_po_dispersion_trmnt$posterior_val]

rs_outside_improved_dispersion_densities <- ggplot(rs_outside_improved_dispersion_prior_posterior)+
  stat_function(aes(col = "Prior"), fun = dhcauchy, n = 101, args = list(sigma = 5))+
  
  geom_density(aes(posterior_val, col = "Posterior"))+
  theme_bw()+
  labs(x = "Value", y = "Density",
       title = "Improved Model rs: Posterior vs Prior of Reciprocal Dispersion Parameter",
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
  xlim(c(0,750))
ggsave(rs_outside_improved_dispersion_densities, file = file.path(out.dir2, "rs_outside_improved_dispersion_densities.pdf"), h = 16, w = 18)
####










#rs1_outside_aB (Improved) Prior vs Posterior----
rs1_outside_B_improved_pois_prior_beta_trmnts_dt <- data.table()
rs1_outside_B_improved_prior_posterior <- copy(rs1_outside_B_improved_pois_prior_beta_trmnts_dt)
draws_rs1_outside_B_improved <- as.data.table(rs1_outside_B_improved)
rs1_outside_B_improved_po_beta_trmnt <- subset(draws_rs1_outside_B_improved, select = "TreatmentProactive")
setnames(rs1_outside_B_improved_po_beta_trmnt, "TreatmentProactive", "posterior_val")
rs1_outside_B_improved_prior_posterior[, posterior_val:= rs1_outside_B_improved_po_beta_trmnt$posterior_val]



rs1_outside_B_improved_densities <- ggplot(rs1_outside_B_improved_prior_posterior)+
  stat_function(aes(col = "Prior"), fun = dnorm, n = 101, size = 1.5, args = list(mean = 0, sd = 1))+
  geom_density(aes(posterior_val, col = "Posterior"), size = 2)+
  theme_bw()+
  scale_x_continuous(limits = c(-2.5, 2.5))+
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
rs1_outside_B_improved_densities
ggsave(rs1_outside_B_improved_densities, file = file.path(out.dir2, "rs1_outside_B_improved_densities.pdf"), h = 16, w = 18)







update_geom_defaults("line", list(size = 2,
                                  width = 2))

update_geom_defaults("pointrange", list(size = 1,
                                        width = 1,
                                        fatten = 1,
                                        linewidth = 1.5))
rs_outside_post_prior <- posterior_vs_prior(rs_outside,
                                            pars = "TreatmentProactive")+ylim(c(-10, 10))+
  labs(title = "Model a.1: Proactive Treatment Parameter")+theme_bw()+
  theme(text = element_text(size = 26),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=30),
        plot.title = element_text(size=30, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size=30, hjust = 0.5),
        axis.title=element_text(size=28), 
        legend.position = "none")+
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank())
rs_outside_post_prior
rs_outside_post_prior$layers[[1]]$aes_params$colour <- c("#2196F3", "red")
rs_outside_post_prior
ggsave(rs_outside_post_prior, file = file.path(out.dir2, "rs_outside_post_prior.pdf"), h = 20, w = 16)



rs_outside_improved_post_prior <- posterior_vs_prior(rs_outside_improved,
                                                     pars = "TreatmentProactive")+ylim(c(-10, 10))+
  labs(title = "Model a.2: Proactive Treatment Parameter")+theme_bw()+
  theme(text = element_text(size = 26),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=30),
        plot.title = element_text(size=30, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size=30, hjust = 0.5),
        axis.title=element_text(size=28), 
        legend.position = "none")+
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank())
rs_outside_improved_post_prior
rs_outside_improved_post_prior$layers[[1]]$aes_params$colour <- c("#2196F3", "red")
rs_outside_improved_post_prior



ggsave(rs_outside_improved_post_prior, file = file.path(out.dir2, "rs_outside_improved_post_prior.pdf"), h = 20, w = 16)

rs_outside_post_prior_grob
grid.arrange(rs_outside_post_prior, rs_outside_improved_post_prior, nrow = 1)

rs_outside_post_prior_grob <- arrangeGrob(rs_outside_post_prior, rs_outside_improved_post_prior, nrow = 1)
ggsave(rs_outside_post_prior_grob, file = file.path(out.dir2, "rs_outside_post_prior_grob.pdf"), h = 12, w = 22)



#rs Dispersion---
rs_outside_post_prior_dispersion <- 
  posterior_vs_prior(rs_outside,
                                                       pars = "reciprocal_dispersion")+ylim(c(0, 1100))+
  labs(title = "Model a.1: Reciprocal Dispersion Parameter")+theme_bw()+
  theme(text = element_text(size = 26),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=30),
        plot.title = element_text(size=30, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size=30, hjust = 0.5),
        axis.title=element_text(size=28), 
        legend.position = "none")+
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank())
rs_outside_post_prior_dispersion$layers[[1]]$aes_params$colour <- c("#2196F3", "red")
rs_outside_post_prior_dispersion
rs_outside_post_prior_dispersion
rs_outside_post_prior_dispersion
ggsave(rs_outside_post_prior_dispersion, file = file.path(out.dir2, "rs_outside_post_prior_dispersion.pdf"), h = 20, w = 16)

rs_outside_improved_post_prior_dispersion <- posterior_vs_prior(rs_outside_improved,
                                                                pars = "reciprocal_dispersion")+ylim(c(0, 1100))+
  labs(title = "Model a.2: Reciprocal Dispersion Parameter")+theme_bw()+
  theme(text = element_text(size = 26),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=30),
        plot.title = element_text(size=30, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size=30, hjust = 0.5),
        axis.title=element_text(size=28), 
        legend.position = "none")+
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank())
rs_outside_improved_post_prior_dispersion$layers[[1]]$aes_params$colour <- c("#2196F3", "red")
rs_outside_improved_post_prior_dispersion

rs_outside_post_prior_dispersion_grob <- arrangeGrob(rs_outside_post_prior_dispersion, rs_outside_improved_post_prior_dispersion, nrow = 1)
ggsave(rs_outside_post_prior_dispersion_grob, file = file.path(out.dir2, "rs_outside_post_prior_dispersion_grob.pdf"), h = 12, w = 22)




rs_pois_outside_after <- stan_glm(Incidence~Treatment+A+B+C+D+E+F+G+H+I+log(Hist3yr)+log(Baseline),
                                  data = rbctconf_outside_after, family = "poisson",
                                  offset = log(Baseline),
                                  prior = normal(0, 0.5, autoscale = TRUE),
                                  prior_intercept = normal(0, 0.5 ,autoscale = TRUE),
                                  diagnostic_file = file.path(tempdir(), "df.csv"),refresh=0)
loo(rs_pois_outside_after, k_threshold = 0.7)
loo(rs_pois_outside_after)
rs_pois_outside_after$prior.info



draws_rs_pois_outside_after <- as.data.table(rs_pois_outside_after)


#Bayesian Model 1 of preprint----
rs_outside <- stan_glm.nb(Incidence~Treatment+A+B+C+D+E+F+G+H+I+log(Hist3yr)+log(Baseline),
                          data = rbctconf_outside, prior_intercept=normal(0, 10),
                          diagnostic_file = file.path(tempdir(), "df.csv"),refresh=0)
rs_outside_posterior_predictive <- posterior_predict(rs_outside)
draws_rs_outside <- as.data.table(rs_outside)
loo_rs_outside <- loo(rs_outside)

rs_outside_improved <- stan_glm.nb(Incidence~Treatment+A+B+C+D+E+F+G+H+I+log(Hist3yr)+log(Baseline),
                                   data = rbctconf_outside, 
                                   prior_aux= cauchy(0, 5, autoscale = TRUE),
                                   prior = normal(0, 1, autoscale = TRUE),
                                   prior_intercept = normal(0, 1, autoscale = TRUE),
                                   diagnostic_file = file.path(tempdir(), "df.csv"),refresh=0)
loo_rs_outside_improved <- loo(rs_outside_improved, k_threshold = 0.7)

#Small values for reciprocal_dispersion correspond to greater dispersion
draws_rs_outside_improved <- as.data.table(rs_outside_improved)

rs_outside_follow <-stan_glm.nb(Incidence~Treatment+A+B+C+D+E+F+G+H+I+log(Hist3yr)+log(Baseline),prior_intercept=normal(0, 10),
                                prior = normal(0, 10),
                                diagnostic_file = file.path(tempdir(), "df.csv"),refresh=0, data = rbctconf_outside_follow)
pp_check(rs_outside_follow, plotfun = "dens_overlay")
loo_rs_outside_follow <- loo(rs_outside_follow, k_threshold = 0.7)
draws_rs_outside_follow <- as.data.table(rs_outside_follow)

rs_outside_follow_improved <- stan_glm.nb(Incidence~Treatment+A+B+C+D+E+F+G+H+I+log(Hist3yr)+log(Baseline),
                                          data = rbctconf_outside_follow,
                                          prior_aux= cauchy(0, 5),
                                          prior = normal(0, 1),
                                          diagnostic_file = file.path(tempdir(), "df.csv"),refresh=0)
loo_rs_outside_follow_improved <- loo(rs_outside_follow_improved, k_threshold = 0.7)
pp_check(rs_outside_follow_improved)
draws_rs_outside_follow_improved <- as.data.table(rs_outside_follow_improved)

rs_outside_after <-stan_glm.nb(Incidence~Treatment+A+B+C+D+E+F+G+H+I+log(Hist3yr)+log(Baseline),prior_intercept=normal(0, 10),
                               prior = normal(0, 10),
                               diagnostic_file = file.path(tempdir(), "df.csv"),refresh=0, data = rbctconf_outside_after)
pp_check(rs_outside_after, plotfun = "dens_overlay")
draws_rs_outside_after <- as.data.table(rs_outside_after)
loo_rs_outside_after <- loo(rs_outside_after, k_threshold = 0.7)

rs_outside_after_improved <- stan_glm.nb(Incidence~Treatment+A+B+C+D+E+F+G+H+I+log(Hist3yr)+log(Baseline),
                                         data = rbctconf_outside_after, 
                                         prior_aux= cauchy(0, 5),
                                         prior = normal(0, 1),
                                         diagnostic_file = file.path(tempdir(), "df.csv"),refresh=0)
loo_rs_outside_after_improved <- loo(rs_outside_after_improved, k_threshold = 0.7)
pp_check(rs_outside_after_improved)
draws_rs_outside_after_improved <- as.data.table(rs_outside_after_improved)






#Bayesian Model 2----
rs1_outside_B<-stan_glm.nb(Incidence~Treatment+log(hdyrsrisk)+log(Hist3yr),
                           prior_intercept=normal(0,10), prior=normal(0,10),
                           data=rbctconf_outside,diagnostic_file =file.path(tempdir(),"df.csv"),refresh=0)

pp_check(rs1_outside_B, plotfun = "dens_overlay")
loo_rs1_outside_B <- loo(rs1_outside_B)
draws_rs1_outside_B <- as.data.table(rs1_outside_B)

rs1_outside_B_improved<-stan_glm.nb(Incidence~Treatment+log(hdyrsrisk)+log(Hist3yr),
                                    prior=normal(0,1),
                                    prior_aux = cauchy(0, 5),
                                    data=rbctconf_outside,diagnostic_file =file.path(tempdir(),"df.csv"),refresh=0)
loo_rs1_outside_B_improved <- loo(rs1_outside_B_improved, k_threshold = 0.7)
pp_check(rs1_outside_B_improved, plotfun = "dens_overlay")
draws_rs1_outside_B_improved <- as.data.table(rs1_outside_B_improved)

rs1_outside_B_follow<-stan_glm.nb(Incidence~Treatment+log(hdyrsrisk)+log(Hist3yr),
                                  prior_intercept=normal(0,10), prior=normal(0,10),
                                  data=rbctconf_outside_follow,diagnostic_file =file.path(tempdir(),"df.csv"),refresh=0)
loo_rs1_outside_B_follow <- loo(rs1_outside_B_follow)
pp_check(rs1_outside_B_follow, plotfun = "dens_overlay")
draws_rs1_outside_B_follow <- as.data.table(rs1_outside_B_follow)

rs1_outside_B_follow_improved<-stan_glm.nb(Incidence~Treatment+log(hdyrsrisk)+log(Hist3yr),
                                           prior_intercept = normal(0, 10, autoscale = TRUE),
                                           prior = normal(0, 1, autoscale = TRUE),
                                           prior_aux= cauchy(0, 5, autoscale = TRUE),
                                           data=rbctconf_outside_follow,diagnostic_file =file.path(tempdir(),"df.csv"),refresh=0)
loo_rs1_outside_B_follow_improved <- loo(rs1_outside_B_follow_improved, k_threshold = 0.7)
draws_rs1_outside_B_follow_improved <- as.data.table(rs1_outside_B_follow_improved)


rs1_outside_B_after<-stan_glm.nb(Incidence~Treatment+log(hdyrsrisk)+log(Hist3yr),
                                 prior_intercept=normal(0,10), prior=normal(0,10),
                                 data=rbctconf_outside_after,diagnostic_file =file.path(tempdir(),"df.csv"),refresh=0)
loo_rs1_outside_B_after <- loo(rs1_outside_B_after)
pp_check(rs1_outside_B_after, plotfun = "dens_overlay")
draws_rs1_outside_B_after <- as.data.table(rs1_outside_B_after)

rs1_outside_B_after_improved<-stan_glm.nb(Incidence~Treatment+log(hdyrsrisk)+log(Hist3yr),
                                          prior_intercept = normal(0, 5, autoscale = TRUE),
                                          prior = normal(0, 1, autoscale = TRUE),
                                          prior_aux= cauchy(0, 5, autoscale = TRUE),
                                          data=rbctconf_outside_after,diagnostic_file =file.path(tempdir(),"df.csv"),refresh=0)
draws_rs1_outside_B_after_improved <- as.data.table(rs1_outside_B_after_improved)
pp_check(rs1_outside_B_after_improved, plotfun = "dens_overlay")
loo_rs1_outside_B_after_improved <- loo(rs1_outside_B_after_improved, k_threshold = 0.7)



#Model rs1_outside_aB----
#Uses herd years as offset
rs1_outside_aB<-stan_glm.nb(Incidence~Treatment+log(hdyrsrisk)+log(Hist3yr),
                            offset = log(hdyrsrisk), 
                            prior_intercept=normal(0,10), prior=normal(0,10),
                            data = rbctconf_outside, refresh=0)
pp_check(rs1_outside_aB, plotfun = "dens_overlay")
loo_rs1_outside_aB <- loo(rs1_outside_aB)
draws_rs1_outside_aB <- as.data.table(rs1_outside_aB)

rs1_outside_aB_improved<-stan_glm.nb(Incidence~Treatment+log(hdyrsrisk)+log(Hist3yr),
                                     prior_intercept = normal(0, 5),
                                     prior = normal(0, 1),
                                     prior_aux= cauchy(0, 5),
                                     offset = log(hdyrsrisk),
                                     data = rbctconf_outside, refresh=0)
loo_rs1_outside_aB_improved <- loo(rs1_outside_aB_improved, k_threshold = 0.7)
draws_rs1_outside_aB_improved <- as.data.table(rs1_outside_aB_improved)


pp_check(rs1_outside_aB_improved, plotfun = "dens_overlay")
pp_check(rs1_outside_aB, plotfun = "dens_overlay")

grid.arrange(rs1_outside_aB_improved_post_plot, rs1_outside_aB_after_improved_post_plot, ncol = 1)
rs1_outside_aB_follow<-stan_glm.nb(Incidence~Treatment+log(hdyrsrisk)+log(Hist3yr),
                                   offset = log(hdyrsrisk),
                                   prior_intercept=normal(0,10), prior=normal(0,10),
                                   data = rbctconf_outside_follow, refresh=0)
pp_check(rs1_outside_aB_follow, plotfun = "dens_overlay")

loo_rs1_outside_aB_follow <- loo(rs1_outside_aB_follow)
draws_rs1_outside_aB_follow <- as.data.table(rs1_outside_aB_follow)


rs1_outside_aB_after<-stan_glm.nb(Incidence~Treatment+log(hdyrsrisk)+log(Hist3yr),
                                   offset = log(hdyrsrisk),
                                   prior_intercept=normal(0,10), prior=normal(0,10),
                                   data = rbctconf_outside_after, refresh=0)
pp_check(rs1_outside_aB_after, plotfun = "dens_overlay")

loo_rs1_outside_aB_after <- loo(rs1_outside_aB_after, k_threshold = 0.7)
draws_rs1_outside_aB_after <- as.data.table(rs1_outside_aB_after)

rs1_outside_aB_improved_post_plot <- mcmc_areas(draws_rs1_outside_aB_improved, pars = "TreatmentProactive", prob = 0.95) + xlab('Log Effect')+
  ggplot2::theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13))+
  theme(plot.subtitle = element_text(hjust = 0.5, size = 12))+
  labs(title = "Posterior Distribution of Treatment Parameter", subtitle = "with medians and 95% credible intervals")
rs1_outside_aB_improved_post_plot


rs1_outside_aB_follow_improved<-stan_glm.nb(Incidence~Treatment+log(hdyrsrisk)+log(Hist3yr),
                                            prior_intercept = normal(0, 5),
                                            prior = normal(0, 1),
                                            prior_aux= cauchy(0, 5),
                                            offset = log(hdyrsrisk),
                                            data = rbctconf_outside_follow, refresh=0)
pp_check(rs1_outside_aB_follow_improved, plotfun = "dens_overlay")
loo_rs1_outside_aB_follow_improved <- loo(rs1_outside_aB_follow_improved)
draws_rs1_outside_aB_follow_improved <- as.data.table(rs1_outside_aB_follow_improved)
rs1_outside_aB_follow_improved_post_plot <- mcmc_areas(draws_rs1_outside_aB_follow_improved, pars = "TreatmentProactive", prob = 0.95) + xlab('Log Effect')+
  ggplot2::theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20))+
  theme(text = element_text(size = 22), 
        plot.title = element_text(hjust = 0.5, face = "bold", size = 24))+
  theme(plot.subtitle = element_text(hjust = 0.5, size = 18),
        axis.text.x = element_text(size = 28),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 28)
  )
rs1_outside_aB_follow_improved_post_plot
ggsave(rs1_outside_aB_follow_improved_post_plot, file = file.path(out.dir2, "rs1_outside_aB_follow_improved_post_plot.pdf"), h = 8, w = 16)








rs1_outside_aB_follow_improved_post_plot <- mcmc_areas(draws_rs1_outside_aB_follow_improved, pars = "TreatmentProactive", prob = 0.95) + xlab('Log Effect')+
  ggplot2::theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))+
  theme(plot.subtitle = element_text(hjust = 0.5, size = 15),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        text = element_text(size = 15))+
  labs(title = "Posterior Distribution of Treatment Parameter", subtitle = "with medians and 95% credible intervals")
rs1_outside_aB_follow_improved_post_plot
ggsave(rs1_outside_aB_follow_improved_post_plot, file = file.path(out.dir2, "rs1_outside_aB_follow_improved_post_plot.pdf"), h = 8, w = 16)







rs1_outside_aB_follow<-stan_glm.nb(Incidence~Treatment+log(hdyrsrisk)+log(Hist3yr),
                                   offset = log(hdyrsrisk),
                                   prior_intercept=normal(0,10), prior=normal(0,10),
                                   data = rbctconf_outside_follow, refresh=0)
loo_rs1_outside_aB_follow <- loo(rs1_outside_aB_follow)
pp_check(rs1_outside_aB_follow)
draws_rs1_outside_aB_follow <- as.data.table(rs1_outside_aB_follow)

rs1_outside_aB_improved_post_plot <- mcmc_areas(draws_rs1_outside_aB_improved, pars = "TreatmentProactive", prob = 0.95) + xlab('Log Effect')+
  ggplot2::theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13))+
  theme(plot.subtitle = element_text(hjust = 0.5, size = 12))+
  labs(title = "Posterior Distribution of Treatment Parameter", subtitle = "with medians and 95% credible intervals")
rs1_outside_aB_improved_post_plot


rs1_outside_aB_after_improved<-stan_glm.nb(Incidence~Treatment+log(hdyrsrisk)+log(Hist3yr),
                                           prior_intercept = normal(0, 5),
                                           prior = normal(0, 1),
                                           prior_aux= cauchy(0, 5),
                                           offset = log(hdyrsrisk),
                                           data = rbctconf_outside_after, refresh=0)
pp_check(rs1_outside_aB_after_improved, plotfun = "dens_overlay")
loo_rs1_outside_aB_after_improved <- loo(rs1_outside_aB_after_improved, k_threshold = 0.7)
loo_rs1_outside_aB_after_improved
draws_rs1_outside_aB_after_improved <- as.data.table(rs1_outside_aB_after_improved)

rs1_outside_aB_after_improved_post_plot <- mcmc_areas(draws_rs1_outside_aB_after_improved, pars = "TreatmentProactive", prob = 0.95) + xlab('Log Effect')+
  ggplot2::theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20))+
  theme(text = element_text(size = 22), 
        plot.title = element_text(hjust = 0.5, face = "bold", size = 24))+
  theme(plot.subtitle = element_text(hjust = 0.5, size = 18),
        axis.text.x = element_text(size = 28),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 28)
  )
rs1_outside_aB_after_improved_post_plot
ggsave(rs1_outside_aB_after_improved_post_plot, file = file.path(out.dir2, "rs1_outside_aB_after_improved_post_plot.pdf"), h = 8, w = 16)


rs1_outside_aB_after_improved_post_plot <- mcmc_areas(draws_rs1_outside_aB_after_improved, pars = "TreatmentProactive", prob = 0.95) + xlab('Log Effect')+
  ggplot2::theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20))+
  theme(text = element_text(size = 22), 
        plot.title = element_text(hjust = 0.5, face = "bold", size = 24))+
  theme(plot.subtitle = element_text(hjust = 0.5, size = 18),
        axis.text.x = element_text(size = 28),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 28)
  )
rs1_outside_aB_after_improved_post_plot
ggsave(rs1_outside_aB_after_improved_post_plot, file = file.path(out.dir2, "rs1_outside_aB_after_improved_post_plot.pdf"), h = 8, w = 16)



rs1_outside_aB_after_improved_post_plot <- mcmc_areas(draws_rs1_outside_aB_after_improved, pars = "TreatmentProactive", prob = 0.95) + xlab('Log Effect')+
  ggplot2::theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))+
  theme(plot.subtitle = element_text(hjust = 0.5, size = 15),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        text = element_text(size = 15))+
  labs(title = "Posterior Distribution of Treatment Parameter", subtitle = "with medians and 95% credible intervals")
rs1_outside_aB_after_improved_post_plot
ggsave(rs1_outside_aB_after_improved_post_plot, file = file.path(out.dir2, "rs1_outside_aB_after_improved_post_plot.pdf"), h = 8, w = 16)


#rs2B ----
rs2_outsideB<-stan_glm.nb(Incidence~log(hdyrsrisk)+log(Hist3yr),
                          prior_intercept=normal(0,10), prior=normal(0,10),
                          data=rbctconf_outside,diagnostic_file = file.path(tempdir(), "df.csv"),refresh=0)
loo_rs2_outsideB <- loo(rs2_outsideB, k_threshold = 0.7)
loo_rs2_outsideB

pp_check(rs2_outsideB, plotfun = "dens_overlay")

rs2_outsideB_improved<-stan_glm.nb(Incidence~log(hdyrsrisk)+log(Hist3yr),
                                   prior=normal(0, 1),
                                   prior_aux = cauchy(0, 5),
                                   data=rbctconf_outside,diagnostic_file = file.path(tempdir(), "df.csv"),refresh=0)
loo_rs2_outsideB_improved <- loo(rs2_outsideB_improved, k_threshold = 0.7)
loo_rs2_outsideB_improved
pp_check(rs2_outsideB_improved)
pp_check(rs2_outsideB_improved, plotfun = "dens_overlay")
pp_check(rs1_outside_B, plotfun = "dens_overlay")


rs2_outsideB_follow<-stan_glm.nb(Incidence~log(hdyrsrisk)+log(Hist3yr),
                                 prior_intercept=normal(0,10), prior=normal(0,10),
                                 data=rbctconf_outside_follow,diagnostic_file = file.path(tempdir(), "df.csv"),refresh=0)
loo_rs2_outsideB_follow <- loo(rs2_outsideB_follow)
pp_check(rs2_outsideB_follow, plotfun = "dens_overlay")



rs2_outsideB_follow_improved<-stan_glm.nb(Incidence~log(hdyrsrisk)+log(Hist3yr),
                                          prior=normal(0,1),
                                          prior_aux = cauchy(0, 5),
                                          data=rbctconf_outside_follow,diagnostic_file = file.path(tempdir(), "df.csv"),refresh=0)
loo_rs2_outsideB_follow_improved <- loo(rs2_outsideB_follow_improved, k_threshold = 0.7)
loo_rs2_outsideB_follow_improved

pp_check(rs2_outsideB_follow_improved, plotfun = "dens_overlay")

rs2_outsideB_after<-stan_glm.nb(Incidence~log(hdyrsrisk)+log(Hist3yr),
                                prior_intercept=normal(0,10), prior=normal(0,10),
                                data=rbctconf_outside_after,diagnostic_file = file.path(tempdir(), "df.csv"),refresh=0)
loo_rs2_outsideB_after <- loo(rs2_outsideB_after)
loo_rs2_outsideB_after
pp_check(rs2_outsideB_after, plotfun = "dens_overlay")



rs2_outsideB_after_improved<-stan_glm.nb(Incidence~log(hdyrsrisk)+log(Hist3yr),
                                         prior_aux = cauchy(0, 5),
                                         prior = normal(0, 1),
                                         data=rbctconf_outside_after,diagnostic_file = file.path(tempdir(), "df.csv"),refresh=0)
loo_rs2_outsideB_after_improved <- loo(rs2_outsideB_after_improved, k_threshold = 0.7)
loo_rs2_outsideB_after_improved
pp_check(rs2_outsideB_after_improved, plotfun = "dens_overlay")




#Bayesian Plotting -----
rs2_outsideB_pp <- pp_check(rs2_outsideB, plotfun = "dens_overlay")+ 
  theme(text = element_text(size = 20),axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        plot.title = element_text(size=22, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size=20, hjust = 0.5),
        axis.title=element_text(size=20), 
        legend.text=element_text(size=18))+
  theme(legend.position = "bottom")+
  labs(title = "Posterior Predictive Check of Model rs2_outsideB", x = "Incidence", y = "Density")
rs2_outsideB_pp

rs1_outside_aB_pp <- pp_check(rs1_outside_aB, plotfun = "dens_overlay")+ 
  theme(text = element_text(size = 20),axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        plot.title = element_text(size=22, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size=20, hjust = 0.5),
        axis.title=element_text(size=20), 
        legend.text=element_text(size=18))+
  theme(legend.position = "bottom")+
  labs(title = "Posterior Predictive Check of Model rs1_outside_aB", x = "Incidence", y = "Density")
rs1_outside_aB_pp
ggsave(rs1_outside_aB_pp, file = file.path(out.dir2, "rs1_outside_aB_pp.pdf"), h = 14, w = 8)

rs1_outside_aB_improved_pp <- pp_check(rs1_outside_aB_improved, plotfun = "dens_overlay")+ 
  theme(text = element_text(size = 20),axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        plot.title = element_text(size=22, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size=20, hjust = 0.5),
        axis.title=element_text(size=20), 
        legend.text=element_text(size=18))+
  theme(legend.position = "bottom")+
  labs(title = "Posterior Predictive Check of Improved Model rs1_outside_aB", x = "Incidence", y = "Density")

ggsave(rs1_outside_aB_improved_pp, file = file.path(out.dir2, "rs1_outside_aB_improved_pp.pdf"), h = 14, w = 8)
rs1_outside_aB_improved_posterior_predictive <- posterior_predict(rs1_outside_aB_improved)
ppc_max_rs1_outside_aB <- ppc_stat(rbctconf_outside$Incidence, rs1_outside_aB_posterior_predictive, stat = "max")+
  xlim(c(0, 400))+ 
  theme(text = element_text(size = 20),axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        plot.title = element_text(size=22, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size=20, hjust = 0.5),
        axis.title=element_text(size=20), 
        legend.text=element_text(size=18))+
  theme(legend.position = "bottom")+
  labs(title = "Model rs1_outside_aB Posterior Predictive Check: Test Statistic (Maximum)", x = "Maximum Incidence", y = "Count")

ppc_max_rs1_outside_aB_improved <-ppc_stat(rbctconf_outside$Incidence, rs1_outside_aB_improved_posterior_predictive, stat = "max")+
  xlim(c(0, 400))+ 
  theme(text = element_text(size = 20),axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        plot.title = element_text(size=22, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size=20, hjust = 0.5),
        axis.title=element_text(size=20), 
        legend.text=element_text(size=18))+
  theme(legend.position = "bottom")+
  labs(title = "Improved Model rs1_outside_aB Posterior Predictive Check: Test Statistic (Maximum)", x = "Maximum Incidence", y = "Count")

rs_pois_outside_posterior_predictive <- posterior_predict(rs_pois_outside)
ppc_max_rs_pois_outside <-ppc_stat(rbctconf_outside$Incidence, rs_pois_outside_posterior_predictive, stat = "max")+
  xlim(c(0, 400))+ 
  theme(text = element_text(size = 20),axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        plot.title = element_text(size=22, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size=20, hjust = 0.5),
        axis.title=element_text(size=20), 
        legend.text=element_text(size=18))+
  theme(legend.position = "bottom")+
  labs(title = "Poisson GLM Posterior Predictive Check: Test Statistic (Maximum)", x = "Maximum Incidence", y = "Count")
ppc_max_rs_pois_outside
grid.arrange(ppc_max_rs1_outside_aB, ppc_max_rs1_outside_aB_improved, ppc_max_rs_pois_outside, nrow = 3)
col_max_ppcs_initial <- arrangeGrob(ppc_max_rs1_outside_aB, ppc_max_rs1_outside_aB_improved, ppc_max_rs_pois_outside, nrow = 3)
ggsave(col_max_ppcs_initial, file = file.path(out.dir2, "col_max_ppcs_initial.pdf"), h = 24, w = 15)
dim(rs_pois_outside_posterior_predictive)
rs_pois_outside_pp <- pp_check(rs_pois_outside, plotfun = "dens_overlay")+ 
  theme(text = element_text(size = 20),axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        plot.title = element_text(size=22, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size=20, hjust = 0.5),
        axis.title=element_text(size=20), 
        legend.text=element_text(size=18))+
  theme(legend.position = "bottom")+
  labs(title = "Posterior Predictive Check of Poisson Model", x = "Incidence", y = "Density")+
  xlim(c(0, 200))
rs_pois_outside_pp
ggsave(rs_pois_outside_pp, file = file.path(out.dir2, "rs_pois_outside_pp.pdf"), h = 14, w = 8)

grid.arrange(rs1_outside_aB_pp, rs1_outside_aB_improved_pp, rs2_outsideB_pp, rs_pois_outside_pp, nrow = 2)
grid_pps_initial_period <- arrangeGrob(rs1_outside_aB_pp, rs1_outside_aB_improved_pp, rs2_outsideB_pp, rs_pois_outside_pp, nrow = 2)
ggsave(grid_pps_initial_period, file = file.path(out.dir2, "grid_pps_initial_period.pdf"), h = 14, w = 16)

draws_rs1B_improved
rs1_outside_B_improved_post_plot <- mcmc_areas(draws_rs1_outside_B_improved, pars = "TreatmentProactive", prob = 0.95) + xlab('Log Effect')+
  ggplot2::theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20))+
  theme(plot.subtitle = element_text(hjust = 0.5, size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_blank(),
        text = element_text(size = 18))
rs1_outside_B_improved_post_plot
ggsave(rs1_outside_B_improved_post_plot, file = file.path(out.dir2, "rs1_outside_B_improved_post_plot.pdf"), h = 8, w = 16)







rs1_outside_aB_improved_post_plot <- mcmc_areas(draws_rs1_outside_aB_improved, pars = "TreatmentProactive", prob = 0.95) + xlab('Log Effect')+
  ggplot2::theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20))+
  theme(plot.subtitle = element_text(hjust = 0.5, size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_blank(),
        text = element_text(size = 18))+
  labs(title = "Improved Model rs1_outside_aB: Posterior Distribution of Treatment Parameter", subtitle = "with medians and 95% credible intervals")
rs1_outside_aB_improved_post_plot
ggsave(rs1_outside_aB_improved_post_plot, file = file.path(out.dir2, "rs_1aBimproved_after_post_plot.pdf"), h = 8, w = 16)
draws_rs1_outside_aB_after_improved

rs_outside_improved_follow_post_plot <- mcmc_areas(draws_rs_follow_improved, pars = "TreatmentProactive", prob = 0.95) + xlab('Log Effect')+
  ggplot2::theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20))+
  theme(plot.subtitle = element_text(hjust = 0.5, size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_blank(),
        text = element_text(size = 18))+
  labs(title = "Improved Model rs: Posterior Distribution of Treatment Parameter", subtitle = "with medians and 95% credible intervals")
rs_outside_improved_follow_post_plot
ggsave(rs_outside_improved_follow_post_plot, file = file.path(out.dir2, "rs_outside_improved_follow_post_plot.pdf"), h = 8, w = 16)
sum(rbctconf_outside$yrs_since_initialcull)
head(rbctconf_outside_follow)
loo_compare(loo_rs_pois_outside, loo_rs1_outside_aB_improved)
loo_rs_pois_outside
loo_rs1_outside_aB_improved
draws_rs1_outside_aB_follow_improved

rs1_outside_aB_improved_follow_post_plot <- mcmc_areas(draws_rs1_outside_aB_follow_improved, pars = "TreatmentProactive", prob = 0.95) + xlab('Log Effect')+
  ggplot2::theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20))+
  theme(plot.subtitle = element_text(hjust = 0.5, size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_blank(),
        text = element_text(size = 18))+
  labs(title = "Improved Model rs1_outside_aB: Posterior Distribution of Treatment Parameter", subtitle = "with medians and 95% credible intervals")
rs1_outside_aB_improved_follow_post_plot
ggsave(rs1_outside_aB_improved_follow_post_plot, file = file.path(out.dir2, "rs_1aBimproved_post_plot.pdf"), h = 8, w = 16)



#rs_pois_outside post plots-----
rs_pois_outside_post_plot <- mcmc_areas(draws_rs_pois_outside, pars = "TreatmentProactive", prob = 0.95) + xlab('Log Effect')+
  ggplot2::theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20))+
  theme(text = element_text(size = 22), 
        plot.title = element_text(hjust = 0.5, face = "bold", size = 24))+
  theme(plot.subtitle = element_text(hjust = 0.5, size = 18),
        axis.text.x = element_text(size = 28),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 28)
  )
rs_pois_outside_post_plot
ggsave(rs_pois_outside_post_plot, file = file.path(out.dir2, "rs_pois_outside_post_plot.pdf"), h = 8, w = 16)




rs_pois_outside_follow_post_plot <- mcmc_areas(draws_rs_pois_outside_follow, pars = "TreatmentProactive", prob = 0.95) + xlab('Log Effect')+
  ggplot2::theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20))+
  theme(text = element_text(size = 22), 
        plot.title = element_text(hjust = 0.5, face = "bold", size = 24))+
  theme(plot.subtitle = element_text(hjust = 0.5, size = 18),
        axis.text.x = element_text(size = 28),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 28)
  )
# +
# labs(title = "Poisson Model: Posterior Distribution of Treatment Parameter", subtitle = "with medians and 95% credible intervals")
rs_pois_outside_follow_post_plot
ggsave(rs_pois_outside_follow_post_plot, file = file.path(out.dir2, "rs_pois_outside_follow_post_plot.pdf"), h = 8, w = 16)




rs_pois_outside_after_post_plot <- mcmc_areas(draws_rs_pois_outside_after, pars = "TreatmentProactive", prob = 0.95) + xlab('Log Effect')+
  ggplot2::theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20))+
  theme(text = element_text(size = 22), 
        plot.title = element_text(hjust = 0.5, face = "bold", size = 24))+
  theme(plot.subtitle = element_text(hjust = 0.5, size = 18),
        axis.text.x = element_text(size = 28),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 28)
  )
# +
#   labs(title = "Poisson Model: Posterior Distribution of Treatment Parameter", subtitle = "with medians and 95% credible intervals")
rs_pois_outside_after_post_plot
ggsave(rs_pois_outside_after_post_plot, file = file.path(out.dir2, "rs_pois_outside_after_post_plot.pdf"), h = 8, w = 16)





#Model 1aB (c.2) and Poisson Effects Together -----
tmp <- copy(draws_rs1_outside_aB_improved)
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

tmp <- copy(draws_rs_pois_outside)
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
# df_pois[, quant:= quant+3]
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
ggsave(pois_1ab_initial_plot, file = file.path(out.dir2, "pois_1ab_initial_plot_outside.pdf"), h = 12, w = 24)  
out.dir2 


#Model 1aB and Poisson Effects Together (follow-up) ----
tmp <- copy(draws_rs1_outside_aB_follow_improved)
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

tmp <- copy(draws_rs_pois_outside_follow)
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
# df_pois[, quant:= quant+3]
df <- rbind(df, df_pois)
df[, variable:= c(rep("nb", nrow(df)/2),
                  rep("pois", nrow(df)/2))]
pois_1ab_follow_plot_outside <- ggplot(
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
pois_1ab_follow_plot_outside <- pois_1ab_follow_plot_outside+ guides(color = guide_legend(override.aes = list(linewidth = 4)))

ggsave(pois_1ab_follow_plot_outside, file = file.path(out.dir2, "pois_1ab_follow_plot_outside.pdf"), h = 12, w = 24) 

rs1_outside
rs1_outside_aB_after_improved
#Model 1aB and Poisson Effects Together (post-trial) ----
tmp <- copy(draws_rs1_outside_aB_after_improved)
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

tmp <- copy(draws_rs_pois_outside_after)
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
# df_pois[, quant:= quant+3]
df <- rbind(df, df_pois)
df[, variable:= c(rep("nb", nrow(df)/2),
                  rep("pois", nrow(df)/2))]
pois_rs_after_plot_outside <- ggplot(
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
pois_rs_after_plot_outside <- pois_rs_after_plot_outside + guides(color = guide_legend(override.aes = list(linewidth = 4)))
pois_rs_after_plot_outside
ggsave(pois_rs_after_plot_outside, file = file.path(out.dir2, "pois_rs_after_plot_outside.pdf"), h = 12, w = 24) 





length(which(draws_rs_pois_outside_after$TreatmentProactive > 0))/nrow(draws_rs_pois_outside_after)
length(which(draws_rs_outside_follow_improved$TreatmentProactive > 0))/nrow(draws_rs_outside)
length(which(draws_rs_outside_after$TreatmentProactive > 0))/nrow(draws_rs_outside)

length(which(draws_rs_outside_after$TreatmentProactive > 0))/nrow(draws_rs_outside)
length(which(draws_rs_outside_after_improved$TreatmentProactive > 0))/nrow(draws_rs_outside)
length(which(draws_rs1_outside_B_after$TreatmentProactive > 0))/nrow(draws_rs_outside)
length(which(draws_rs1_outside_B_after_improved$TreatmentProactive > 0))/nrow(draws_rs_outside)

length(which(draws_rs1_outside_aB_after$TreatmentProactive > 0))/nrow(draws_rs1_outside_aB_after)
length(which((exp(draws_rs_pois_outside_after$TreatmentProactive) - 1) < 0.15))/nrow(draws_rs_outside)
length(which((exp(draws_rs1_outside_aB_after_improved$TreatmentProactive) - 1) > 0))/nrow(draws_rs_outside)

quantile(exp(draws_rs_outside$TreatmentProactive), probs = c(0.025, 0.5, 0.975))-1
quantile(exp(draws_rs_outside_follow$TreatmentProactive), probs = c(0.025, 0.5, 0.975))-1
quantile(exp(draws_rs_outside_after$TreatmentProactive), probs = c(0.025, 0.5, 0.975))-1
loo(rs_outside_after, k_threshold = 0.7)
loo(rs_outside_after_improved, k_threshold = 0.7)
loo(rs1_outside_B_after, k_threshold = 0.7)
loo(rs1_outside_B_after_improved, k_threshold = 0.7)
loo(rs1_outside_aB_after, k_threshold = 0.7)
loo(rs1_outside_aB_after_improved, k_threshold = 0.7)
loo(rs2_outsideB_after, k_threshold = 0.7)
loo(rs2_outsideB_after_improved, k_threshold = 0.7)


loo(rs_outside_follow, k_threshold = 0.7)
rs_outside_follow_improved$call
loo(rs_outside_follow_improved, k_threshold = 0.7)
loo(rs1_outside_B_follow, k_threshold = 0.7)
loo(rs1_outside_B_follow_improved, k_threshold = 0.7)
loo(rs1_outside_aB_follow, k_threshold = 0.7)
loo(rs1_outside_aB_follow_improved, k_threshold = 0.7)
loo(rs2_outsideB_follow, k_threshold = 0.7)
loo(rs2_outsideB_follow_improved, k_threshold = 0.7)






#Posterior Median Estimates of exponentiated (log-scale) treatment parameter -----
quantile(exp(draws_rs_outside$TreatmentProactive), probs = c(0.025, 0.5, 0.975))-1
quantile(exp(draws_rs_outside_after$TreatmentProactive), probs = c(0.025, 0.5, 0.975))-1

quantile(exp(draws_rs_outside_improved$TreatmentProactive), probs = c(0.025, 0.5, 0.975))-1
quantile(exp(draws_rs_outside_follow$TreatmentProactive), probs = c(0.025, 0.5, 0.975))-1
quantile(exp(draws_rs_outside_follow_improved$TreatmentProactive), probs = c(0.025, 0.5, 0.975))-1
quantile(exp(draws_rs_outside_after_improved$TreatmentProactive), probs = c(0.025, 0.5, 0.975))-1

quantile(exp(draws_rs1_outside_B$TreatmentProactive), probs = c(0.025, 0.5, 0.975))-1
quantile(exp(draws_rs1_outside_B_follow$TreatmentProactive), probs = c(0.025, 0.5, 0.975))-1
quantile(exp(draws_rs1_outside_B_after$TreatmentProactive), probs = c(0.025, 0.5, 0.975))-1

quantile(exp(draws_rs1_outside_B_after_improved$TreatmentProactive), probs = c(0.025, 0.5, 0.975))-1
quantile(exp(draws_rs1_outside_B_follow_improved$TreatmentProactive), probs = c(0.025, 0.5, 0.975))-1
quantile(exp(draws_rs1_outside_B_improved$TreatmentProactive), probs = c(0.025, 0.5, 0.975))-1


quantile(exp(draws_rs1_outside_aB$TreatmentProactive), probs = c(0.025, 0.5, 0.975))-1
quantile(exp(draws_rs1_outside_aB_follow$TreatmentProactive), probs = c(0.025, 0.5, 0.975))-1
quantile(exp(draws_rs1_outside_aB_after$TreatmentProactive), probs = c(0.025, 0.5, 0.975))-1
plot(density(exp(draws_rs1_outside_aB_after$TreatmentProactive)- 1))
quantile(exp(draws_rs1_outside_aB_improved$TreatmentProactive), probs = c(0.025, 0.5, 0.975))-1
quantile(exp(draws_rs1_outside_aB_follow_improved$TreatmentProactive), probs = c(0.025, 0.5, 0.975))-1
quantile(exp(draws_rs1_outside_aB_after_improved$TreatmentProactive), probs = c(0.025, 0.5, 0.975))-1

quantile(exp(draws_rs_pois_outside$TreatmentProactive), probs = c(0.025, 0.5, 0.975))-1
quantile(exp(draws_rs_pois_outside_follow$TreatmentProactive), probs = c(0.025, 0.5, 0.975))-1
quantile(exp(draws_rs_pois_outside_after$TreatmentProactive), probs = c(0.025, 0.5, 0.975))-1

