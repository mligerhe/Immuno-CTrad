# install.packages('ResourceSelection')
# install.packages('glmnet')
# install.packages('doParallel')
# install.packages('pROC')
# install.packages('verification')
# install.packages('fastDummies')
# install.packages('survival')
# install.packages('survminer')
# install.packages('doParallel')
# install.packages('car')
# install.packages('cutpointr')
# install.packages('epiR')
# install.packages('dplyr')
library(ResourceSelection)
library(glmnet)
library(doParallel)
library(pROC)
library(verification)
library(fastDummies)
library(survival)
library(survminer)
library(doParallel)
library(car)
library(cutpointr)
library(epiR)
library(dplyr)
# Training model cohort 1 ####

# Load Phase I radiomics features #
features = read.csv('.../Data/PhaseI_RAD-clin_features.csv')

# Apply elastic-net selection to improve model performance and avoid multicolinearity
mdlY1 <-as.factor(as.matrix(features["outcome"]))
mdlX1 <- as.matrix(data.frame(features[2:9]))

set.seed(97)
a <- seq(0.1, 1, 0.1)
search <- foreach(i = a, .combine = rbind)%do%{
  cv <- cv.glmnet(mdlX1, mdlY1, family = "binomial", nfold = 4 , type.measure = "class", parallel = TRUE, alpha = 1)
  data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se], lambda.1se = cv$lambda.1se, alpha = 1)
}
cv3 <- search[search$cvm == min(search$cvm), ]
cv3
md_fs <- glmnet(mdlX1,mdlY1, family = "binomial", lambda =cv3$lambda.1se[1],alpha =cv3$alpha[1])
coef(md_fs) #model same kVP no location
myCoefs <- coef(md_fs)
vars_nonzero = myCoefs@Dimnames[[1]][ which(myCoefs != 0 ) ][-1] #select variables with non 0 coef

md_rad_clin = glm(outcome~., data = features[c('outcome',vars_nonzero)], family = binomial('logit'))
summary(md_rad_clin) 


df_cutoff = data.frame(outcome = features$outcome, predicted = as.numeric(predict(md_rad_clin,features, type = "response")))

p <- cutpointr(df_cutoff, predicted,outcome, 
               metric = sens_constrain)


p$optimal_cutpoint
R_1 = pROC::roc(features$outcome, md_rad_clin$fitted.values, ci=TRUE)
#RESULTS
R_1$auc;R_1$ci;
coords_1 = coords(R_1,input = 'threshold', x = p$optimal_cutpoint,
                      ret=c("threshold", "sensitivity","specificity", "npv" , "ppv","tp", "tn","fp", "fn"),
                      transpose = FALSE);roc.area(as.numeric(R_1$response)-1, 
                                                  R_1$predictor)$p.value

features$radscore = R_1$predictor
features$rad_class =R_1$predictor>p$optimal_cutpoint
pred_table = table(features$rad_class, features$outcome);pred_table
pred_table = DescTools::Rev(pred_table)
rval <- epi.tests(pred_table, conf.level = 0.95); rval
# Validation GU ####

GU_db = read.csv('.../Data/GU_RAD-clin_features.csv')
names(GU_db)[10] = 'outcome'

R_2= pROC::roc(GU_db$outcome, as.numeric(predict(md_rad_clin,GU_db, type = "response")), ci=TRUE)
R_2$auc;R_2$ci;
coords_2 = coords(R_2,input = 'threshold', x = p$optimal_cutpoint,
                      ret=c("threshold", "sensitivity","specificity", "npv" , "ppv","tp", "tn","fp", "fn",'youden'),
                      transpose = FALSE);roc.area(as.numeric(R_2$response)-1, R_2$predictor)$p.value

GU_db$radclinscore = R_2$predictor
GU_db$classifier = R_2$predictor>p$optimal_cutpoint
pred_table = table(GU_db$classifier, GU_db$outcome)
pred_table = DescTools::Rev(pred_table)
rval2 <- epi.tests(pred_table, conf.level = 0.95); rval

# Validation Lung ####

Lung_db = read.csv('.../Data/Lung_RAD-clin_features.csv')
names(Lung_db)[10] = 'outcome'

#Rad+clin validation Lung

R_3= pROC::roc(Lung_db$outcome, as.numeric(predict(md_rad_clin,Lung_db, type = "response")), ci=TRUE)
R_3$auc;R_3$ci;
coords_3 = coords(R_3,input = 'threshold', x = p$optimal_cutpoint,
                      ret=c("threshold", "sensitivity","specificity", "npv" , "ppv","tp", "tn","fp", "fn"),
                      transpose = FALSE);roc.area(as.numeric(R_3$response)-1, 
                                                  R_3$predictor)$p.value


Lung_db$radclinscore = R_3$predictor
Lung_db$classifier = R_3$predictor>p$optimal_cutpoint
pred_table = table(Lung_db$classifier, Lung_db$outcome)
pred_table = DescTools::Rev(pred_table)
rval3 <- epi.tests(pred_table, conf.level = 0.95); rval

# plot ROCs Clin+Rad #### 
setwd('.../Figures')
tiff("Figure_3_B.tiff", width = 7, height = 7, units = 'in', res = 600)
par(pty = 's')
col_1 = rgb(red =254,green = 154,blue = 1,maxColorValue = 255)
col_2 = rgb(red =0,green = 128,blue = 128,maxColorValue = 255)
col_3 = rgb(red =153,green = 51,blue = 101,maxColorValue = 255)

plot(R_1, xaxt = 'n', legacy.axes = TRUE, ci=!is.null(R_1$ci),
     col = col_1,type='l', xlab = "", ylab = "",
     print.auc = FALSE, print.thres=FALSE, lwd=3.5,
     main= 'Radiomics-Clinical Signature',cex.main = 1.5,
     print.thres.pch=19,print.thres.adj=c(-1.05,0.3),
     print.thres.cex=1.3,print.thres.col="black",
     add=FALSE,cex.lab = 1.8, cex.axis = 1.8, ylim = c(0,1.1))
legend('bottomright', legend = c(paste('Training: ',format(round(R_1$auc,2),nsmall = 2),'(',format(round(R_1$ci[1],2),nsmall = 2),'-',format(round(R_1$ci[3],2),nsmall = 2),')',sep=''),
                                 paste('Bladder: ',format(round(R_2$auc,2),nsmall = 2),'(',format(round(R_2$ci[1],2),nsmall = 2),'-',format(round(R_2$ci[3],2),nsmall = 2),')',sep = ''),
                                 paste('Lung:',format(round(R_3$auc,2),nsmall = 2),'(',format(round(R_3$ci[1],2),nsmall = 2),'-',format(round(R_3$ci[3],2),nsmall = 2),')',sep='')),
       lty = c(1,5,3),lwd = c(4,4,4),col=c(col_1,col_2,col_3), xjust = 4, yjust = 1,
       title = "AUC(95% CI)", pch = c(NA,NA,NA), cex = 1.5)

axis(side = 1,seq(1.0,0.0,-0.2),pos = -0.044,c('0.0','0.2','0.4','0.6','0.8','1.0'),cex.axis = 1.8)

mtext(text = "1 - Specificity",side = 1,line = 3.5, cex = 2)

mtext(text = "Sensitivity",
      side = 2, line = 4, cex = 2)

plot(R_2, xaxt = 'n', legacy.axes = TRUE, ci=!is.null(R_2$ci),
     col = col_2, type = 'l',lty=5,lwd=3.5,
     print.auc = FALSE, print.thres=FALSE, print.thres.col="black",
     print.thres.pch=19,print.thres.adj=c(1,-0.7),print.thres.cex=1.3,add=TRUE,ylim = c(0,1.1))

plot(R_3, xaxt = 'n', legacy.axes = TRUE, ci=!is.null(R_3$ci),
     col = col_3,type = 's',lty=3,lwd=3.5,
     print.auc = FALSE, print.thres=FALSE,print.thres.pch=19,print.thres.adj=c(0.7,-1),
     print.thres.cex=1.3,print.thres.col="black",add=TRUE,ylim=c(0,1.05))
dev.off()
roc.area(as.numeric(R_1$response)-1,R_1$predictor)$p.value
roc.area(as.numeric(R_2$response)-1,R_2$predictor)$p.value
roc.area(as.numeric(R_3$response)-1,R_3$predictor)$p.value

# Clinical-only model training ####
md1_clin <- glm(outcome~baseline_lymphos+baseline_albumin, data = features, family = binomial(link='logit'))
coef(md1_clin) #model same kVP no location
df_cutoff_clin = data.frame(outcome = features["outcome"], predicted = as.numeric(predict(md1_clin, features, type = "response")))

p_clin <- cutpointr(df_cutoff_clin, predicted,outcome, 
                    metric = sens_constrain)


R_1_clin= pROC::roc(features[,10], as.numeric(predict(md1_clin,features[c('baseline_lymphos','baseline_albumin','baseline_num_met_organs')], type = "response")), ci=TRUE)

p_clin$optimal_cutpoint
R_1_clin$auc;R_1_clin$ci;
coords_clin1 = coords(R_1_clin,input = 'threshold', x = p_clin$optimal_cutpoint,
                                ret=c("threshold", "sensitivity","specificity", "npv" , "ppv","tp", "tn","fp", "fn"),
                                transpose = FALSE);roc.area(as.numeric(R_1_clin$response)-1, 
                                                            R_1_clin$predictor)$p.value

features$clinradscore = R_1_clin$predictor
features$clinrad_class =R_1_clin$predictor>p_clin$optimal_cutpoint
pred_table = table(features$clinrad_class, features$outcome);pred_table
pred_table = DescTools::Rev(pred_table)
rval_clin1 <- epi.tests(pred_table, conf.level = 0.95); rval_clin1

#Clin-only GU
R_2_clin = pROC::roc(GU_db$outcome, as.numeric(predict(md1_clin,GU_db[names(md1_clin$coefficients)[-1]], type = "response")), ci=TRUE)
R_2_clin$auc;R_2_clin$ci;
coords_clin2 = coords(R_2_clin,input = 'threshold', x = p_clin$optimal_cutpoint,
                                ret=c("threshold", "sensitivity","specificity", "npv" , "ppv","tp", "tn","fp", "fn"),
                                transpose = FALSE);roc.area(as.numeric(R_2_clin$response)-1,R_2_clin$predictor)$p.value

GU_db$clinradscore = R_2_clin$predictor
GU_db$clinrad_class =R_2_clin$predictor>p_clin$optimal_cutpoint
pred_table = table(GU_db$clinrad_class, GU_db$outcome);pred_table
pred_table = DescTools::Rev(pred_table)
rval_clin2 <- epi.tests(pred_table, conf.level = 0.95); rval_clin2


#Clin-only Lung
R_3_clin = pROC::roc(Lung_db$outcome, as.numeric(predict(md1_clin,Lung_db[names(md1_clin$coefficients)[-1]], type = "response")), ci=TRUE)
R_3_clin$auc;R_3_clin$ci;
coords_clin3 = coords(R_3_clin,input = 'threshold', x = p_clin$optimal_cutpoint,
                                ret=c("threshold", "sensitivity","specificity", "npv" , "ppv","tp", "tn","fp", "fn"),
                                transpose = FALSE);roc.area(as.numeric(R_3_clin$response)-1,R_3_clin$predictor)$p.value

Lung_db$clinradscore = R_3_clin$predictor
Lung_db$clinrad_class =R_3_clin$predictor>p_clin$optimal_cutpoint
pred_table = table(Lung_db$clinrad_class, Lung_db$outcome);pred_table
pred_table = DescTools::Rev(pred_table)
rval_clin3 <- epi.tests(pred_table, conf.level = 0.95); rval_clin3

# Compare clinical+radiomics and radiomics model####
clinonly.model <- glm(outcome ~baseline_lymphos+baseline_albumin+baseline_LDH, data = features, family = binomial)
radonly.model <- glm(outcome ~Radscore, data = features, family = binomial)
radclin.model <- glm(outcome ~Radscore+baseline_lymphos+baseline_albumin+baseline_LDH, data = features, family = binomial)


AIC_c = AIC(clinonly.model)
AIC_rc = AIC(radclin.model)
AIC_rad = AIC(radonly.model)
AIC_c
AIC_rad
AIC_rc

BIC_c = BIC(clinonly.model)
BIC_rc = BIC(radclin.model)
BIC_rad = BIC(radonly.model)
BIC_c
BIC_rad
BIC_rc

anova(radonly.model,radclin.model, test='Chisq')

hl <- hoslem.test(md_rad_clin$y, fitted(md_rad_clin), g=5)
hl

# plot ROCs Clin only #### 
setwd('.../Figures')
tiff("Figure_E1.tiff", width = 7, height = 7, units = 'in', res = 600)
par(pty = 's')
col_1 = rgb(red =254,green = 154,blue = 1,maxColorValue = 255)
col_2 = rgb(red =0,green = 128,blue = 128,maxColorValue = 255)
col_3 = rgb(red =153,green = 51,blue = 101,maxColorValue = 255)

plot(R_1_clin, xaxt = 'n', legacy.axes = TRUE, ci=!is.null(R_1$ci),
     col = col_1,type='l', xlab = "", ylab = "",
     print.auc = FALSE, print.thres=FALSE, lwd=3.5,
     main= 'Clinical Model',cex.main = 1.5,
     print.thres.pch=19,print.thres.adj=c(-1.05,0.3),
     print.thres.cex=1.3,print.thres.col="black",
     add=FALSE,cex.lab = 1.8, cex.axis = 1.8, ylim = c(0,1.1))
legend('bottomright', legend = c(paste('Training: ',format(round(R_1_clin$auc,2),nsmall = 2),'(',format(round(R_1_clin$ci[1],2),nsmall = 2),'-',format(round(R_1_clin$ci[3],2),nsmall = 2),')',sep=''),
                                 paste('Bladder: ',format(round(R_2_clin$auc,2),nsmall = 2),'(',format(round(R_2_clin$ci[1],2),nsmall = 2),'-',format(round(R_2_clin$ci[3],2),nsmall = 2),')',sep = ''),
                                 paste('Lung:',format(round(R_3_clin$auc,2),nsmall = 2),'(',format(round(R_3_clin$ci[1],2),nsmall = 2),'-',format(round(R_3_clin$ci[3],2),nsmall = 2),')',sep='')),
       lty = c(1,5,3),lwd = c(4,4,4),col=c(col_1,col_2,col_3), xjust = 4, yjust = 1,
       title = "AUC(95% CI)", pch = c(NA,NA,NA), cex = 1.5)

R_1_clin$auc;R_1_clin$ci;coords(R_1_clin,input = 'threshold', x = 0.4219269,
                                ret=c("threshold", "specificity", "sensitivity", "tp", "tn","fp", "fn"),
                                transpose = FALSE);roc.area(as.numeric(R_1_clin$response)-1, 
                                                            R_1_clin$predictor)$p.value



axis(side = 1,seq(1.0,0.0,-0.2),pos = -0.044,c('0.0','0.2','0.4','0.6','0.8','1.0'),cex.axis = 1.8)

mtext(text = "1 - Specificity",side = 1,line = 3.5, cex = 2)

mtext(text = "Sensitivity",
      side = 2, line = 4, cex = 2)

plot(R_2_clin, xaxt = 'n', legacy.axes = TRUE, ci=!is.null(R_2_clin$ci),
     col = col_2, type = 'l',lty=5,lwd=3.5,
     print.auc = FALSE, print.thres=FALSE, print.thres.col="black",
     print.thres.pch=19,print.thres.adj=c(1,-0.7),print.thres.cex=1.3,add=TRUE,ylim = c(0,1.1))
R_2_clin$auc;R_2_clin$ci;coords(R_2_clin,input = 'threshold', x = 0.4219269,
                                ret=c("threshold", "specificity", "sensitivity", "tp", "tn","fp", "fn"),
                                transpose = FALSE);roc.area(as.numeric(R_2_clin$response)-1, 
                                                            R_2_clin$predictor)$p.value


plot(R_3_clin, xaxt = 'n', legacy.axes = TRUE, ci=!is.null(R_3_clin$ci),
     col = col_3,type = 's',lty=3,lwd=3.5,
     print.auc = FALSE, print.thres=FALSE,print.thres.pch=19,print.thres.adj=c(0.7,-1),
     print.thres.cex=1.3,print.thres.col="black",add=TRUE,ylim=c(0,1.05))
R_3_clin$auc;R_3_clin$ci;coords(R_3_clin,input = 'threshold', x = 0.4219269,
                                ret=c("threshold", "specificity", "sensitivity", "tp", "tn","fp", "fn"),
                                transpose = FALSE);roc.area(as.numeric(R_3_clin$response)-1, 
                                                            R_3_clin$predictor)$p.value

dev.off()
roc.area(as.numeric(R_1_clin$response)-1,R_1_clin$predictor)$p.value
roc.area(as.numeric(R_2_clin$response)-1,R_2_clin$predictor)$p.value
roc.area(as.numeric(R_3_clin$response)-1,R_3_clin$predictor)$p.value



# OS PhaseI analysis ####
get_hr = function(mdc){
  HR = exp(coef(mdc))
  HRCIlow = exp(confint(mdc))[,1]
  HRCIhigh = exp(confint(mdc))[,2]
  p_val_model = summary(mdc)$coefficients[,5]
  aic= AIC(mdc)
  c_idx= mdc$concordance[6]
  c_idx_low= mdc$concordance[6] - 1.96*mdc$concordance[7]
  c_idx_upp= mdc$concordance[6] +1.96*mdc$concordance[7]
  temp = cox.zph(mdc)
  data.frame(HR = paste(round(HR,2),'[',round(HRCIlow,2),'-',round(HRCIhigh,2),']',sep=''), p_val_model, aic,
             CI = paste(round(c_idx,2),'[',round(c_idx_low,2),'-',round(c_idx_upp,2),']',sep=''),p_val = temp$table[c(1:length(temp$table[,1])-1),3])
}
features$status_at_last_FU = as.numeric(features$status_at_last_FU)-1

#Multivariate and univariate studies with Overall survival
formula_radclinmodel= Surv(OS_months, status_at_last_FU) ~ Radscore+baseline_lymphos+baseline_albumin
model_clinradmodel = coxph(formula_radclinmodel, features)
HR_multivarmd = get_hr(model_clinradmodel)

formula_radscore = Surv(OS_months, status_at_last_FU) ~ Radscore
model_clinscore = coxph(formula_radscore, features)
HR_radmodel = get_hr(model_clinscore)

features$radclinscore = R_1$predictor
formula_clinscore = Surv(OS_months, status_at_last_FU) ~ radclinscore
model_clinscore = coxph(formula_clinscore, features)
HR_radclinscore = get_hr(model_clinscore)

univar_db= foreach(i = 2:9, .combine = 'rbind')%do%{
  var_f = names(features)[i]
  form =as.formula(sprintf('Surv(OS_months, status_at_last_FU) ~%s', var_f))
  mdc = coxph(form,features)
  cbind(var_f,get_hr(mdc))
}

# OS GU analysis ####
GU_db$status_at_last_FU = as.numeric(GU_db$status_at_last_FU)-1
model_clinradmodel_GU = coxph(formula_radclinmodel, GU_db)
HR_clinradmodel_GU = get_hr(model_clinradmodel_GU)

model_radmodel_GU = coxph(formula_radscore, GU_db)
HR_radmodel_GU = get_hr(model_radmodel_GU)

GU_db$radclinscore = R_2$predictor
model_clinscore_GU = coxph(formula_clinscore, GU_db)
HR_radclinscore_GU = get_hr(model_clinscore_GU)


col_1 = rgb(red =52,green = 51,blue = 154,maxColorValue = 255)
col_2 = rgb(red =254,green = 0,blue = 0,maxColorValue = 255)
setwd('.../Figures')
tiff("Figure_5_A.tiff", width = 7, height = 7, units = 'in', res = 600)
GU_db$classifier = GU_db$Radscore>p$optimal_cutpoint


wilcox.test(GU_db$OS_months~GU_db$classifier)
fit <- survfit(Surv(OS_months, status_at_last_FU) ~ classifier,
               data = GU_db)

print(summary(fit, times = 24))

ggsurvplot(fit, censor.shape="|", censor.size = 4, palette = c("blue", "red"),
                         ggtheme = theme_bw(),
                         xlim = c(0,43), title = 'Kaplan-Meier Bladder Patients',
                         font.title = c(20,'plain', 'black'), 
                         font.caption = c(14,'plain', 'black'), 
                         font.x = c(20,'plain', 'black'),
                         font.y = c(20,'plain', 'black'),
                         font.tickslab  = c(13,'plain', 'black'), 
                         font.legend  = c(20,'plain', 'black'),
                         legend = 'top',
                         legend.size = 2,
                         legend.title = element_blank(),
                         legend.labs = c("Low RadScore","High RadScore"),
                         xlab = 'Time (months)', size = 1,  
                         data = GU_db, 
                         risk.table = T, 
                         pval = 0.002)


dev.off()


# OS Lung analysis ####
Lung_db$status_at_last_FU = as.numeric(Lung_db$status_at_last_FU)-1
model_clinradmodel_Lung = coxph(formula_radclinmodel, Lung_db)
HR_clinradmodel_Lung = get_hr(model_clinradmodel_Lung)

model_radmodel_Lung = coxph(formula_radscore, Lung_db)
HR_radmodel_Lung = get_hr(model_radmodel_Lung)

Lung_db$radclinscore = R_3$predictor
model_clinscore_Lung = coxph(formula_clinscore, Lung_db)
HR_radclinscore_Lung = get_hr(model_clinscore_Lung)



Lung_db$classifier = Lung_db$radclinscore> p$optimal_cutpoint

col_1 = rgb(red =52,green = 51,blue = 154,maxColorValue = 255)
col_2 = rgb(red =254,green = 0,blue = 0,maxColorValue = 255)
tiff("Figure_5_B.tiff", width = 7, height = 7, units = 'in', res = 600)

Lung_db[Lung_db$OS_months >40,]$OS_months = 40
fit <- survfit(Surv(OS_months, status_at_last_FU) ~ classifier,
               data = Lung_db)

print(summary(fit, times = 24))


ggsurvplot(fit, censor.shape="|", censor.size = 4, palette = c("blue", "red"),
           ggtheme = theme_bw(),
           xlim = c(0,43), title = 'Kaplan-Meier Lung Patients',
           font.title = c(20,'plain', 'black'), 
           font.caption = c(14,'plain', 'black'), 
           font.x = c(20,'plain', 'black'),
           font.y = c(20,'plain', 'black'),
           font.tickslab  = c(13,'plain', 'black'), 
           font.legend  = c(20,'plain', 'black'),
           legend = 'top',
           legend.size = 2,
           legend.title = element_blank(),
           legend.labs = c("Low RadScore","High RadScore"),
           xlab = 'Time (months)', size = 1,  
           data = Lung_db, 
           risk.table = T, 
           pval = T)

dev.off()

HR_models = rbind(Radscore = HR_radmodel,Radscore_GU = HR_radmodel_GU,Radscore_Lung = HR_radmodel_Lung, radclinscore = HR_radclinscore, 
                  radclinscore_GU = HR_radclinscore_GU,radclinscore_Lung = HR_radclinscore_Lung)
