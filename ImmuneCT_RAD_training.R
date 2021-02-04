# install.packages('ResourceSelection')
# install.packages('glmnet')
# install.packages('doParallel')
# install.packages('pROC')
# install.packages('fastDummies')
# install.packages('verification')
# install.packages('ggplot2')
# install.packages('cutpointr')
# install.packages('DescTools')
# install.packages('lme4')
# install.packages('epiR')
library(ResourceSelection)
library(glmnet)
library(doParallel)
library(pROC)
library(fastDummies)
library(verification)
library(ggplot2)
library(cutpointr)
library(DescTools)
library(lme4)
library(epiR)
#Load radiomics features
features = read.csv('.../Data/PhaseI_radfeatures_outcome.csv')
#normalized radiomics features
features[3:106] = scale(features[3:106])
# Reproducibility filter (ICC) ####
icc_file1 = ".../Data/IAP_standardize_Repeatability.csv"
icc_file2 = ".../IAP_nonstand_Repeatability.csv"
icc_file3 = ".../Data/interobserver_Repeatability.csv"

icc1 = read.csv2(icc_file1, header = TRUE,sep = ';',row.names=1)
icc2 = read.csv2(icc_file2, header = TRUE,sep = ';',row.names=1)
icc3 = read.csv2(icc_file3, header = TRUE,sep = ';',row.names=1)

var1 = rownames(icc1[icc1$icc > 0.7,])
var2 = rownames(icc2[icc2$icc > 0.7,])
var3 = rownames(icc3[icc3$icc > 0.7,])


final_var = var1[var1%in%var3]
final_var = final_var[final_var%in%var2]
final_var = final_var[final_var %in% names(features)]

features = features[c('Patient', 'LesionNumber', final_var, 'location', 'outcome')]


# Add dummy variables ######
dummy_vars = dummy_cols(na.omit(features$location))

features = cbind(features, dummy_vars[-1])


# Define data and output matrices ####
mdlY1 <-as.factor(as.matrix(features["outcome"]))
mdlX1 <- as.matrix(data.frame(features[c(final_var, names(dummy_vars[-1]))]))

# Elastic net feaure selection ####
set.seed(97)
a <- seq(0.1, 1, 0.1)
search <- foreach(i = a, .combine = rbind)%do%{
  cv <- cv.glmnet(mdlX1, mdlY1, family = "binomial", nfold = 4 , type.measure = "class", parallel = TRUE, alpha = i)
  data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se], lambda.1se = cv$lambda.1se, alpha = i)
}
cv3 <- search[search$cvm == min(search$cvm), ]
cv3

md1 <- glmnet(mdlX1,mdlY1, family = "binomial", lambda =cv3$lambda.1se[2],alpha =cv3$alpha[2] )
coef(md1) #model same kVP no location

vars_nonzero = md1$beta@Dimnames[[1]][md1$beta@i+1]
# Hosmer-Lemeshow Goodness of Fit (GOF) Test
mda_fit = as.numeric(predict(md1, mdlX1, type = "response"))
hl <- hoslem.test(as.numeric(mdlY1)-1, mda_fit, g=15)
hl


# mixed model ####
formula =  as.formula(paste("outcome ~", paste(vars_nonzero, collapse = " + "),"+(1|Patient)"))
formula

md_glmer = glmer(formula, data = features[c('Patient',vars_nonzero, 'outcome')], family = binomial("logit"))
coefs_sd = summary(md_glmer) 
CI95_high = coefs_sd$coefficients[,1]+1.96*coefs_sd$coefficients[,2]
CI95_low = coefs_sd$coefficients[,1]-1.96*coefs_sd$coefficients[,2]
coef = coefs_sd$coefficients[,1]
coefs_md = data.frame(coef = coef[-1], CI_low = CI95_low[-1], ci_high=CI95_high[-1])
coefs_md = data.frame(names = as.character(names(coef)[-1]),vars = paste(round(coef,2), '[',round(CI95_low,2),'-',round(CI95_high,2),']', sep = '')[-1])

#High sensitivity cutoff

df_cutoff = data.frame(outcome = features$outcome, predicted = as.numeric(predict(md_glmer, type = 'response',re.form=NA)))

p <- cutpointr(df_cutoff, predicted,outcome, 
               metric = sens_constrain)
p$optimal_cutpoint

#Roc curve
R_1 = pROC::roc(features$outcome, as.numeric(predict(md_glmer, type = "response",re.form = NA)), ci=TRUE)
ci  = R_1$ci
coords_1 = coords(R_1,input = 'threshold', x = p$optimal_cutpoint,
                  ret=c("threshold","sensitivity","specificity","npv","ppv","tp","tn","fp","fn"),
                  transpose = FALSE)
features$radscore = R_1$predictor
features$rad_class =R_1$predictor>p$optimal_cutpoint
pred_table = table(features$rad_class, features$outcome);pred_table
pred_table = DescTools::Rev(pred_table)
rval <- epi.tests(pred_table, conf.level = 0.95); rval #Confidence intervals SE



plot(R_1, xaxt = 'n', legacy.axes = TRUE, ci=!is.null(R_1$ci),col = 'black',
     print.auc = FALSE,  print.thres=p$optimal_cutpoint,print.thres.pch=19,print.thres.adj=c(1.05,0.4),
     print.thres.cex=1.3,add=FALSE,cex.lab = 1.5, cex.axis = 1.3,ylim = c(0,1.1))
axis(side = 1,seq(1.0,0.0,-0.2),pos = -0.044,c('0.0','0.2','0.4','0.6','0.8','1.0'),cex.lab = 1.5, cex.axis = 1.3)

# Validation GU ####

GU_db = read.csv(".../Data/GU_radfeatures_outcome.csv")
GU_db[4:110] = scale(GU_db[4:110])


# Add dummy variables#
dummy_vars = dummy_cols(GU_db$LesionType)
GU_db = cbind(GU_db, dummy_vars[-c(1)])


# Validate Model #
GU_db$Patient = factor(paste('GU-',GU_db$Patient,sep = ''))

library(fastDummies)
dummy_vars = dummy_cols(GU_db$LesionType)

GU_db = cbind(GU_db, dummy_vars[-c(1)])

GU_db$outcome = as.integer(GU_db$outcome)-1

R_2= pROC::roc(GU_db$outcome, as.numeric(predict(md_glmer,GU_db, type = "response",re.form=NA)), ci=TRUE)
R_2$auc;R_2$ci;
coords_2 = coords(R_2,input = 'threshold', x = p$optimal_cutpoint,
                  ret=c("threshold", "sensitivity","specificity", "npv" , "ppv","tp", "tn","fp", "fn"),
                  transpose = FALSE);roc.area(R_2$response,R_2$predictor)$p.value


GU_db$radscore = R_2$predictor
GU_db$rad_class =R_2$predictor>p$optimal_cutpoint
pred_table = table(GU_db$rad_class, GU_db$outcome);pred_table
pred_table = DescTools::Rev(pred_table)
rval2 <- epi.tests(pred_table, conf.level = 0.95); rval

plot(R_2, xaxt = 'n', legacy.axes = TRUE, ci=!is.null(R_2$ci),col = 'black',lty = 6,
     print.auc = FALSE, grid=c(0.1, 0.2), print.thres=p$optimal_cutpoint,print.thres.pch=19,print.thres.adj=c(-0.05,1.2),
     print.thres.cex=1.3,add=FALSE)
legend(x=-0.9, y=1.08, legend = c(paste('Training: ',round(R_1$auc,2),' (',round(R_1$ci[1],2),' - ',round(R_1$ci[3],2),')',sep=''),
                                  paste('External Validation: ',round(R_2$auc,2),' (',round(R_2$ci[1],2),' - ',round(R_2$ci[3],2),')',sep='')),
       lty = c(1,6),lwd = c(3,3),col=c('black','black'), xjust = 4, yjust = 1,
       title = "AUC(95% CI)", pch = c(NA,NA), cex = 1.3)


# Validation Lung ####
Lung_db = read.csv(".../Data/Lung_radfeatures_outcome.csv")
Lung_db[4:110] = scale(Lung_db[4:110])
Lung_db$Patient = factor(paste('Lung-',Lung_db$Patient,sep = ''))

dummy_vars = dummy_cols(Lung_db$LesionType)
Lung_db = cbind(Lung_db, dummy_vars[-c(1)])

Lung_db$Patient = factor(Lung_db$Patient)
Lung_db$outcome = as.integer(Lung_db$outcome)-1

# Validate Model #

R_3 = pROC::roc(Lung_db$outcome, as.numeric(predict(md_glmer,Lung_db, type = "response",re.form = NA)), ci=TRUE)
R_3$auc;R_3$ci;
coords_3= coords(R_3,input = 'threshold', x = p$optimal_cutpoint,
                 ret=c("threshold", "sensitivity","specificity", "npv" , "ppv","tp", "tn","fp", "fn"),
                 transpose = FALSE);roc.area(as.numeric(R_3$response)-1,R_3$predictor)$p.value
Lung_db$radscore = R_3$predictor
Lung_db$rad_class =R_3$predictor>p$optimal_cutpoint
pred_table = table(Lung_db$rad_class, Lung_db$outcome);pred_table
pred_table = DescTools::Rev(pred_table)
rval3 <- epi.tests(pred_table, conf.level = 0.95); rval

plot(R_3, xaxt = 'n', legacy.axes = TRUE, ci=!is.null(R_3$ci),col = 'black',lty = 6,
     print.auc = FALSE, grid=c(0.1, 0.2), print.thres=0.379677,print.thres.pch=19,print.thres.adj=c(-0.05,1.2),
     print.thres.cex=1.3,add=FALSE)

# Plot ROCs ####
setwd('.../Figures')
tiff("Figure_3_A.tiff", width = 7, height = 7, units = 'in', res = 600)
par(pty = 's')
col_1 = rgb(red =254,green = 154,blue = 1,maxColorValue = 255)
col_2 = rgb(red =0,green = 128,blue = 128,maxColorValue = 255)
col_3 = rgb(red =153,green = 51,blue = 101,maxColorValue = 255)

plot(R_1, xaxt = 'n', legacy.axes = TRUE, ci=!is.null(R_1$ci),
     col = col_1,type='l', xlab = "", ylab = "",
     print.auc = FALSE, print.thres=FALSE, lwd=3.5,
     main= 'Radiomics Signature',cex.main = 1.5,
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
roc.area(R_1$response,R_1$predictor)$p.value
roc.area(R_2$response,R_2$predictor)$p.value
roc.area(R_3$response,R_3$predictor)$p.value
