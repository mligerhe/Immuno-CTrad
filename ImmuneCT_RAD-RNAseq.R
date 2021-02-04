# Libraries ####
library(glmnet)
library(fastDummies)
library(dplyr)
library(ggpubr)
library(beeswarm)
# Load features ######
RNA_db = read.csv('.../Data/RNAseq_CTrad.csv')
# Normalize data ####

RNA_db[3:109] = scale(RNA_db[3:109])
# Add dummy variables  ####
dummy_vars = dummy_cols(RNA_db$LesionType)
RNA_db = cbind(RNA_db, dummy_vars[-c(1)])


# Validate Model ####
load(".../Data/glmner_md.RData")

RNA_rad_db =data.frame(Patient = RNA_db$Patient, LesionType = RNA_db$LesionType,
                 RadScore = predict(md_glmer, RNA_db,re.form = NA, type = 'response'))

# Immune signature ####

rna_db = read.csv('.../Data/RNA_seq_data.csv')

RNA_rad_db = merge(RNA_rad_db,rna_db, by = 'Patient')
RNA_rad_db = na.omit(RNA_rad_db)
RNA_rad_db$RadScore_d = factor(ifelse(RNA_rad_db$RadScore<=0.29,'Low','High'))

RNA_rad_db%>%
  group_by(RadScore_d)%>%
  summarize(mean_OS = mean(CYT),
            sd_OS = sd(CYT))

# CYTOTOXIC boxplot ####
setwd('.../Figures')
tiff("Figure_6.tiff", width = 7, height = 7, units = 'in', res = 600)
par(pty = 's')
boxplot(CYT~RadScore_d,data = RNA_rad_db, names = c('High','Low'),xlab = 'Radiomics Score', 
        ylab = 'Cytotoxic enrichment',add=FALSE, cex.lab = 1.8, cex.axis = 1.8,  main = '', cex.main = 2)
beeswarm(CYT~RadScore_d,data = RNA_rad_db,method = 'swarm',pch = 16, cex = 2, add = T)
pval = wilcox.test(RNA_rad_db$CYT~RNA_rad_db$RadScore_d)$p.value
text(x=1.1, y= 0.4,labels = paste('p=',round(pval,2),sep=''),cex = 1.5)

dev.off()




