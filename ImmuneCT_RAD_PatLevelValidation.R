library(ggplot2)
features = read.csv('.../Data/PhaseI_RAD-clin_features.csv')
GU_db = read.csv('.../Data/GU_RAD-clin_features.csv')
Lung_db = read.csv('.../Data/Lung_RAD-clin_features.csv')
features$Radscore = features$Radscore*10
GU_db$Radscore = GU_db$Radscore*10
Lung_db$Radscore = Lung_db$Radscore*10

md_1 = glm(outcome~Radscore, data =features, family = binomial(link='logit'))
summary(md_1)
n1_od = unname(exp(coef(md_1))[2])
n1_cil = exp(confint(md_1))[2,1]
n1_cih = exp(confint(md_1))[2,2]

md_2 = glm(clinical_benefit_rate~Radscore, data =GU_db, family = binomial(link='logit'))
summary(md_2)
n2_od = unname(exp(coef(md_2))[2])
n2_cil = exp(confint(md_2))[2,1]
n2_cih = exp(confint(md_2))[2,2]

md_3 = glm(clinical_benefit_rate~Radscore, data =Lung_db, family = binomial(link='logit'))
summary(md_3)
n3_od = unname(exp(coef(md_3))[2])
n3_cil = exp(confint(md_3))[2,1]
n3_cih = exp(confint(md_3))[2,2]


df = data.frame(labels = c(paste('Cohort 1 \n',round(n1_od,2),'(', round(n1_cil,2),'-',round(n1_cih,2),', p=' ,round(summary(md_1)$coefficients[2,4],3),')',sep=''),
                           paste('Cohort 2 \n',round(n2_od,2),'(', round(n2_cil,2),'-',round(n2_cih,2),', p=' ,round(summary(md_2)$coefficients[2,4],3),')',sep=''),
                           paste('Cohort 3 \n',round(n3_od,2),'(', round(n3_cil,2),'-',round(n3_cih,2),', p=' ,round(summary(md_3)$coefficients[2,4],3),')',sep='')),
                odds= c(n1_od,n2_od,n3_od),ci_l= c(n1_cil,n2_cil,n3_cil),ci_h= c(n1_cih,n2_cih,n3_cih))
df
setwd('.../Figures')
tiff("Figure_4.tiff", width = 10, height = 7, units = 'in', res = 600)

ggplot(data = df,
       mapping = aes(y = forcats::fct_inorder(f = rev(x = labels)))) +
  geom_vline(aes(xintercept = 1), size = .5, linetype = "dashed") +
  geom_point(mapping = aes(x = rev(odds)),shape = 'square',size = 6, color = "black") +
  geom_errorbarh(mapping = aes(xmin = rev(ci_l),
                               xmax = rev(ci_h)), height=.1) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        plot.margin = margin(t = 1, b = 1, r=15,l=0)) +
  coord_trans(x = "exp") +
  scale_x_continuous(breaks = c(0.5,seq(1,1.8,0.2)), limits = c(0.5,1.8)) +
  labs(x = "OR",
       y = "")
dev.off()
