#### Script for behavioral/brain analysis of MKM Study: Becker, Repantis, Dresler & K??hn, 2021 ##################
# (c) almaxi@gmail.com

### 0) load libraries ################ 
library(lmerTest)
library(corrplot)
library(sjPlot)
library(sjmisc)
library(knitr)
library(magrittr)
library(sjlabelled)      
library(sjmisc)                                                                                    
library(sjstats) 
library(ggeffects)
library(performance)
library(parameters)
library(betareg)
library(MKinfer)
library(rcompanion)
library(lavaan)
library(psych)
library(GPArotation)
library(glmmTMB)
#library(emmeans)

rm(list = ls())

### 1) load data and rename variables ################
setwd('C:/Users/Maxi/Google Drive/personal/Mx/MKM/data/')
MKM2 <- read.table("BeckerRepantis_et_al_2021_MKM_long.csv", sep = ";", dec = ",", header=TRUE, na.strings=c("", " " , "NA", "NAN" )) 

#rename variables:
MKM2$early_visRecall  = MKM2$recallScanR
MKM2$late_visRecall = MKM2$post_recallScanR 
MKM2$early_audRecall = MKM2$recallAudio1r 
MKM2$late_audRecall = MKM2$post_recallAudioR 
MKM2$early_false_audRecall = MKM2$recallAudio1l
MKM2$late_false_audRecall = MKM2$post_recallAudioL 
MKM2$group = MKM2$t1_med
MKM2[MKM2$group == "Methylphenidat",]$group = "MPH"
MKM2[MKM2$group == "Modafinil",]$group = "MOD"
MKM2[MKM2$group == "Koffein",]$group = "CAF"
MKM2[MKM2$condition == "enhancer",]$condition = "stimulant"
MKM2$group = as.factor(MKM2$group)  
MKM2$med_order = as.factor(MKM2$med_order) 
MKM2$PIN = as.factor(MKM2$PIN)

### 2) Results: Effects of stimulants on memory  ############

###  early_visRecall
mean(MKM2[MKM2$condition != "placebo",]$early_visRecall, na.rm = T)
  sd(MKM2[MKM2$condition != "placebo",]$early_visRecall, na.rm = T)
mean(MKM2[MKM2$condition == "placebo",]$early_visRecall, na.rm = T)
  sd(MKM2[MKM2$condition == "placebo",]$early_visRecall, na.rm = T)
  
M1a <- glmmTMB( early_visRecall ~   group +med_order +(1|PIN),
             data=MKM2 ,family = poisson)
M1b <- glmmTMB( early_visRecall ~   group+condition+med_order  +(1|PIN),
                data=MKM2 ,family = poisson)
anova(M1a, M1b)
plot(check_distribution(M1b))
summary(M1b)
ggM1 <-ggpredict(M1b, c(  "condition", "group")) %>% plot(show.title = F) + ggplot2::theme_classic()

#### late_visRecall 
mean(MKM2[MKM2$condition != "placebo",]$late_visRecall, na.rm = T)
  sd(MKM2[MKM2$condition != "placebo",]$late_visRecall, na.rm = T)
mean(MKM2[MKM2$condition == "placebo",]$late_visRecall, na.rm = T)
  sd(MKM2[MKM2$condition == "placebo",]$late_visRecall, na.rm = T)

M2a <- glmmTMB( late_visRecall ~   group+med_order  +(1|PIN),
                data=MKM2  ,ziformula=~1, family = poisson)
M2b <- glmmTMB( late_visRecall ~   group+condition+med_order  +(1|PIN),
             data=MKM2   ,ziformula=~1,family = poisson) 
  anova(M2a, M2b)
  plot(check_distribution(M2b))
  summary(M2b)
  ggM2 <-ggpredict(M2b, c(  "condition", "group")) %>% plot(show.title = F) + ggplot2::theme_classic()

#### early_audRecall
mean(MKM2[MKM2$condition != "placebo",]$early_audRecall, na.rm = T)
  sd(MKM2[MKM2$condition != "placebo",]$early_audRecall, na.rm = T)
mean(MKM2[MKM2$condition == "placebo",]$early_audRecall, na.rm = T)
  sd(MKM2[MKM2$condition == "placebo",]$early_audRecall, na.rm = T)

M3a <- glmmTMB( early_audRecall ~  group+med_order  +(1|PIN),
             data=MKM2 , family=poisson ) 
M3b <- glmmTMB( early_audRecall ~  group+condition+med_order  +(1|PIN),
               data=MKM2 , family=poisson ) 
  anova(M3a, M3b)
  plot(check_distribution(M3b))
  summary(M3b)
  ggM3 <-ggpredict(M3b, c(  "condition", "group")) %>% plot(show.title = F) + ggplot2::theme_classic()

#### late_audRecall 
mean(MKM2[MKM2$condition != "placebo",]$late_audRecall, na.rm = T)
  sd(MKM2[MKM2$condition != "placebo",]$late_audRecall, na.rm = T)
mean(MKM2[MKM2$condition == "placebo",]$late_audRecall, na.rm = T)
  sd(MKM2[MKM2$condition == "placebo",]$late_audRecall, na.rm = T)

hist((MKM2$late_audRecall))
M5a <- glmmTMB( late_audRecall ~  group+med_order +(1|PIN),
             data=MKM2,  family = poisson  ) 
M5b <- glmmTMB( late_audRecall ~   group+condition+med_order  +(1|PIN),
               data=MKM2,  family = poisson  ) 
  anova(M5a, M5b)
  plot(check_distribution(M5b))
  summary(M5b)
  ggM5 <-ggpredict(M5b, c(  "condition", "group")) %>% plot(show.title = F) + ggplot2::theme_classic()

#### early_false_audRecall - lures
mean(MKM2[MKM2$condition != "placebo",]$early_false_audRecall, na.rm = T)
  sd(MKM2[MKM2$condition != "placebo",]$early_false_audRecall, na.rm = T)
mean(MKM2[MKM2$condition == "placebo",]$early_false_audRecall, na.rm = T)
  sd(MKM2[MKM2$condition == "placebo",]$early_false_audRecall, na.rm = T)

hist((MKM2$early_false_audRecall))
M6a <- glmmTMB( early_false_audRecall ~   group+med_order  +(1|PIN),
             data=MKM2  ,family = poisson)
M6b <- glmmTMB( early_false_audRecall ~  group+condition+med_order   +(1|PIN),
               data=MKM2  ,family = poisson) 
  anova(M6a, M6b)
  plot(check_distribution(M6b))
  summary(M6b)
  ggM6 <- ggpredict(M6b, c(  "condition", "group")) %>% plot(show.title = F) + ggplot2::theme_classic()

#### late_false_audRecall - lures
mean(MKM2[MKM2$condition != "placebo",]$late_false_audRecall, na.rm = T)
  sd(MKM2[MKM2$condition != "placebo",]$late_false_audRecall, na.rm = T)
mean(MKM2[MKM2$condition == "placebo",]$late_false_audRecall, na.rm = T)
  sd(MKM2[MKM2$condition == "placebo",]$late_false_audRecall, na.rm = T)

hist((MKM2$late_false_audRecall))
M7a <- glmmTMB( late_false_audRecall ~   group+med_order   +(1|PIN), data=MKM2 )
M7b <- glmmTMB( late_false_audRecall ~   group+condition +med_order  +(1|PIN), data=MKM2 )
  anova(M7a, M7b)
  plot(check_distribution(M7b))
  summary(M7b)
  ggM7 <- ggpredict(M7b, c(  "condition", "group")) %>% plot(show.title = F) + ggplot2::theme_classic()


#### d_prime
  mean(MKM2[MKM2$condition != "placebo",]$d_prime, na.rm = T)
    sd(MKM2[MKM2$condition != "placebo",]$d_prime, na.rm = T)
  mean(MKM2[MKM2$condition == "placebo",]$d_prime, na.rm = T)
    sd(MKM2[MKM2$condition == "placebo",]$d_prime, na.rm = T)
  
  M4a <- glmmTMB( d_prime ~  group +med_order  +(1|PIN),data=MKM2)
  M4b <- glmmTMB( d_prime ~  group+condition +med_order  +(1|PIN), data=MKM2)
  anova(M4a, M4b)

  plot(check_distribution(M4b))
  summary(M4b)
  ggM4 <-ggpredict(M4b, c(  "condition", "group")) %>% plot(show.title = F) + ggplot2::theme_classic()


###### ---- Results: plot single memory tests #####
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(ggpubr)

fig1 <- ggarrange(ggM1, ggM2, ggM3,ggM5, ggM4, ggM6,
                  common.legend = TRUE, legend = "bottom",
                  #labels = c("A", "B", "C", 'D'),
                  ncol = 3, nrow = 3)
fig1

## ----Results: create table for single  memory tests ####
tab_model(M1b,M2b, M4b)
tab_model( M3b, M5b, M6b  )

## ----Results: correct all Chi?? tests of individual memory parameters for multiple comparison ###################
p= c(0.003321 , 2.759e-08, 0.007641 ,1.234e-05, 0.03532 , 0.1411, 0.01532) # 
p.adjust(p, method = c( "fdr"), n = length(p))

########################################################################### #
############################################################################# #
#### 3) Results: Build latent memory change factor  #######
MKM3 <- read.table("BeckerRepantis_et_al_2021_MKM_short.csv", sep = ";", dec = ",", header=TRUE, na.strings=c("", " " , "NA", "NAN" )) 

#create difference values (Stimulant - Placebo)
MKM3$Diff_recallScanR      = (MKM3$Enhancer_recallScanR  - MKM3$Placebo_recallScanR) # memory im scanner(70 W??rter merken)
MKM3$Diff_recallAudio1r    = (MKM3$Enhancer_recallAudio1r  - MKM3$Placebo_recallAudio1r) #24h recall audio
MKM3$Diff_d_prime          = (MKM3$Enhancer_d_prime  - MKM3$Placebo_d_prime) # implicit memory (andere Ged??chtnisaufgabe, Zahlen)
MKM3$Diff_post_recallAudioR= (MKM3$Enhancer_post_recallAudioR  - MKM3$Placebo_post_recallAudioR)
MKM3$Diff_post_recallScanR = (MKM3$Enhancer_post_recallScanR  - MKM3$Placebo_post_recallScanR)
# note, Diff_post_recallScanR & the lures did not work with the model and are discarded


memoryfac1 <-  ' diffmemory  =~ 1*Diff_d_prime + Diff_recallAudio1r + Diff_post_recallAudioR + Diff_recallScanR
                                 Diff_recallAudio1r ~~ 1*Diff_recallAudio1r
                                 '
  fit <- lavaan(memoryfac1, data=MKM3
                ,auto.var=T, auto.fix.first=T ,
                auto.cov.lv.x=T, estimator = "MLR")
  summary(fit, fit.measures=TRUE)
  #modindices(fit, sort = TRUE)
  #standardizedSolution(fit)
  predict(fit)
  inspect(fit,what="std")$lambda

######################################################################## #
#### 4) Predict memory change factor from  2 hierarchical clusters Drug>Plc contrast #####
MKM4 <- read.table("BeckerRepantis_et_al_2021_MKM_FNC.csv", sep = ";", dec = ",", header=TRUE, na.strings=c("", " " , "NA", "NAN" )) 
  
MKM4$scan = as.factor(MKM4$scan)
#MKM4$med_order = as.factor(MKM4$med_order)
MKM4$PIN = as.factor(MKM4$PIN)
MKM4[MKM4$group == "Methylphenidat",]$group = "MPH"
MKM4[MKM4$group == "Modafinil",]$group = "MOD"
MKM4[MKM4$group == "Koffein",]$group = "CAF"
MKM4$group = as.factor(MKM4$group )

MKM4$Diff_Cluster1_within= MKM4$Cluster1_within_drug - MKM4$Cluster1_within_plc
MKM4$Diff_Cluster2_within= MKM4$Cluster2_within_drug - MKM4$Cluster2_within_plc

#check for collinearity
spearmanRho(x=MKM4$Diff_Cluster1_within, y= MKM4$Diff_Cluster2_within, ci = TRUE, method = "pearson", conf = 0.95, histogram = T, R=10000) 

# regression model for cluster 1
M8a <- lm( Final_MemFaktor ~ scan + med_order, data=MKM4 ) 
M8b <- lm( Final_MemFaktor ~ Diff_Cluster1_within + scan + med_order , data=MKM4 )
  plot(check_distribution(M8b))
  qqnorm(residuals(M8b)); qqline(residuals(M8b))
  anova(M8a, M8b)

# regression model for cluster 2
M8c <- lm( Final_MemFaktor ~ scan + med_order, data=MKM4 )
M8d <- lm( Final_MemFaktor ~ Diff_Cluster2_within + scan + med_order, data=MKM4 )
  plot(check_distribution(M8d))
  qqnorm(residuals(M8d)); qqline(residuals(M8d))
  anova(M8c, M8d)
  summary(M8b)

## ----Results: correct for multiple comparison for p-values of both clusters  ###################
p1= c( 0.008097 , 0.7994    )
p.adjust(p1, method = c( "fdr"), n = length(p1))

## display both regression models 
tab_model(M8b,M8d)

################################################################################# ##
### Figure 2: plot raw values  ####
fcon_c1 =  ggplot(MKM4, aes(x =Diff_Cluster1_within,y= Final_MemFaktor )) + #0.008097
  geom_point() +
  stat_smooth(method = "lm", se= TRUE) +
  labs(y="latent memory change", x= "connectivity change - cluster 1")+
  geom_jitter() +theme_classic() +  theme(text = element_text(size=15))

fcon_c2 =  ggplot(MKM4, aes(x =Diff_Cluster2_within,y= Final_MemFaktor )) + #, colour =group 
  geom_point() +
  stat_smooth(method = "lm", se= TRUE) +
  labs(y="latent memory change", x= "connectivity change - cluster 2" )+
  geom_jitter() +theme_classic() +  theme(text = element_text(size=15))

library(ggpubr)
figure2 <- ggarrange(fcon_c1,fcon_c2,
                     ncol = 2, nrow = 1,  common.legend = TRUE, legend = "right")
figure2
  # correlation values
  spearmanRho(x=MKM4$Final_MemFaktor, y= MKM4$Diff_Cluster1_within, ci = TRUE, method = "pearson", conf = 0.95, histogram = T, R=10000) 
  spearmanRho(x=MKM4$Final_MemFaktor, y= MKM4$Diff_Cluster2_within, ci = TRUE, method = "pearson", conf = 0.95, histogram = T, R=10000) 


 ### Figure S5: plot raw values  split for individual stimulants  ####
fcon_S4a =  ggplot(MKM4, aes(x =Diff_Cluster1_within,y= Final_MemFaktor, colour =group  )) + #0.008097
   geom_point() +
   stat_smooth(method = "lm", se= TRUE) +
   labs(y="latent memory change", x= "connectivity change - cluster 1")+
   geom_jitter() +theme_classic() +  theme(text = element_text(size=15))

fcon_S5b =  ggplot(MKM4, aes(x =Diff_Cluster2_within,y= Final_MemFaktor , colour =group )) + #, colour =group 
  geom_point() +
  stat_smooth(method = "lm", se= TRUE) +
  labs(y="latent memory change", x= "connectivity change - cluster 2" )+
  geom_jitter() +theme_classic() +  theme(text = element_text(size=15))

#library(ggpubr)
figureS5 <- ggarrange(fcon_S4a,fcon_S4b,
           ncol = 2, nrow = 1,  common.legend = TRUE, legend = "right")
figureS5



