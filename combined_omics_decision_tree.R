# base_dir = "C:\\Users\\rgarr\\Documents\\poly-omics-risk\\"
# setwd("/Users/johnsterrett/Research-Projects/Team-rotation/poly-omics-scores/viromics/")
# base_dir = "/Users/chris/Documents/GRADSCHOOL/PolyOmicsRotation/poly-omics-risk"
base_dir = "C:/Users/chris/OneDrive/Documents/poly-omics-risk"

# list.files("./viromics/")

# install.packages("data.table")
# install.packages("ggplot2")
# install.packages("R.utils")
# install.packages("tidyverse")
# install.packages("UpSetR")
# install.packages("cowplot")
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("biomformat")
# install.packages("bestglm")
# install.packages("tree")
# install.packages("compositions")
# install.packages("lme4")
# install.packages("glmmLasso")
# install.packages("pROC")
# install.packages("ggpubr")
# install.packages("MixRF")
# install.packages("caret")
# install.packages("randomForest")
# install.packages("glue") 
# install.packages("sjPlot")
library(data.table)
library(ggplot2)
library(R.utils)
library(tidyverse)
library(UpSetR)
library(cowplot)
library(biomformat)
library(compositions)

library(bestglm)
library(MASS)
library(tree)
library(lme4)
library(glmmLasso)
library(pROC)
library(ggpubr)
library(MixRF)
library(caret)
library(randomForest)
library(glue)
library(fmsb)
library(sjPlot)
library(sjmisc)
library(sjlabelled)
library(rcompanion)
`%ni%` <- Negate(`%in%`)


## import test data metadata and extract patient ID, external ID and diagnosis--####

setwd(base_dir)
metadata <- fread("./testing_metadata.txt", header=T)
metadata <- as.data.frame(metadata[,c("External ID","Participant ID","diagnosis", "consent_age", "race","sex", "site_name")])
colnames(metadata) <- c("External_ID","Participant_ID", "diagnosis","consent_age", "race", "sex","site_name")
metadata$diagnosis <- as.character(metadata$diagnosis)
metadata$diagnosis[metadata$diagnosis == "UC"] <- 1
metadata$diagnosis[metadata$diagnosis == "CD"] <- 1
metadata$diagnosis[metadata$diagnosis == "nonIBD"] <- 0
metadata$diagnosis <- as.numeric(metadata$diagnosis)
metadata$External_ID <- as.character(metadata$External_ID)
metadata$Participant_ID <- as.character(metadata$Participant_ID)
metadata$sex <- as.factor(metadata$sex)
metadata$site_name <- as.factor(metadata$site_name)
metadata$consent_age <- as.numeric(metadata$consent_age)
metadata$race <- as.factor(metadata$race)
metadata


## import MixRF predictions  --####
# will duplicate, one for lasso one for MixRF

# viromics
# setwd(glue("{base_dir}viromics"))
viromic_glmer <- fread("./viromics/viromics_features_scores.txt", header=T)
viromic_glmer <- as.data.frame(viromic_glmer[,c("External_ID","predicted")])
viromic_glmer <- viromic_glmer %>% rename("virome_pred" ="predicted")
viromic_glmer

# metagenomics
# setwd(glue("{base_dir}metagenomics"))
metagenomic_glmer <- fread("./metagenomics/metagenomics_features_scores.txt", header=T)
metagenomic_glmer <- as.data.frame(metagenomic_glmer[,c("External_ID","predicted")])
metagenomic_glmer <- metagenomic_glmer %>% rename("metagen_pred" = "predicted")
metagenomic_glmer

# metatranscriptomics
# setwd(glue("{base_dir}metatranscriptomics"))
metatranscriptomic_glmer <- fread("./metatranscriptomics/metatranscriptomics_features_scores.txt", header=T)
metatranscriptomic_glmer <- as.data.frame(metatranscriptomic_glmer[,c("External_ID","predicted")])
metatranscriptomic_glmer <- metatranscriptomic_glmer %>% rename("metatrans_pred" ="predicted")
metatranscriptomic_glmer

# # metabolomics
# setwd(glue("{base_dir}metabolomics"))
metabolomic_glmer <- fread("./metabolomics/metabolomics_features_scores.txt", header=T)
metabolomic_glmer <- as.data.frame(metabolomic_glmer[,c("External_ID","predicted")])
metabolomic_glmer <- metabolomic_glmer %>% rename("metabol_pred" ="predicted")
metabolomic_glmer


## merge datasets by external ID  --#################################################################################
# get dataset with cols: extID,Diagnosis,Score1,Score2,Score3,Score4

# df_list <- list(metadata,viromic_MixRF,metagenomic_MixRF,metatranscriptomic_MixRF)
df_list <- list(metadata,viromic_glmer,metagenomic_glmer,metatranscriptomic_glmer,metabolomic_glmer)
rf_scores <- df_list %>% reduce(full_join, by="External_ID")
rf_scores <- as.data.frame(rf_scores[complete.cases(rf_scores),])
rf_scores 



## average individual scores by participant ID --########################################################################

avg_par_scores <- rf_scores %>% 
  group_by(Participant_ID) %>%
  summarize( virome_pred = mean(virome_pred), 
            metagen_pred = mean(metagen_pred), 
            metatrans_pred = mean(metatrans_pred),
            metabol_pred = mean(metabol_pred))
avg_par_scores <- as.data.frame(avg_par_scores)


par_meta <- subset(metadata, select = -(External_ID)) 
par_meta <- par_meta %>% distinct()
par_meta

df_for_model <- merge(par_meta, avg_par_scores, by ='Participant_ID')
df_for_model


rownames(df_for_model) <- df_for_model$`Participant_ID`
df_for_model$`Participant_ID` <- NULL

# view(df_for_model)
colnames(df_for_model) <- c("diagnosis","consent_age","race","sex","site_name","VRM","MGN","MTS","MBL")


## combined regression --#####

combomod <- glm(as.formula(paste0("diagnosis ~ .")), data = df_for_model, family = "binomial")


combomod_sum <- summary(combomod)
combomod_sum


combomod_nocovar <- glm(as.formula(paste0("diagnosis ~ VRM + MGN + MTS + MBL")), data = df_for_model, family = "binomial")

combomod_nocovar_sum <- summary(combomod_nocovar)
combomod_nocovar_sum

combomod_somecovar <- glm(as.formula(paste0("diagnosis ~ VRM + MGN + MTS + MBL + consent_age + sex")), data = df_for_model, family = "binomial")

combomod_somecovar_sum <- summary(combomod_somecovar)
combomod_somecovar_sum
plot_model(combomod_somecovar, vline.color = 'gray') + theme_bw() + ggtitle("")
ggsave("combined_model_forest_OR.png", width=5, height=3, units="in", dpi=320)

# leave one out
loomod <- glm(as.formula(paste0("diagnosis ~ MGN + MTS + MBL + consent_age + sex")), data = df_for_model, family = "binomial")
plot_model(loomod, vline.color = 'gray') + theme_bw()
loomod <- glm(as.formula(paste0("diagnosis ~ VRM + MTS + MBL + consent_age + sex")), data = df_for_model, family = "binomial")
plot_model(loomod, vline.color = 'gray') + theme_bw()
loomod <- glm(as.formula(paste0("diagnosis ~ VRM + MGN + MBL + consent_age + sex")), data = df_for_model, family = "binomial")
plot_model(loomod, vline.color = 'gray') + theme_bw()
loomod <- glm(as.formula(paste0("diagnosis ~ VRM + MGN + MTS + consent_age + sex")), data = df_for_model, family = "binomial")
plot_model(loomod, vline.color = 'gray') + theme_bw()

combomod_only_somecovar <- glm(as.formula(paste0("diagnosis ~ consent_age + sex")), data = df_for_model, family = "binomial")
# Nagelkerke (Cragg and Uhler)        0.285719 (using only age sex and race)
nagelkerke(combomod_only_somecovar, null = NULL, restrictNobs = FALSE)

combomod_only_somecovar <- summary(combomod_only_somecovar)
combomod_only_somecovar
## nagelkerke R #####


NagelkerkeR2(combomod_nocovar) #0.3657229 (using only the scores)
NagelkerkeR2(combomod)
NagelkerkeR2(combomod_somecovar)
NagelkerkeR2(combomod_only_somecovar)

library(rcompanion)
nagelkerke(combomod_only_somecovar, null = NULL, restrictNobs = FALSE)
# Nagelkerke (Cragg and Uhler)         0.464423 (using scores + age + sex + race)
nagelkerke(combomod_somecovar, null = NULL, restrictNobs = FALSE)

library(corrplot)
mtcor <- as.matrix(df_for_model[,c("VRM","MGN","MTS","MBL")])
M = cor(mtcor)
testRes = cor.mtest(mtcor, conf.level = 0.95)

corrplot.mixed(M)
cp <- corrplot(M, addCoef.col = 'black', tl.pos = 'd',
         cl.pos = 'n', col = COL2('PiYG'), type = "lower")
ggsave(plot = print(cp), filename = "correlation_of_scores.png", width=3, height=3, units="in", dpi=320)
########## make auc plot

auc_df = df_for_model

auc_df$diagnosis <- as.factor(auc_df$diagnosis)
auc_df$predicted <- as.numeric(scale(auc_df$VRM + auc_df$MGN + auc_df$MTS + auc_df$MBL))
PredPlot <- ggplot(data = auc_df, aes(x = diagnosis, y = predicted))+
  scale_fill_viridis_d( option = "D")+
  theme_dark()+
  geom_violin(fill = "gray70",alpha=0.4, position = position_dodge(width = .5),size=1,color="gray22",width=.5,lwd=.2) +
  geom_boxplot(fill = "gray95",notch = F, shape=21, outlier.size = -1, color="gray32",lwd=.5, alpha = .75)+
  theme(plot.title = element_text(hjust = 0.5))+
  # theme(axis.title.x = element_text(size=14))+
  theme(axis.text.x = element_text(colour = "black"))+
  theme(axis.text.y = element_text(colour = "black"))+
  theme( axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"))+
  theme(panel.grid.major.x = element_blank())+
  theme(panel.background = element_rect(fill = 'white'), panel.grid = element_line(color='gray80'))+
  # labs(title = addToTitle)+
  ylab("Score (VRM + MGN + MTS + MBL)")+
  xlab("Actual Diagnosis")


caseControlGLM <- glm(as.formula(auc_df$diagnosis ~ auc_df$predicted), data = auc_df, family = "binomial", na.action = na.omit)
# caseControlGLM <- glm(as.formula(auc_df$diagnosis ~ auc_df$predicted + auc_df$consent_age + auc_df$sex), data = auc_df, family = "binomial", na.action = na.omit)


predpr <- predict(caseControlGLM, auc_df, allow.new.levels = T, type = c("response"))
caseControlroccurve <- pROC::roc(auc_df$diagnosis ~ predpr, quiet=T, plot=T)
caseControlroccurveCI <- pROC::roc(auc_df$diagnosis ~ predpr, ci=T, quiet=T)
nr2 <- NagelkerkeR2(caseControlGLM)$R2
print(paste0("Nagelkerke's R2: ", nr2))
# caseControlplot <- plot(caseControlroccurve, main=paste("Case vs Control AUC =", round(caseControlroccurve$auc, 3)))

caseControlp <- formatC(coef(summary(caseControlGLM))[,4][2], format = "e", digits = 0)
caseControlOR <- exp(cbind("Odds ratio" = coef(caseControlGLM), confint.default(caseControlGLM, level = 0.95)))
if (as.numeric(strsplit(caseControlp,"")[[1]][1]) == 1){
  caseControlpthresh <- paste0("OR [95% CI] = ", format(round(as.numeric(caseControlOR[2,1]), 2), nsmall = 2), " [", format(round(as.numeric(caseControlOR[2,2]), 2),nsmal = 2), ", ", format(round(as.numeric(caseControlOR[2,3]), 2),nsmal = 2),"], ","p < ", as.numeric(caseControlp)/as.numeric(strsplit(caseControlp,"")[[1]][1]))
}else{
  caseControlpthresh <- paste0("OR [95% CI] = ", format(round(as.numeric(caseControlOR[2,1]), 2), nsmall = 2), " [", format(round(as.numeric(caseControlOR[2,2]), 2),nsmal = 2), ", ", format(round(as.numeric(caseControlOR[2,3]), 2),nsmal = 2),"], ","p < ", 10*as.numeric(caseControlp)/as.numeric(strsplit(caseControlp,"")[[1]][1]))
}  
if (as.numeric(caseControlp) >= 0.001){
  caseControlpthresh <- paste0("OR [95% CI] = ", format(round(as.numeric(caseControlOR[2,1]), 2), nsmall = 2), " [", format(round(as.numeric(caseControlOR[2,2]), 2),nsmal = 2), ", ", format(round(as.numeric(caseControlOR[2,3]), 2),nsmal = 2),"], ","p = ", formatC(as.numeric(caseControlp), format = "g"))
}
caseControlAUCsummary <- paste0("AUC [95% CI] = ", format(round(as.numeric(caseControlroccurve$auc), 2), nsmall = 2), " [", format(round(as.numeric(caseControlroccurveCI$ci[1]), 2),nsmal = 2), ", ", format(round(as.numeric(caseControlroccurveCI$ci[3]), 2),nsmal = 2),"]")
caseControlORsummary <- caseControlpthresh

label_disty = .13
minPRS <- min(auc_df$predicted)-(abs(max(auc_df$predicted)-min(auc_df$predicted)))*(label_disty*1.5*1.02)
maxPRS <- max(auc_df$predicted) + (abs(max(auc_df$predicted)-min(auc_df$predicted)))*(label_disty*1.5*1.02)
botLabLoc <- min(auc_df$predicted)-(abs(max(auc_df$predicted)-min(auc_df$predicted)))*(label_disty*1.5)
topLabLoc <- max(auc_df$predicted) + (abs(max(auc_df$predicted)-min(auc_df$predicted)))*(label_disty*0.21)

plot_df <- as.data.frame(cbind(caseControlGLM$fitted.values, auc_df$diagnosis))

colnames(plot_df) <- c("predicted", "diagnosis")
plot_df$diagnosis <- as.factor(plot_df$diagnosis)


PredPlot <- PredPlot +
  geom_signif(textsize = 2.25, comparisons = list(c("0", "1")), annotations=caseControlAUCsummary, color="black", y_position = topLabLoc,tip_length=.03)+
  geom_signif(textsize = 2.25, comparisons = list(c("0", "1")), annotations=caseControlpthresh, color="black", y_position = botLabLoc,tip_length=-.03) +
  scale_y_continuous(breaks = seq(-100,100, by=2), limits = c(minPRS,maxPRS))

print(caseControlAUCsummary)
print(caseControlpthresh)
PredPlot
ggsave("combinde_scores_auc.png", width=3, height=3, units="in", dpi=320)
stop()



# #---make a regression tree---#########################################################################################################################################################################
# require(tree)
mytree <- tree(diagnosis ~ VRM + MGN + MTS + MBL, data = df_for_model, control = tree.control(mindev =0, minsize=3, nobs=29) )
plot(mytree)
text(mytree, pretty = 1, cex = .8,  digits = 1)


tab_model(combomod_nocovar)
plot_model(combomod_nocovar, vline.color = 'gray') + theme_bw()


## looking at variance for participant -################
# setwd(base_dir)
# metadata2 <- fread("./testing_metadata.txt", header=T)
# metadata2 <- as.data.frame(metadata2[,c("External ID","Participant ID","diagnosis", "consent_age", "race","sex", "site_name","Antibiotics")])
# colnames(metadata2) <- c("External_ID","Participant_ID", "diagnosis","consent_age", "race", "sex","site_name", "Antibiotics")
# metadata2$diagnosis <- as.character(metadata2$diagnosis)
# metadata2$diagnosis[metadata2$diagnosis == "UC"] <- 1
# metadata2$diagnosis[metadata2$diagnosis == "CD"] <- 1
# metadata2$diagnosis[metadata2$diagnosis == "nonIBD"] <- 0
# metadata2$diagnosis <- as.numeric(metadata2$diagnosis)
# metadata2$External_ID <- as.character(metadata2$External_ID)
# metadata2$Participant_ID <- as.character(metadata2$Participant_ID)
# metadata2$sex <- as.factor(metadata2$sex)
# metadata2$site_name <- as.factor(metadata2$site_name)
# metadata2$consent_age <- as.numeric(metadata2$consent_age)
# metadata2$race <- as.factor(metadata2$race)
# metadata2$Antibiotics <- as.factor(metadata2$Antibiotics)
# metadata2
# 
# df_list <- list(metadata2,viromic_MixRF,metagenomic_MixRF,metatranscriptomic_MixRF,metabolomic_MixRF)
# rf_scores <- df_list %>% reduce(full_join, by="External_ID")
# rf_scores <- as.data.frame(rf_scores[complete.cases(rf_scores),])
# rf_scores 
# 
# par_variance = rf_scores %>%
#   subset(select= -c(External_ID, consent_age, race, sex, site_name))
# 
# hist(group_by) 

par_variance <- rf_scores %>%
  group_by(Participant_ID) %>%
  summarize( virome_pred_avg = mean(virome_pred),
             virome_var = var(virome_pred),
             virome_range = max(virome_pred)-min(virome_pred),
             #virome_min = min(virome_pred),
             #virome_max = max(virome_pred),
             metagen_pred_avg = mean(metagen_pred), 
             metagen_range = max(metagen_pred)-min(metagen_pred),
             metatrans_pred_avg = mean(metatrans_pred),
             metatrans_var = var(metatrans_pred),
             metatrans_range = max(metatrans_pred)-min(metatrans_pred),
             metabol_pred_avg = mean(metabol_pred),
             metabol_var = var(metabol_pred),
             metabol_range = max(metabol_pred)-min(metabol_pred)
             )
  


par_variance  


hist(par_variance$virome_var)
hist(par_variance$metagen_var)
hist(par_variance$metatrans_var)
hist(par_variance$metabol_var)

hist(par_variance$virome_range)
hist(par_variance$metagen_range)
hist(par_variance$metatrans_range)
hist(par_variance$metabol_range)


vir_by_par <- rf_scores %>%
  transmute(Participant_ID=Participant_ID, virome_pred=virome_pred)

vir_by_par$virome_pred <- as.numeric(vir_by_par$virome_pred)
vir_by_par$Participant_ID <- as.factor(vir_by_par$Participant_ID)



# vir_by_par <- vir_by_par %>%
#   group_by(Participant_ID) %>%
#   hist(virome_pred)


vir_by_par <- reshape(
  data = vir_by_par,
  direction = "wide",
  idvar = "Participant_ID",
  
)





for (par in rf_scores$Participant_ID) {
  hist(rf_scores$virome_pred)
}
