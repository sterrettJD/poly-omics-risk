# base_dir = "C:\\Users\\rgarr\\Documents\\poly-omics-risk\\"
#setwd("/Users/johnsterrett/Research-Projects/Team-rotation/poly-omics-scores/viromics/")
base_dir = "/Users/chris/Documents/GRADSCHOOL/PolyOmicsRotation/poly-omics-risk"

list.files("./viromics/")

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
`%ni%` <- Negate(`%in%`)


## import test data metadata and extract patient ID, external ID and diagnosis--####

setwd(base_dir)
metadata <- fread("./testing_metadata.txt", header=T)
metadata <- as.data.frame(metadata[,c("External ID","Participant ID","diagnosis")])
colnames(metadata) <- c("External_ID","Participant_ID", "diagnosis")
metadata$diagnosis <- as.character(metadata$diagnosis)
metadata$diagnosis[metadata$diagnosis == "UC"] <- 1
metadata$diagnosis[metadata$diagnosis == "CD"] <- 1
metadata$diagnosis[metadata$diagnosis == "nonIBD"] <- 0
metadata$diagnosis <- as.numeric(metadata$diagnosis)

metadata


## import MixRF predictions  --####
# will duplicate, one for lasso one for MixRF

# viromics
# setwd(glue("{base_dir}viromics"))
viromic_MixRF <- fread("./viromics/viromics_MixRF_scores.txt", header=T)
viromic_MixRF <- as.data.frame(viromic_MixRF[,c("External_ID","predicted")])
viromic_MixRF <- viromic_MixRF %>% rename("virome_pred" ="predicted")
viromic_MixRF

# metagenomics
# setwd(glue("{base_dir}metagenomics"))
metagenomic_MixRF <- fread("./metagenomics/metagenomics_MixRF_scores.txt", header=T)
metagenomic_MixRF <- as.data.frame(metagenomic_MixRF[,c("External_ID","predicted")])
metagenomic_MixRF <- metagenomic_MixRF %>% rename("metagen_pred" ="predicted")
metagenomic_MixRF

# metatranscriptomics
# setwd(glue("{base_dir}metatranscriptomics"))
metatranscriptomic_MixRF <- fread("./metatranscriptomics/metatranscriptomics_MixRF_scores.txt", header=T)
metatranscriptomic_MixRF <- as.data.frame(metatranscriptomic_MixRF[,c("External_ID","predicted")])
metatranscriptomic_MixRF <- metatranscriptomic_MixRF %>% rename("metatrans_pred" ="predicted")
metatranscriptomic_MixRF

# # metabolomics
# setwd(glue("{base_dir}metabolomics"))
metabolomic_MixRF <- fread("./metabolomics/metabolomics_MixRF_scores.txt", header=T)
metabolomic_MixRF <- as.data.frame(metabolomic_MixRF[,c("External_ID","predicted")])
metabolomic_MixRF <- metabolomic_MixRF %>% rename("metabol_pred" ="predicted")
metabolomic_MixRF


## merge datasets by external ID  --#################################################################################
# get dataset with cols: extID,Diagnosis,Score1,Score2,Score3,Score4

# df_list <- list(metadata,viromic_MixRF,metagenomic_MixRF,metatranscriptomic_MixRF)
df_list <- list(metadata,viromic_MixRF,metagenomic_MixRF,metatranscriptomic_MixRF,metabolomic_MixRF)
rf_scores <- df_list %>% reduce(full_join, by="External_ID")
rf_scores <- as.data.frame(rf_scores[complete.cases(rf_scores),])
rf_scores 


## average individual scores by participant ID --####

avg_par_scores <- rf_scores %>% 
  group_by(Participant_ID) %>%
  summarize(diagnosis = mean(diagnosis), 
            virome_pred = mean(virome_pred), 
            metagen_pred = mean(metagen_pred), 
            metatrans_pred = mean(metatrans_pred),
            metabol_pred = mean(metabol_pred))
avg_par_scores <- as.data.frame(avg_par_scores)
rownames(avg_par_scores) <- avg_par_scores$`Participant_ID`
avg_par_scores$`Participant_ID` <- NULL

view(avg_par_scores)

## combined regression --#####



combomod <- glm(as.formula(paste0("diagnosis ~ .")), data = avg_par_scores, family = "binomial")

combomod_sum <- summary(combomod)
combomod_sum

# #---make a regression tree---#########################################################################################################################################################################
# require(tree)
mytree <- tree(diagnosis ~ ., data = avg_par_scores)
plot(mytree)
text(mytree, pretty = 0, cex = .8)

