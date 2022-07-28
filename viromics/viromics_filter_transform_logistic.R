 # setwd("C:\\Users\\rgarr\\Documents\\poly-omics-risk")
#setwd("/Users/johnsterrett/Research-Projects/Team-rotation/poly-omics-scores/viromics/")
setwd("/Users/chris/Documents/GRADSCHOOL/PolyOmicsRotation/poly-omics-risk/viromics/")

list.files()

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
# install.packages("fmsb")
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
library(fmsb)
`%ni%` <- Negate(`%in%`)

# ## metadata --###################################
# metadata <- fread("https://ibdmdb.org/tunnel/products/HMP2/Metadata/hmp2_metadata.csv", header=T, stringsAsFactors=T)
# str(metadata)
# 
# 
# # make upset plot for sample overlap by datatype
# listInput <- list(
#   # biopsy_16S = c(as.character(subset(metadata, data_type == "biopsy_16S")[["External ID"]])),
#   # host_genome = c(as.character(subset(metadata, data_type == "host_genome")[["External ID"]])),
#   # host_transcriptomics  = c(as.character(subset(metadata, data_type == "host_transcriptomics")[["External ID"]])),
#   # methylome = c(as.character(subset(metadata, data_type == "methylome")[["External ID"]])),
#   # proteomics = c(as.character(subset(metadata, data_type == "proteomics")[["External ID"]])),
#   # serology  = c(as.character(subset(metadata, data_type == "serology")[["External ID"]])),
#   # stool_16S = c(as.character(subset(metadata, data_type == "stool_16S")[["External ID"]])),
#   metabolomics  = c(as.character(subset(metadata, data_type == "metabolomics")[["External ID"]])), 
#   metagenomics = c(as.character(subset(metadata, data_type == "metagenomics")[["External ID"]])),  
#   metatranscriptomics = c(as.character(subset(metadata, data_type == "metatranscriptomics")[["External ID"]])), 
#   viromics  = c(as.character(subset(metadata, data_type == "viromics")[["External ID"]]))
# )
# upset(fromList(listInput), order.by = "freq")
# 
# # make a list of the ids that are in the intersection of metabolomics, metagenomics, viromics, metatranscriptomics
# i1 <- intersect(c(as.character(subset(metadata, data_type == "metabolomics")[["External ID"]])), 
#                 c(as.character(subset(metadata, data_type == "metagenomics")[["External ID"]]))
# )
# 
# i2 <- intersect(i1, c(as.character(subset(metadata, data_type == "viromics")[["External ID"]])))
# i3 <- intersect(i2, c(as.character(subset(metadata, data_type == "metatranscriptomics")[["External ID"]])))
# length(i3)
# idlist <- unique(i3)
# # idlist <- as.data.frame(idlist); colnames(idlist) <- c("idlist")
# # head(idlist)
# # idlistsep <- idlist %>% separate(idlist,into=c("idlist","junk"),convert=TRUE,sep="_profi")
# # idlistsep <- idlistsep$idlist
# 
# # subset to jsut one datatype to explore diagnosis counts
# metadata <- subset(metadata, data_type == "viromics")
# dim(metadata)
# # metadata <- subset(metadata, `External ID` %in% idlist)
# 
# 
# # fill missinig enries with NA and view the columns that are mostly non-NA
# metadata[metadata==""] <- NA
# isna <- sapply(metadata, function(x) sum(is.na(x)))
# isna[isna < 100]
# summary(metadata$diagnosis)
# metadata1 <- as.data.frame(metadata[,c("External ID","diagnosis","Participant ID","site_name","consent_age","sex")])
# metadata1$diagnosis <- as.character(metadata1$diagnosis)
# metadata1$diagnosis[metadata1$diagnosis == "UC"] <- 1
# metadata1$diagnosis[metadata1$diagnosis == "CD"] <- 1
# metadata1$diagnosis[metadata1$diagnosis == "nonIBD"] <- 0
# metadata1$diagnosis <- as.numeric(metadata1$diagnosis)
# metadata1$'External ID' <- as.character(metadata1$'External ID')
# metadata1$'Participant ID' <- as.character(metadata1$'Participant ID')
# metadata1$sex <- as.factor(metadata1$sex)
# metadata1$site_name <- as.factor(metadata1$site_name)
# metadata1$consent_age <- as.numeric(metadata1$consent_age)

## Viromics taxonomic_profiles --###################################
# #MetaPhlAn2
# fname <- fread("https://ibdmdb.org/tunnel/products/HMP2/Viromics/1732/taxonomic_profiles.tsv.gz", header=T)
# fname <- fread("taxonomic_profiles_3.tsv.gz", header=T)
# fname <- as.data.frame(fname)

## load VirMap --#######################
fname <- fread("https://ibdmdb.org/tunnel/products/HMP2/Viromics/1732/virome_virmap_analysis.tsv.gz", header=T)
#fname <- fread("virome_virmap_analysis.tsv.gz", header=T)
fname <- as.data.frame(fname)

df_3 <- fname
#View(df_3)
print("data dimensions:")
print(dim(df_3))
print("data glimpse:")
str(df_3[1:5,1:4])
rownames(df_3) <- df_3$`Virus`
df_3$`Virus` <- NULL
#View(as.data.frame(rownames(df_3)))
# separate dataframe out by rows which were annotated vs rows which weren't
commentRow <- grep("# ", rownames(df_3), invert = F)
mygrep <- grep("UNKNOWN", rownames(df_3), invert = F)
grep_species <- grep("species=", rownames(df_3), invert = F)
mygrepni <- which(1:length(rownames(df_3)) %ni% mygrep)
if(length(commentRow) > 0){
  mygrepni <- mygrepni[which(mygrepni != commentRow)]
}



# Preprocessing, make data compositional for species --###################################
grouped_df_3 <- df_3[mygrepni,]
rownames(grouped_df_3) <- rownames(df_3)[mygrepni]
grouped_df_3 <- df_3[grep_species,]
rownames(grouped_df_3) <- rownames(df_3)[grep_species]
num_grouped_df_3 <- mutate_all(grouped_df_3, function(x) as.numeric(as.character(x)))
#View(num_grouped_df_3)
#View(as.data.frame(rownames(num_grouped_df_3)))


# is it compositional?
summary(colSums(num_grouped_df_3))
hist(colSums(num_grouped_df_3))


# percent missing filter
percMissing <- rowSums(num_grouped_df_3==0)/ncol(num_grouped_df_3)*100
hist(percMissing, xlim = c(80,100),  plot = TRUE)
length(percMissing)
length(which(percMissing>5))
length(which(percMissing>1))

num_grouped_df_3 <- num_grouped_df_3[which(percMissing<95),]

summary(colSums(num_grouped_df_3))
hist(colSums(num_grouped_df_3))
# view(num_grouped_df_3)
## add epsilon to all entries to set up for center log transform --##############
pepsi <- 1E-06
num_grouped_df_3 <- num_grouped_df_3 + pepsi
#View(num_grouped_df_3)


## CLT transform --###################################
# transform by center log transform across rows (features)
clr_num_grouped_df_3 <- compositions::clr(num_grouped_df_3,)
# clr_num_grouped_df_3 <- as.data.frame(apply(num_grouped_df_3, 1, clr))

alldf_pretransform <- as.data.frame(as.numeric(array(as.matrix(num_grouped_df_3)))); colnames(alldf_pretransform) <- c("alldf_pretransform")
hist(alldf_pretransform$alldf_pretransform)
alldf_posttransform <- as.data.frame(as.numeric(array(as.matrix(clr_num_grouped_df_3)))); colnames(alldf_posttransform) <- c("alldf_posttransform")
hist(alldf_posttransform$alldf_posttransform)

## rename the columns for the merge with diagnosis --##############
namesies <- as.data.frame(colnames(clr_num_grouped_df_3)); colnames(namesies) <- c("namesies")
head(namesies)
# namesiessep <- namesies %>% separate(namesies,into=c("namesies","junk"),convert=TRUE,sep="_profi")
namesiessep <- namesies
colnames(clr_num_grouped_df_3) <- namesiessep$namesies

# intersecty <- intersect(c(as.character(namesiessep$namesies)), 
#                         c(as.character(metadata1$`External ID`))
# )

#View(intersecty)
dim(namesiessep)
# dim(metadata1)
# length(intersecty)



## transpose transformed data --##############
transp_clr_num_grouped_df_3 <- as.data.frame(t(clr_num_grouped_df_3))
transp_clr_num_grouped_df_3$`External ID` <- rownames(transp_clr_num_grouped_df_3)

## read in and Merge with metadata --###################################
#read in training set and merge with viromics data
train_meta= fread("../training_metadata.txt", sep='\t', header=FALSE)
train_meta=as.data.frame(train_meta)
colnames(train_meta)= c("External ID", "Participant ID","race","data_type", "consent_age", "site_name", "diagnosis", "sex", "Antibiotics")
train_meta=train_meta %>% filter(data_type == 'viromics')

train_meta$diagnosis <- as.character(train_meta$diagnosis)
train_meta$diagnosis[train_meta$diagnosis == "UC"] <- 1
train_meta$diagnosis[train_meta$diagnosis == "CD"] <- 1
train_meta$diagnosis[train_meta$diagnosis == "nonIBD"] <- 0
train_meta$diagnosis <- as.numeric(train_meta$diagnosis)
train_meta$'External ID' <- as.character(train_meta$'External ID')
train_meta$'Participant ID' <- as.character(train_meta$'Participant ID')
train_meta$sex <- as.factor(train_meta$sex)
train_meta$site_name <- as.factor(train_meta$site_name)
train_meta$consent_age <- as.numeric(train_meta$consent_age)
train_meta$race <- as.factor(train_meta$race)

# view(train_meta)


#read in testing metadata and merge with viromics
test_meta= fread("../testing_metadata.txt", sep='\t', header=FALSE)
test_meta=as.data.frame(test_meta)
colnames(test_meta)= c("External ID", "Participant ID","race","data_type", "consent_age", "site_name", "diagnosis", "sex", "Antibiotics")
test_meta=test_meta %>% filter(data_type == 'viromics')

test_meta$diagnosis <- as.character(test_meta$diagnosis)
test_meta$diagnosis[test_meta$diagnosis == "UC"] <- 1
test_meta$diagnosis[test_meta$diagnosis == "CD"] <- 1
test_meta$diagnosis[test_meta$diagnosis == "nonIBD"] <- 0
test_meta$diagnosis <- as.numeric(test_meta$diagnosis)
test_meta$'External ID' <- as.character(test_meta$'External ID')
test_meta$'Participant ID' <- as.character(test_meta$'Participant ID')
test_meta$sex <- as.factor(test_meta$sex)
test_meta$site_name <- as.factor(test_meta$site_name)
test_meta$consent_age <- as.numeric(test_meta$consent_age)
test_meta$race <- as.factor(test_meta$race)


# Reclassify the one individual in testing that has a unique race class:
summary(test_meta$race)
summary(train_meta$race)
test_meta$race <- as.character(test_meta$race)
# clump "American Indian or Alaska Native" into "Other"
test_meta$race[which(test_meta$race == "American Indian or Alaska Native")] <- "Other"
test_meta$race <- as.factor(test_meta$race)
summary(test_meta$race)
summary(train_meta$race)


## combine data with metadata and reformat --#####
mergey <- merge(transp_clr_num_grouped_df_3, train_meta, by = "External ID")
mergeytest <- merge(transp_clr_num_grouped_df_3, test_meta, by = "External ID")
# make the ids the rownames for each dataframe, and then remove that column
rownames(mergeytest) <- mergeytest$`External ID`
mergeytest$`External ID` <- NULL
rownames(mergey) <- mergey$`External ID`
mergey$`External ID` <- NULL

## remove any columns that have no/little variation between samples --####
mergeytest_colsd <- apply(mergeytest[complete.cases(mergeytest), 
                                     grep("species=", names(mergeytest), invert = F)], 2, sd, na.rm=T)
mergey_colsd <- apply(mergey[complete.cases(mergey),
                             grep("species=", names(mergey), invert = F)], 2, sd, na.rm=T)
# mergeytest_colsd <- apply(mergeytest, 2, sd, na.rm=T)
# mergey_colsd <- apply(mergey, 2, sd, na.rm=T)

# qthresh <- quantile(colsd, 0.05, na.rm=T)
hist(as.numeric(mergeytest_colsd))
hist(as.numeric(mergey_colsd))
toKeep <- c(intersect(c(names(mergeytest_colsd)[which(mergeytest_colsd > 0.1)]),
                      c(names(mergey_colsd)[which(mergey_colsd > 0.1)])),
            c("External ID", "Participant ID", "site_name", "diagnosis", 
              "consent_age", "sex", "race", "Antibiotics")
)
# toKeep <- intersect(c(names(mergeytest_colsd)[which(mergeytest_colsd > 0)]), c(names(mergey_colsd)[which(mergey_colsd > 0)]))

length(toKeep)
length(mergey_colsd)
mergeytest <- mergeytest[,which(names(mergeytest) %in% toKeep)]
mergey <- mergey[,which(names(mergey) %in% toKeep)]

## rename the column names with just the species and remove spaces--####
cn <- as.data.frame(colnames(mergey)); colnames(cn) <- c("cn")
cn <- cn %>% separate(cn,into=c("junk","species_taxid"),convert=TRUE,sep="species=")
species_taxid = cn$species_taxid
cn <- cn %>% separate(species_taxid,into=c("species","therest"),convert=TRUE,sep=";")
cn <- cn$species
cn[is.na(cn)] <- c("Participant_ID","race","consent_age","site_name","diagnosis","sex","Antibiotics")
#cn[is.na(cn)] <- c("Participant_ID", "site_name", "diagnosis", "consent_age", "sex", "race", "Antibiotics")
cn <- gsub(" ","_",cn)
cn <- gsub("-","_",cn)
cn <- gsub(pattern = "\\[", replacement = "", cn)
cn <- gsub(pattern = "\\]", replacement = "", cn)

head(cn,20)
tail(cn,20)

colnames(mergeytest) <- cn
colnames(mergey) <- cn


## run LASSO regression for all features --###################################

## linear mixed model
mergey$diagnosis <- as.integer(as.character(mergey$diagnosis))
# mergey$diagnosis <- as.numeric(mergey$diagnosis)
mergey$site_name <- as.factor(mergey$site_name)
mergey$consent_age <- as.numeric(mergey$consent_age)
mergey$`Participant_ID` <- as.factor(mergey$`Participant_ID`)
mergey$sex <- as.factor(mergey$sex)
mergey$race <- as.factor(mergey$race)
mergey$Antibiotics <- as.factor(mergey$Antibiotics)


# str(mergey[,1:30])
# str(mergey[,(ncol(mergey)-30):ncol(mergey)])

traindf <- as.data.frame(mergey[complete.cases(mergey),])

traindf$diagnosis <- as.integer(as.character(traindf$diagnosis))
# traindf$diagnosis <- as.numeric(traindf$diagnosis)
traindf$site_name <- as.factor(traindf$site_name)
traindf$consent_age <- as.numeric(traindf$consent_age)
traindf$`Participant_ID` <- as.factor(traindf$`Participant_ID`)
traindf$sex <- as.factor(traindf$sex)
traindf$race <- as.factor(traindf$race)
traindf$Antibiotics <- as.factor(traindf$Antibiotics)


dim(traindf)
dim(mergey)
varlist <- cn[which(cn %ni% c("diagnosis", "Participant_ID"))]

for(i in 1:length(varlist)){
  if(varlist[i] %in% c("site_name", "sex", "race", "Antibiotics")){
    varlist[i] <- paste0("as.factor(",varlist[i],")")
  }
}
# varstring <- paste0(varlist[sample(1:length(varlist),200)], collapse = " + ", sep = "")
varstring <- paste0(varlist, collapse = " + ", sep = "")

# ## LASSO Lambda Search ######################################################

numvariables <- c()
lambdavec <- c(seq(from = 0, to = 10, by = 0.5), seq(from = 12, to = 40, by = 2))
for(lambdy in lambdavec){
  lm1 <- glmmLasso(as.formula(paste0("diagnosis ~ ",varstring)),
                   data = traindf,
                   rnd = list(Participant_ID=~1),
                   lambda=lambdy,
                   family = binomial(link = "logit"))
  summary(lm1)
  lassoFeatures <- names(lm1$coefficients[which(lm1$coefficients != 0)])
  lassoFeatures <- lassoFeatures[lassoFeatures %ni% c("(Intercept)")]
  lassoFeatures <- unique(c(lassoFeatures, "Participant_ID", "site_name", "diagnosis", "consent_age", "sex", "race", "Antibiotics"))
  numvariables <- c(numvariables, length(lassoFeatures))
}
plot(x = lambdavec, y = numvariables)

ggplot() +
  geom_point(aes(x = lambdavec, y = numvariables)) +
  geom_vline(xintercept = 5, color = "blue", linetype = "dashed") +
  theme_bw() +
  xlab("Penalty Coefficient (Lambda)") +
  ylab("Number of Included Variables")
ggsave("lambda_elbow.png", width=6, height=4, units="in", dpi=320)

stop()
## LASSO regression with lambda=5 -##########################################################################

lm1 <- glmmLasso(as.formula(paste0("diagnosis ~ ",varstring)),
                 data = traindf, 
                 rnd = list(Participant_ID=~1),
                 lambda= 5,
                 family = binomial(link = "logit"))

summary(lm1)
lassoFeatures <- names(lm1$coefficients[which(lm1$coefficients != 0)])
lassoFeatures <- lassoFeatures[grep("Participant_ID|site_name|diagnosis|consent_age|sex|race|Antibiotics|Intercept", lassoFeatures, invert = T)]
lassoFeatures <- unique(c(lassoFeatures, "Participant_ID", "site_name", "diagnosis", "consent_age", "sex", "race", "Antibiotics"))



# summary(lm1)
# lassoFeatures <- names(lm1$coefficients[which(lm1$coefficients != 0)])
# lassoFeatures <- lassoFeatures[grep("Participant_ID|site_name|diagnosis|consent_age|sex|race|Antibiotics|Intercept", lassoFeatures, invert = T)]
# lassoFeatures <- unique(c(lassoFeatures, "Participant_ID", "site_name", "diagnosis", "consent_age", "sex", "race", "Antibiotics"))
# # lassoFeatures <- unique(c(lassoFeatures, "Participant_ID", "site_name", "diagnosis", "consent_age", "sex", "race"))

## Prediction Boxplot and AUC for LASSO ######################################################

df_bestglm <- as.data.frame(traindf[,lassoFeatures])
df_bestglm$diagnosis <- as.factor(as.character(df_bestglm$diagnosis))
summary(df_bestglm$diagnosis)




varlist2 <- names(df_bestglm)[which(names(df_bestglm) %ni% c("diagnosis", "Participant_ID"))]
# varstring <- paste0(varlist[sample(1:length(varlist),200)], collapse = " + ", sep = "")
varstring2 <- paste0(varlist2, collapse = " + ", sep = "")

mymod <- lme4::glmer(as.formula(paste0("diagnosis ~ ",varstring2, " + (1|Participant_ID)")), 
                     data = df_bestglm, 
                     family = binomial)
mymodsum <- summary(mymod)
mymodsum


##--Extract feature weights--###########
# make data frame of coefficients 
mod_coef_df <- coef(mymodsum) %>% data.frame()

covar_cols <- c("site_nameCincinnati",
                "site_nameEmory",                                                                  
                "site_nameMGH",
                "site_nameMGH Pediatrics",
                "consent_age",
                "sexMale",
                "raceMore than one race",
                "raceOther",
                "raceWhite",
                "AntibioticsYes")

# remove covariate rows from coefficient dataframe
mod_coef_df_nocovar <- mod_coef_df[which(rownames(mod_coef_df) %ni% covar_cols),]


# prediction on reserved validation samples
# filter prediction dataframe to only complete cases
predictionDF <- mergeytest[complete.cases(mergeytest),]

# grab only the columns that we have coefs for (starting at index 2 removes the intercept)
predictionDF_vars <- predictionDF %>% 
  dplyr::select(rownames(mod_coef_df_nocovar)[2:nrow(mod_coef_df_nocovar)])

# make a column for the intercept
predictionDF_vars$Intercept <- 1

# make intercept the first column
predictionDF_vars <-predictionDF_vars %>% 
  dplyr::select(Intercept, everything())

##--predict the risk scores without covariates--########
pred_risk_scores <- as.matrix(predictionDF_vars) %*% as.matrix(mod_coef_df_nocovar$Estimate)

# combine the predicted and actual
pred_df <- as.data.frame(cbind(predictionDF$diagnosis, scale(pred_risk_scores))); colnames(pred_df) <- c("actual", "predicted")
pred_df$actual <- as.factor(pred_df$actual)




# null model with only covariates
mymod_ONLYcovar <- lme4::glmer(as.formula("diagnosis ~ (1|site_name) + consent_age + sex + race + Antibiotics + (1|Participant_ID)"),
                               data = df_bestglm, #[sample(1:nrow(df_bestglm),round(nrow(df_bestglm)*.6)),], 
                               family = binomial)

# library(brglm2)
# library(detectseparation)
# glm(as.formula("diagnosis ~ site_name + consent_age + sex + race + Antibiotics"), 
#                                   data = df_bestglm, 
#                        family = binomial, method="detect_separation")

# predict based on only covariates
null_model_predictions <- as.data.frame(cbind(predictionDF$diagnosis, scale(predict(mymod_ONLYcovar, predictionDF, allow.new.levels = T))))
colnames(null_model_predictions) <- c("actual", "predicted")




# make a violin plot of the prediction
boxViolinPlot <- function(auc_df = avg_par_scores, covars = "", covars_only=F){
  mycovarstring <- covars
  auc_df$actual <- as.factor(auc_df$actual)
  auc_df$predicted <- as.numeric(auc_df$predicted)
  PredPlot <- ggplot(data = auc_df, aes(x = actual, y = predicted))+
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
    ylab("Score")+
    xlab("Actual Diagnosis")
  
  
  
  
  if(covars==""){
    caseControlGLM <- glm(as.formula(auc_df$actual ~ auc_df$predicted), data = auc_df, family = "binomial", na.action = na.omit)
  } else if(covars_only==F){
    caseControlGLM <- glm(as.formula(paste0("actual ~ auc_df$predicted + ",mycovarstring)), data = auc_df, family = "binomial", na.action = na.omit)
  } else if(covars_only==T){
    caseControlGLM <- glm(as.formula(paste0("actual ~ ", mycovarstring)), data = auc_df, family = "binomial", na.action = na.omit)
  }
  
  
  
  predpr <- predict(caseControlGLM, avg_par_scores, allow.new.levels = T, type = c("response"))
  caseControlroccurve <- pROC::roc(avg_par_scores$actual ~ predpr, quiet=T, plot=T)
  caseControlroccurveCI <- pROC::roc(avg_par_scores$actual ~ predpr, ci=T, quiet=T)
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
  
  plot_df <- as.data.frame(cbind(caseControlGLM$fitted.values, auc_df$actual))
  
  colnames(plot_df) <- c("predicted", "actual")
  plot_df$actual <- as.factor(plot_df$actual)
  
  
  
  
  
  PredPlot <- PredPlot +
    geom_signif(textsize = 2.25, comparisons = list(c("0", "1")), annotations=caseControlAUCsummary, color="black", y_position = topLabLoc,tip_length=.03)+
    geom_signif(textsize = 2.25, comparisons = list(c("0", "1")), annotations=caseControlpthresh, color="black", y_position = botLabLoc,tip_length=-.03) +
    scale_y_continuous(breaks = seq(-100,100, by=2), limits = c(minPRS,maxPRS))
  
  print(caseControlAUCsummary)
  print(caseControlpthresh)
  return(PredPlot)
}

m_pred_df <- merge(pred_df,predictionDF[,c("Participant_ID","diagnosis", "site_name", "consent_age", "sex", "race", "Antibiotics")], by="row.names")
rownames(m_pred_df) <- m_pred_df$Row.names

avg_par_scores <- m_pred_df %>% 
  group_by(Participant_ID) %>%
  summarize(actual = first(diagnosis), 
            predicted = mean(predicted),
            site_name = first(site_name),
            consent_age = first(consent_age),
            sex = first(sex),
            race = first(race),
            Antibiotics = first(Antibiotics)
  )
avg_par_scores <- as.data.frame(avg_par_scores)
rownames(avg_par_scores) <- avg_par_scores$Participant_ID
avg_par_scores$Antibiotics <- as.factor(avg_par_scores$Antibiotics)


# diagnosis ~ score
print("MODEL WITH ONLY FEATURES, basic COVARIATES")
PredPlot <- boxViolinPlot(auc_df = avg_par_scores, covars = "consent_age + sex + race", covars_only=F)
PredPlot
ggsave("pred_features.png", width=2.5, height=2.5, units="in", dpi=320)

#make plot to see variation within each individual
library(ggridges)
sort_m_pred_df <- m_pred_df[order(m_pred_df$actual),]
myorder <- unique(sort_m_pred_df$Participant_ID)
sort_m_pred_df <- sort_m_pred_df %>% 
  mutate(Participant_ID = factor(Participant_ID, levels = rev(myorder)))
precolor <- sort_m_pred_df %>% 
  group_by(Participant_ID) %>%
  summarize(actual = first(diagnosis), 
            consent_age = first(consent_age),
            sex = first(sex),
            race = first(race),
            Antibiotics = first(Antibiotics)
  )
yaxiscoloring <- ifelse(precolor$actual == 1, "red", "blue")
ggplot(sort_m_pred_df, aes(x = predicted, y = Participant_ID, fill = stat(x))) +
  geom_density_ridges_gradient(scale = 2, rel_min_height = 0.05,
                               jittered_points = TRUE,
                               position = position_points_jitter(width = 0.05, height = 0),
                               point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.5) + 
  theme_minimal() + coord_cartesian(clip = "off") + # To avoid cut off
  scale_fill_viridis_c(name = NULL, option = "H", alpha = 0.5) +
  labs(title = 'Score distribution per individual') +
  theme(axis.text.y = element_text(angle = 30, hjust = 1, colour = yaxiscoloring)) +
  xlab("Score") + ylab("Participant cases (red) & controls (blue)") +
  theme(legend.position="bottom", legend.key.width = unit(1.7, 'cm'), legend.text = element_blank())
ggsave("scores_per_individual.png", width=4.31, height=5.7, units="in", dpi=320, bg='#ffffff')
stop()

# null model

null_model_predictions <- null_model_predictions %>% as.data.frame()
colnames(null_model_predictions) <- c("actual", "predicted")
null_model_predictions$actual <- as.factor(null_model_predictions$actual)

m_pred_df <- merge(null_model_predictions,predictionDF[,c("Participant_ID","diagnosis", "site_name", "consent_age", "sex", "race", "Antibiotics")], by="row.names")
rownames(m_pred_df) <- m_pred_df$Row.names

avg_par_scores <- m_pred_df %>% 
  group_by(Participant_ID) %>%
  summarize(actual = first(diagnosis), 
            predicted = mean(predicted),
            site_name = first(site_name),
            consent_age = first(consent_age),
            sex = first(sex),
            race = first(race),
            Antibiotics = first(Antibiotics)
  )
avg_par_scores <- as.data.frame(avg_par_scores)
rownames(avg_par_scores) <- avg_par_scores$Participant_ID
avg_par_scores$Antibiotics <- as.factor(avg_par_scores$Antibiotics)


print("NULL COVARIATE MODEL")
PredPlot <- boxViolinPlot(auc_df = avg_par_scores, covars = "", covars_only=F)
PredPlot
#ggsave("pred_null.png", width=2.5, height=2.5, units="in", dpi=320)

pred_df
viromics_pred_score_featuresonly <- tibble::rownames_to_column(pred_df, "External_ID")
write.table(viromics_pred_score_featuresonly, "viromics_features_scores.txt", sep="\t", col.names=T, row.names=F, quote=F)

# remove the intercept
featureplot_df <- mod_coef_df_nocovar[2:nrow(mod_coef_df_nocovar),]
featureplot_df$Feature <- rownames(featureplot_df)
featureplot_df <- featureplot_df %>% dplyr::arrange(Estimate)


featureplot_df %>% 
  ggplot(mapping = aes(y=Feature, x=Estimate)) +
  geom_point()

# make our plotting dataframe
df2 <- featureplot_df %>%
  tibble::rownames_to_column() %>%
  dplyr::rename("variable" = rowname) %>%
  dplyr::arrange(Estimate) %>%
  dplyr::mutate(variable = forcats::fct_inorder(variable))

plot_varimp2 <- ggplot2::ggplot(df2) +
  geom_segment(
    aes(
      x = variable,
      y = 0,
      xend = variable,
      yend = Estimate
    ),
    size = 1.5,
    alpha = 0.7
  ) +
  geom_point(aes(x = variable, y = Estimate, col = variable),
             size = 4,
             show.legend = F) +
  coord_flip() +
  labs(y = "Weight", x = NULL, title = "") +
  theme_bw() +
  theme(legend.title = element_text(size = 14)) +
  theme(
    axis.text.x = element_text(
      color = "black",
      size = 13,
      angle = 0,
      hjust = .5,
      vjust = .5
    ),
    axis.text.y = element_text(
      color = "black",
      size = length(unique(df2$variable))*1,
      angle = 0
    ),
    axis.title.x = element_text(
      color = "black",
      size = 13,
      angle = 0
    ),
    axis.title.y = element_text(
      color = "black",
      size = 13,
      angle = 90
    )
  )
plot_varimp2
# ## Prediction Boxplot and AUC for LASSO ######################################################
# df_bestglm <- as.data.frame(traindf[,c(lassoFeatures)])
# df_bestglm$diagnosis <- as.factor(as.character(df_bestglm$diagnosis))
# summary(df_bestglm$diagnosis)
# 
# varlist2 <- names(df_bestglm)[which(names(df_bestglm) %ni% c("diagnosis", "Participant_ID"))]
# # varstring <- paste0(varlist[sample(1:length(varlist),200)], collapse = " + ", sep = "")
# varstring2 <- paste0(varlist2, collapse = " + ", sep = "")
# 
# mymod <- lme4::glmer(as.formula(paste0("diagnosis ~ ",varstring2, " + (1|Participant_ID)")),
#                      data = df_bestglm,
#                      family = binomial)
# mymodsum <- summary(mymod)
# # prediction on reserved validation samples
# # filter prediction dataframe to only complete cases
# predictionDF <- mergeytest[complete.cases(mergeytest),]
# pred_df <- as.data.frame(cbind(predictionDF$diagnosis, scale(predict(mymod, predictionDF, allow.new.levels = T)))); colnames(pred_df) <- c("actual", "predicted")
# pred_df$actual <- as.factor(pred_df$actual)
# # make a violin plot of the prediction
# boxViolinPlot <- function(pred_df = pred_df, predictionDF = predictionDF){
#   PredPlot <- ggplot(data = pred_df, aes(x = actual, y = predicted))+
#     scale_fill_viridis_d( option = "D")+
#     theme_dark()+
#     geom_violin(fill = "gray70",alpha=0.4, position = position_dodge(width = .5),size=1,color="gray22",width=.5,lwd=.2) +
#     geom_boxplot(fill = "gray95",notch = F, shape=21, outlier.size = -1, color="gray32",lwd=.5, alpha = .75)+
#     theme(plot.title = element_text(hjust = 0.5))+
#     # theme(axis.title.x = element_text(size=14))+
#     theme(axis.text.x = element_text(colour = "black"))+
#     theme(axis.text.y = element_text(colour = "black"))+
#     theme( axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"))+
#     theme(panel.grid.major.x = element_blank())+
#     theme(panel.background = element_rect(fill = 'white'), panel.grid = element_line(color='gray80'))+
#     # labs(title = addToTitle)+
#     ylab("Score")+
#     xlab("Actual Diagnosis")
# 
#   pred_df$actual <- as.factor(pred_df$actual)
#   pred_df$predicted <- as.numeric(pred_df$predicted)
#   caseControlGLM <- glm(as.formula(pred_df$actual ~ pred_df$predicted), data = pred_df, family = "binomial", na.action = na.omit)
#   predpr <- predict(caseControlGLM, predictionDF, allow.new.levels = T, type = c("response"))
#   caseControlroccurve <- pROC::roc(predictionDF$diagnosis ~ predpr, quiet=T, plot=F)
#   caseControlroccurveCI <- pROC::roc(predictionDF$diagnosis ~ predpr, ci=T, quiet=T)
#   # caseControlplot <- plot(caseControlroccurve, main=paste("Case vs Control AUC =", round(caseControlroccurve$auc, 3)))
# 
#   caseControlp <- formatC(coef(summary(caseControlGLM))[,4][2], format = "e", digits = 0)
#   caseControlOR <- exp(cbind("Odds ratio" = coef(caseControlGLM), confint.default(caseControlGLM, level = 0.95)))
#   if (as.numeric(strsplit(caseControlp,"")[[1]][1]) == 1){
#     caseControlpthresh <- paste0("OR [95% CI] = ", format(round(as.numeric(caseControlOR[2,1]), 2), nsmall = 2), " [", format(round(as.numeric(caseControlOR[2,2]), 2),nsmal = 2), ", ", format(round(as.numeric(caseControlOR[2,3]), 2),nsmal = 2),"], ","p < ", as.numeric(caseControlp)/as.numeric(strsplit(caseControlp,"")[[1]][1]))
#   }else{
#     caseControlpthresh <- paste0("OR [95% CI] = ", format(round(as.numeric(caseControlOR[2,1]), 2), nsmall = 2), " [", format(round(as.numeric(caseControlOR[2,2]), 2),nsmal = 2), ", ", format(round(as.numeric(caseControlOR[2,3]), 2),nsmal = 2),"], ","p < ", 10*as.numeric(caseControlp)/as.numeric(strsplit(caseControlp,"")[[1]][1]))
#   }
#   if (as.numeric(caseControlp) >= 0.001){
#     caseControlpthresh <- paste0("OR [95% CI] = ", format(round(as.numeric(caseControlOR[2,1]), 2), nsmall = 2), " [", format(round(as.numeric(caseControlOR[2,2]), 2),nsmal = 2), ", ", format(round(as.numeric(caseControlOR[2,3]), 2),nsmal = 2),"], ","p = ", formatC(as.numeric(caseControlp), format = "g"))
#   }
#   caseControlAUCsummary <- paste0("AUC [95% CI] = ", format(round(as.numeric(caseControlroccurve$auc), 2), nsmall = 2), " [", format(round(as.numeric(caseControlroccurveCI$ci[1]), 2),nsmal = 2), ", ", format(round(as.numeric(caseControlroccurveCI$ci[3]), 2),nsmal = 2),"]")
#   caseControlORsummary <- caseControlpthresh
# 
#   label_disty = .13
#   minPRS <- min(pred_df$predicted)-(abs(max(pred_df$predicted)-min(pred_df$predicted)))*(label_disty*1.5*1.02)
#   maxPRS <- max(pred_df$predicted) + (abs(max(pred_df$predicted)-min(pred_df$predicted)))*(label_disty*1.5*1.02)
#   botLabLoc <- min(pred_df$predicted)-(abs(max(pred_df$predicted)-min(pred_df$predicted)))*(label_disty*1.5)
#   topLabLoc <- max(pred_df$predicted) + (abs(max(pred_df$predicted)-min(pred_df$predicted)))*(label_disty*0.21)
# 
#   PredPlot <- PredPlot +
#     geom_signif(textsize = 2.25, comparisons = list(c("0", "1")), annotations=caseControlAUCsummary, color="black", y_position = topLabLoc,tip_length=.03)+
#     geom_signif(textsize = 2.25, comparisons = list(c("0", "1")), annotations=caseControlpthresh, color="black", y_position = botLabLoc,tip_length=-.03) +
#     scale_y_continuous(breaks = seq(-100,100, by=2), limits = c(minPRS,maxPRS))
#   return(PredPlot)
# }
# 
# PredPlotLasso <- boxViolinPlot(pred_df = pred_df, predictionDF = predictionDF)
# PredPlotLasso
# 
# 
# 
# # ## Mixed Effects Random Forests via MixRF ######################################################
# 
# myX <- as.data.frame(df_bestglm[,names(df_bestglm) %ni% c("diagnosis", "Participant_ID")])
# 
# myMixRF <- MixRF::MixRF(
#   Y = as.numeric(pull(df_bestglm, diagnosis)),
#   X = myX,
#   random = "(1|Participant_ID)",
#   MaxIterations = 10,
#   data = df_bestglm
# )
# 
# # ## RF Variable Importance ######################################################
# caret::varImp(myMixRF$forest)
# df_tmp <- caret::varImp(myMixRF$forest)
# newdf <- df_tmp; colnames(newdf) <- c("imp")
# # newdf <- data.frame(cbind(df_tmp$Overall)); colnames(newdf) <- c("imp"); rownames(newdf) <- rownames(df_tmp)
# newdf$imp <- as.numeric(newdf$imp)
# df2 <- newdf %>%
#   tibble::rownames_to_column() %>%
#   dplyr::rename("variable" = rowname) %>%
#   dplyr::arrange(imp) %>%
#   dplyr::mutate(variable = forcats::fct_inorder(variable))
# plot_varimp2 <- ggplot2::ggplot(df2) +
#   geom_segment(
#     aes(
#       x = variable,
#       y = 0,
#       xend = variable,
#       yend = imp
#     ),
#     size = 1.5,
#     alpha = 0.7
#   ) +
#   geom_point(aes(x = variable, y = imp, col = variable),
#              size = 4,
#              show.legend = F) +
#   coord_flip() +
#   labs(y = "Importance", x = NULL, title = "") +
#   theme_bw() +
#   theme(legend.title = element_text(size = 14)) +
#   theme(
#     axis.text.x = element_text(
#       color = "black",
#       size = 13,
#       angle = 0,
#       hjust = .5,
#       vjust = .5
#     ),
#     axis.text.y = element_text(
#       color = "black",
#       size = length(unique(df2$variable))*1/3,
#       angle = 0
#     ),
#     axis.title.x = element_text(
#       color = "black",
#       size = 13,
#       angle = 0
#     ),
#     axis.title.y = element_text(
#       color = "black",
#       size = 13,
#       angle = 90
#     )
#   )
# plot_varimp2
# 
# 
# levels(mergeytest[["sex"]])
# levels(myX[["sex"]])
# levels(mergeytest[["site_name"]])
# levels(myX[["site_name"]])
# levels(mergeytest[["race"]])
# levels(myX[["race"]])
# 
# # ## RF Prediction ######################################################
# # filter prediction dataframe to only complete cases
# predictionDF <- mergeytest[complete.cases(mergeytest),]
# predictionDF$consent_age <- as.factor(as.character(predictionDF$consent_age))
# common <- names(myX)[names(myX) %in% names(predictionDF)]
# predictionDF[common] <- lapply(common, function(x) {
#   match.fun(paste0("as.", class(myX[[x]])))(predictionDF[[x]])
# })
# str(predictionDF[,common])
# str(myX[,common])
# 
# # make the factors match between testing and training
# myX$diagnosis <- as.numeric(pull(df_bestglm, diagnosis))
# predictionDF <- predictionDF[,c(common,"diagnosis")]
# predictionDF <- rbind(myX[1, ] , predictionDF)
# predictionDF <- predictionDF[-1,]
# 
# # predictionDF$diagnosis <- as.factor(predictionDF$diagnosis)
# postPred_RF <- as.data.frame(scale(predict(myMixRF$forest, newdata = predictionDF)))
# # postPred_RF <- as.data.frame(predict(myglm, newdata = mergeytest)); postPred_RF <- ifelse(postPred_RF[,1] > 0.25, "1", "0")
# 
# comparePrediction <- cbind(as.character(predictionDF$diagnosis), postPred_RF)
# comparePrediction <- as.data.frame(comparePrediction)
# colnames(comparePrediction) <- c("actual", "predicted")
# (comparePrediction$actual == comparePrediction$prediction)
# summary((comparePrediction$actual == comparePrediction$prediction))
# # View(comparePrediction)
# 
# PredPlotRF <- boxViolinPlot(pred_df = comparePrediction, predictionDF = predictionDF)
# PredPlotRF


## Output predictions from lasso and MixRF --######
# Rerun at your own risk
# viromics_glmer_scores <- pred_df
# viromics_glmer_scores <- tibble::rownames_to_column(viromics_glmer_scores, "External_ID")
# write.table(viromics_glmer_scores, "./viromics/viromics_glmer_scores.txt", sep="\t", col.names=T, row.names=F, quote=F)
# 
# 
# viromics_MixRF_scores <- pred_df
# viromics_MixRF_scores <- tibble::rownames_to_column(viromics_MixRF_scores, "External_ID")
# write.table(viromics_MixRF_scores, "./viromics/viromics_MixRF_scores.txt", sep="\t", col.names=T, row.names=F, quote=F)



# ## run logistic regression for each feature --###################################
# featureNames <- c()
# betas <- c()
# pvals <- c()
# for(i in 1:(ncol(mergey)-1)){
#   # for(i in 1:1){
#   # randName <- names(mergey)[sample(1:length(names(mergey)),1)]
#   randName <- names(mergey)[i]
#   mergeysub <- mergey[,c(randName, "Participant ID", "race", 
#                          "consent_age", "site_name", 
#                          "sex", "Antibiotics", "diagnosis")]
#   
#   colnames(mergeysub) <- c(randName, "Participant ID", "race", 
#                            "consent_age", "site_name", 
#                            "sex", "Antibiotics", "diagnosis")
#   mergeysub$`Participant ID` <- as.factor(mergeysub$`Participant ID`)
#   
#   mymod <- glmer(as.formula(
#     paste0("diagnosis ~ `", randName,"` + consent_age + sex + race + Antibiotics + (1|`Participant ID`) + (1|`site_name`)")), 
#     data = mergeysub, family = "binomial", 
#     #control = glmerControl(optimizer ="Nelder_Mead")
#     control = glmerControl(optimizer ="bobyqa")
#   )
#   #mymod <- glm(as.formula(paste0("diagnosis ~ `", randName,"` + consent_age + sex + race + Antibiotics + (1|`Participant ID`) + (1|`site_name`)")), data = mergeysub, family = "binomial")
#   #mymod <- glm(as.formula(paste0("diagnosis ~ `",randName,"`")), data = mergeysub, family = "binomial")
#   mymodsum <- summary(mymod)
#   featureNames <- c(featureNames, randName)
#   betas <- c(betas, mymodsum$coefficients[2,1])
#   pvals <- c(pvals, mymodsum$coefficients[2,4])
#   # validate on testing set
#   # if(mymodsum$coefficients[2,4] < 0.05){
#   #   print(randName)
#   #   print(cor(mergeytest$diagnosis, predict(mymod, mergeytest)))
#   # }
#   # Sys.sleep(1)
#   # summary(mymod)
#   # print(randName)
#   # print(mymodsum$aic)
# }
# 
# 
# 
# 
# # multiple test p value for significance (bonferoni)
# bonfsigthresh <- 0.05/(ncol(mergey)-1)
# # make dataframe
# tenResults <- as.data.frame(cbind(featureNames,betas,pvals))
# tenResults$featureNames <- as.character(tenResults$featureNames)
# tenResults$betas <- as.numeric(as.character(tenResults$betas))
# tenResults$pvals <- as.numeric(as.character(tenResults$pvals))
# 
# 
# # column for plot colors
# tenResults$mycolors <- NA
# tenResults$mycolors[tenResults$pvals < bonfsigthresh] <- "sig"
# tenResults$mycolors[tenResults$pvals >= bonfsigthresh] <- "notsig"
# tenResults$mycolors <- as.factor(tenResults$mycolors)
# 
# 
# # column for plot labels
# tenResults$delabels <- NA
# tenResults$delabels[tenResults$pvals < bonfsigthresh] <- tenResults$featureNames[tenResults$pvals < bonfsigthresh]
# tenResults$delabels[tenResults$pvals >= bonfsigthresh] <- NA
# tenResults$delabels <- as.character(tenResults$delabels)
# 
# # make volcano plot
# bplot <- ggplot(aes(x = betas, y = -log10(pvals), col=mycolors, label=delabels), data = tenResults) +
#   # geom_bar(stat="identity", fill = "steelblue") + theme_minimal()
#   geom_point() +
#   theme_minimal() +
#   geom_text()
# bplot
# 
# 
# 
# # NO BONFERONNI
# # column for plot colors
# tenResults$significant <- NA
# tenResults$significant[tenResults$pvals < 0.05] <- "sig"
# tenResults$significant[tenResults$pvals >= 0.05] <- "notsig"
# tenResults$significant <- as.factor(tenResults$significant)
# 
# 
# # column for plot labels
# tenResults$labels <- NA
# tenResults$labels[tenResults$pvals < 0.05] <- tenResults$featureNames[tenResults$pvals < 0.05]
# tenResults$labels[tenResults$pvals >= 0.05] <- NA
# tenResults$labels <- as.character(tenResults$labels)
# 
# # make volcano plot
# bplot <- ggplot(aes(x = betas, y = -log10(pvals), col=significant, label=labels), data = tenResults) +
#   # geom_bar(stat="identity", fill = "steelblue") + theme_minimal()
#   geom_point() +
#   theme_minimal() +
#   geom_text()
# bplot
# 
# 
# 
# # extract the features that were significant and run a glm on the full model
# sigFeatures <- tenResults$featureNames[tenResults$mycolors == "sig"]
# df_bestglm <- mergey[,c(sigFeatures,"diagnosis")]
# mymod <- glm(as.formula(paste0("diagnosis ~ .")), data = df_bestglm, family = "binomial")
# mymodsum <- summary(mymod)
# # prediction on reserved validation samples
# pred_df <- as.data.frame(cbind(mergeytest$diagnosis, predict(mymod, mergeytest))); colnames(pred_df) <- c("actual", "predicted")
# pred_df$actual <- as.factor(pred_df$actual)
# # make a violin plot of the prediction
# ggplot(data = pred_df, aes(x = actual, y = predicted))+
#   scale_fill_viridis_d( option = "D")+
#   theme_dark()+
#   geom_violin(fill = "gray70",alpha=0.4, position = position_dodge(width = .5),size=1,color="gray22",width=.5,lwd=.2) +
#   geom_boxplot(fill = "gray95",notch = F, shape=21, outlier.size = -1, color="gray32",lwd=.5, alpha = .75)+
#   theme(plot.title = element_text(hjust = 0.5))+
#   # theme(axis.title.x = element_text(size=14))+
#   theme(axis.text.x = element_text(colour = "black"))+
#   theme(axis.text.y = element_text(colour = "black"))+
#   theme( axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"))+
#   theme(panel.grid.major.x = element_blank())+
#   theme(panel.background = element_rect(fill = 'white'), panel.grid = element_line(color='gray80'))+
#   # labs(title = addToTitle)+
#   ylab("Predicted Diagnosis")+
#   xlab("Actual Diagnosis")
# 
# 
# # bestglm
# # function begins with a data frame containing explanatory variables and response variables. 
# # The response variable should be in the last column. 
# # Varieties of goodness-of-fit criteria can be specified in the IC argument.
# # bestglm can only handle 15 variables.
# subfeats <- sample(names(df_bestglm), 14)
# subfeats <- subfeats[which(subfeats != "diagnosis")]; subfeats <- c(subfeats, "diagnosis")
# bglm <- bestglm::bestglm(df_bestglm[,subfeats], family = binomial, IC = "BIC")
# # prediction on reserved validation samples
# pred_df <- as.data.frame(cbind(mergeytest$diagnosis, predict(bglm$BestModel, mergeytest))); colnames(pred_df) <- c("actual", "predicted")
# pred_df$actual <- as.factor(pred_df$actual)
# # make a violin plot of the prediction
# ggplot(data = pred_df, aes(x = actual, y = predicted))+
#   scale_fill_viridis_d( option = "D")+
#   theme_dark()+
#   geom_violin(fill = "gray70",alpha=0.4, position = position_dodge(width = .5),size=1,color="gray22",width=.5,lwd=.2) +
#   geom_boxplot(fill = "gray95",notch = F, shape=21, outlier.size = -1, color="gray32",lwd=.5, alpha = .75)+
#   theme(plot.title = element_text(hjust = 0.5))+
#   # theme(axis.title.x = element_text(size=14))+
#   theme(axis.text.x = element_text(colour = "black"))+
#   theme(axis.text.y = element_text(colour = "black"))+
#   theme( axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"))+
#   theme(panel.grid.major.x = element_blank())+
#   theme(panel.background = element_rect(fill = 'white'), panel.grid = element_line(color='gray80'))+
#   # labs(title = addToTitle)+
#   ylab("Predicted Diagnosis")+
#   xlab("Actual Diagnosis")
