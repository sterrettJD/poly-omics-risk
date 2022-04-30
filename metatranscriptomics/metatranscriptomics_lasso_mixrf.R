 
setwd("/Users/johnsterrett/Research-Projects/Team-rotation/poly-omics-scores/metatranscriptomics/")

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
# install.packages("bestglm)
# install.packages("tree")
library(data.table)
library(ggplot2)
library(R.utils)
library(tidyverse)
library(UpSetR)
library(cowplot)
library(compositions)
library(bestglm)
library(MASS)
library(tree)
library(lme4)
library(pROC)
library(glmmLasso)
library(ggpubr)
library(MixRF)
library(caret)
library(fmsb)
'%ni%' <- Negate('%in%')

## metadata --###################################
metadata <- fread("https://ibdmdb.org/tunnel/products/HMP2/Metadata/hmp2_metadata.csv", header=T, stringsAsFactors=T)
str(metadata)
# metadata <- subset(metadata, data_type == "metagenomics")

# make barplt by data type
par(mar = c(9, 4, 2, 2) + 1)
barplot(summary(metadata$data_type),las=2)

# make upset plot for sample overlap by datatype
listInput <- list(
  # biopsy_16S = c(as.character(subset(metadata, data_type == "biopsy_16S")[["External ID"]])),
  # host_genome = c(as.character(subset(metadata, data_type == "host_genome")[["External ID"]])),
  # host_transcriptomics  = c(as.character(subset(metadata, data_type == "host_transcriptomics")[["External ID"]])),
  # methylome = c(as.character(subset(metadata, data_type == "methylome")[["External ID"]])),
  # proteomics = c(as.character(subset(metadata, data_type == "proteomics")[["External ID"]])),
  # serology  = c(as.character(subset(metadata, data_type == "serology")[["External ID"]])),
  # stool_16S = c(as.character(subset(metadata, data_type == "stool_16S")[["External ID"]])),
  metabolomics  = c(as.character(subset(metadata, data_type == "metabolomics")[["External ID"]])), 
  metagenomics = c(as.character(subset(metadata, data_type == "metagenomics")[["External ID"]])),  
  metatranscriptomics = c(as.character(subset(metadata, data_type == "metatranscriptomics")[["External ID"]])), 
  viromics  = c(as.character(subset(metadata, data_type == "viromics")[["External ID"]]))
)
upset(fromList(listInput), order.by = "freq")

# make a list of the ids that are in the intersection of metabolomics, metagenomics, viromics, metatranscriptomics
i1 <- intersect(c(as.character(subset(metadata, data_type == "metabolomics")[["External ID"]])), 
                c(as.character(subset(metadata, data_type == "metagenomics")[["External ID"]]))
)

i2 <- intersect(i1, c(as.character(subset(metadata, data_type == "viromics")[["External ID"]])))
i3 <- intersect(i2, c(as.character(subset(metadata, data_type == "metatranscriptomics")[["External ID"]])))
length(i3)
idlist <- unique(i3)
# idlist <- as.data.frame(idlist); colnames(idlist) <- c("idlist")
# head(idlist)
# idlistsep <- idlist %>% separate(idlist,into=c("idlist","junk"),convert=TRUE,sep="_profi")
# idlistsep <- idlistsep$idlist

# subset to jsut one datatype to explore diagnosis counts
metadata <- subset(metadata, data_type == "metatranscriptomics")
# metadata <- subset(metadata, `External ID` %in% idlist)

# fill missinig enries with NA and view the columns that are mostly non-NA
metadata[metadata==""] <- NA
isna <- sapply(metadata, function(x) sum(is.na(x)))
isna[isna < 100]
summary(metadata$diagnosis)
metadata1 <- as.data.frame(metadata[,c("External ID","diagnosis")])
metadata1$diagnosis <- as.character(metadata1$diagnosis)
metadata1$diagnosis[metadata1$diagnosis == "UC"] <- 1
metadata1$diagnosis[metadata1$diagnosis == "CD"] <- 1
metadata1$diagnosis[metadata1$diagnosis == "nonIBD"] <- 0
metadata1$diagnosis <- as.numeric(metadata1$diagnosis)

# training_metadata <- fread("https://raw.githubusercontent.com/sterrettJD/poly-omics-risk/main/training_metadata.txt?token=GHSAT0AAAAAABQ2LF4FWL6XDQDJLZ4C4F5SYS5WH4Q", sep = "\t")
training_metadata <- fread("../training_metadata.txt", sep = "\t")
training_metadata <- subset(training_metadata, data_type == "metatranscriptomics")
training_metadata[training_metadata==""] <- NA
isna <- sapply(training_metadata, function(x) sum(is.na(x)))
isna[isna < 100]
summary(training_metadata$diagnosis)
training_metadata <- as.data.frame(training_metadata[,c("External ID", "Participant ID", "site_name", "diagnosis", "consent_age", "sex", "race", "Antibiotics")])
training_metadata$diagnosis <- as.character(training_metadata$diagnosis)
training_metadata$diagnosis[training_metadata$diagnosis == "UC"] <- 1
training_metadata$diagnosis[training_metadata$diagnosis == "CD"] <- 1
training_metadata$diagnosis[training_metadata$diagnosis == "nonIBD"] <- 0
training_metadata$diagnosis <- as.numeric(training_metadata$diagnosis)
training_metadata$site_name <- as.factor(training_metadata$site_name)
training_metadata$sex <- as.factor(training_metadata$sex)
training_metadata$race <- as.factor(training_metadata$race)

# testing_metadata <- fread("https://raw.githubusercontent.com/sterrettJD/poly-omics-risk/main/testing_metadata.txt?token=GHSAT0AAAAAABQ2LF4EWNW2TVN2PVOHJO6MYS5WHLQ", sep = "\t")
testing_metadata <- fread("../testing_metadata.txt", sep = "\t")
testing_metadata <- subset(testing_metadata, data_type == "metatranscriptomics")
testing_metadata[testing_metadata==""] <- NA
isna <- sapply(testing_metadata, function(x) sum(is.na(x)))
isna[isna < 100]
summary(testing_metadata$diagnosis)
testing_metadata <- as.data.frame(testing_metadata[,c("External ID", "Participant ID", "site_name", "diagnosis", "consent_age", "sex", "race", "Antibiotics")])
testing_metadata$diagnosis <- as.character(testing_metadata$diagnosis)
testing_metadata$diagnosis[testing_metadata$diagnosis == "UC"] <- 1
testing_metadata$diagnosis[testing_metadata$diagnosis == "CD"] <- 1
testing_metadata$diagnosis[testing_metadata$diagnosis == "nonIBD"] <- 0
testing_metadata$diagnosis <- as.numeric(testing_metadata$diagnosis)
testing_metadata$site_name <- as.factor(testing_metadata$site_name)
testing_metadata$sex <- as.factor(testing_metadata$sex)
testing_metadata$race <- as.factor(testing_metadata$race)

# Reclassify the one individual in testing that has a unique race class:
summary(testing_metadata$race)
summary(training_metadata$race)
testing_metadata$race <- as.character(testing_metadata$race)
# clump "American Indian or Alaska Native" into "Other"
testing_metadata$race[which(testing_metadata$race == "American Indian or Alaska Native")] <- "Other"
testing_metadata$race <- as.factor(testing_metadata$race)
summary(testing_metadata$race)
summary(training_metadata$race)

## Metagenomes taxonomic_profiles_3 --###################################
# fname <- fread("https://ibdmdb.org/tunnel/products/HMP2/WGS/1818/taxonomic_profiles_3.tsv.gz", header=T)
fname <- fread("https://ibdmdb.org/tunnel/products/HMP2/MTX/1750/pathabundances_3.tsv.gz", header=T)
fname <- as.data.frame(fname)

df_3 <- fname
print("data dimensions:")
print(dim(df_3))
print("data glimpse:")
str(df_3[1:5,1:4])
rownames(df_3) <- df_3$`Feature\\Sample`
df_3$`Feature\\Sample` <- NULL
# View(as.data.frame(rownames(df_3)))
# separate dataframe out by rows which were annotated vs rows which weren't
commentRow <- grep("# ", rownames(df_3), invert = F)

mygrep <- grep("UNINTEGRATED|UNGROUPED|UNMAPPED|unclassified", rownames(df_3), invert = F) # get any unmapped data
mygrepni <- which(1:length(rownames(df_3)) %ni% mygrep) # grep the lines that aren't unmapped

if(length(commentRow) > 0){
  mygrepni <- mygrepni[which(mygrepni != commentRow)]
}

grouped_df_3 <- df_3[mygrepni,]

# Filter out ones that contain taxonomy
mygrep <- grep("\\|g_", rownames(grouped_df_3), invert = F) 
mygrepni <- which(1:length(rownames(grouped_df_3)) %ni% mygrep)

grouped_df_3 <- grouped_df_3[mygrepni,]

view(grouped_df_3)

## Preprocessing, make data compositional for species --###################################
# grouped_df_3 <- df_3[mygrepni,]
# rownames(grouped_df_3) <- rownames(df_3)[mygrepni]

num_grouped_df_3 <- mutate_all(grouped_df_3, function(x) as.numeric(as.character(x)))

# is it compositional?
summary(colSums(num_grouped_df_3))
hist(colSums(num_grouped_df_3))
length(which(colSums(num_grouped_df_3) < 5))
colnames(num_grouped_df_3)[which(colSums(num_grouped_df_3) < 5)]
summary(rowSums(num_grouped_df_3))
hist(rowSums(num_grouped_df_3))

# remove the 11 individuals that have really low colSums....
num_grouped_df_3 <- num_grouped_df_3[,which(colSums(num_grouped_df_3) > 5)]
hist(colSums(num_grouped_df_3))
summary(colSums(num_grouped_df_3))

# filter to rows where rowmean == 100
# View(num_grouped_df_3[rowMeans(num_grouped_df_3) > 10,])

# normalize by sample sum
# if(normalize == T){
#     num_grouped_df_3 <- as.data.frame(scale(num_grouped_df_3, center=FALSE, scale=colSums(num_grouped_df_3)))
# }

# add epsilon to all entries to set up for center log transform
pepsi <- 1E-06
num_grouped_df_3 <- num_grouped_df_3 + pepsi

## CLT transform --###################################
# transform by center log transform across rows (features)
clr_num_grouped_df_3 <- compositions::clr(num_grouped_df_3,)
# clr_num_grouped_df_3 <- as.data.frame(apply(num_grouped_df_3, 1, clr))

alldf_pretransform <- as.data.frame(as.numeric(array(as.matrix(num_grouped_df_3)))); colnames(alldf_pretransform) <- c("alldf_pretransform")
hist(alldf_pretransform$alldf_pretransform)
alldf_posttransform <- as.data.frame(as.numeric(array(as.matrix(clr_num_grouped_df_3)))); colnames(alldf_posttransform) <- c("alldf_posttransform")
hist(alldf_posttransform$alldf_posttransform)

# rename the columns for the merge with metadata
namesies <- as.data.frame(colnames(clr_num_grouped_df_3)); colnames(namesies) <- c("namesies")
head(namesies)
namesiessep <- namesies %>% separate(namesies,into=c("namesies","junk"),convert=TRUE,sep="_patha")
colnames(clr_num_grouped_df_3) <- namesiessep$namesies

intersecty <- intersect(c(as.character(namesiessep$namesies)), 
                c(as.character(metadata1$`External ID`))
)

dim(namesiessep)
dim(metadata1)
length(intersecty)

# transpose transformed data
transp_clr_num_grouped_df_3 <- as.data.frame(t(clr_num_grouped_df_3))
transp_clr_num_grouped_df_3$`External ID` <- rownames(transp_clr_num_grouped_df_3)

## Merge with metadata --###################################
# merge with diagnosis
mergey <- merge(transp_clr_num_grouped_df_3, training_metadata, by = "External ID")
dim(transp_clr_num_grouped_df_3)
dim(mergey)
length(unique(mergey[,c("Participant ID")]))
# make validation dataset that isn't used for training (the ids that are in the polyomic list)
mergeytest <- merge(transp_clr_num_grouped_df_3, testing_metadata, by = "External ID")
dim(transp_clr_num_grouped_df_3)
dim(mergeytest)
length(unique(mergeytest[,c("Participant ID")]))

# mergeytest <- mergey[which(as.character(mergey$`External ID`) %in% as.character(idlist)),]
# dim(mergeytest)
# make training dataset (the ids that are NOT in the polyomic list)
# mergey <- subset(mergey, `External ID` %ni% idlist)
# dim(mergey)
# make the ids the rownames for each dataframe, and then remove that column
rownames(mergeytest) <- mergeytest$`External ID`
mergeytest$`External ID` <- NULL
rownames(mergey) <- mergey$`External ID`
mergey$`External ID` <- NULL

# remove any columns that have no/little variation between samples...
mergeytest_colsd <- apply(mergeytest[complete.cases(mergeytest),], 2, sd, na.rm=T)
mergey_colsd <- apply(mergey[complete.cases(mergey),], 2, sd, na.rm=T)
# qthresh <- quantile(colsd, 0.05, na.rm=T)
hist(as.numeric(mergeytest_colsd))
hist(as.numeric(mergey_colsd))
toKeep <- c(intersect(c(names(mergeytest_colsd)[which(mergeytest_colsd > 1)]), c(names(mergey_colsd)[which(mergey_colsd > 1)])),
            c("External ID", "Participant ID", "site_name", "diagnosis", "consent_age", "sex", "race", "Antibiotics")
            )
length(toKeep)
length(mergey_colsd)
mergeytest <- mergeytest[,which(names(mergeytest) %in% toKeep)]
mergey <- mergey[,which(names(mergey) %in% toKeep)]

# rename the column names with just the species
cn <- colnames(mergey)
cn
# 
# cn <- cn %>% separate(cn,into=c("junk","species"),convert=TRUE,sep="\\|s__")
# cn <- cn$species
# cn[is.na(cn)] <- c("Participant_ID", "site_name", "diagnosis", "consent_age", "sex", "race", "Antibiotics")
# head(cn,20)
# tail(cn,20)
cn <- gsub(pattern = "\\[", replacement = "", cn)
cn <- gsub(pattern = "\\]", replacement = "", cn)
cn <- gsub(pattern = "-", replacement = "_", cn)
cn <- gsub(pattern = " ", replacement = "_", cn)
cn <- gsub(pattern = ":", replacement = "__", cn)
cn <- gsub(pattern = ",", replacement = "_", cn)
cn <- gsub(pattern = ";", replacement = "_", cn)
cn <- gsub(pattern = "&", replacement = "_and_", cn)
cn <- gsub(pattern = "\\(", replacement = "", cn)
cn <- gsub(pattern = "\\)", replacement = "", cn)
cn <- gsub(pattern = "'", replacement = "prime", cn)
cn <- gsub(pattern = "\\+", replacement = "plus", cn)
cn <- gsub(pattern = "\\/", replacement = "_", cn)


cn[which(sapply(strsplit(cn,""),"[[",1) %in% 0:9)] <- paste0("P_", cn[which(sapply(strsplit(cn,""),"[[",1) %in% 0:9)])


cn

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


str(mergey[,1:30])
str(mergey[,(ncol(mergey)-30):ncol(mergey)])

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

## LASSO Lambda Search ######################################################


# remove colinear columns

pwy_df <- traindf%>% dplyr::select( -Participant_ID,
                 -site_name,
                 -diagnosis,
                 -consent_age,
                 -sex,
                 -race,
                 -Antibiotics)

tmp <- cor(pwy_df)
tmp[upper.tri(tmp)] <- 0
diag(tmp) <- 0

# Above two commands can be replaced with 
# tmp[!lower.tri(tmp)] <- 0

data.new <- 
  pwy_df[, !apply(tmp, 2, function(x) any(abs(x) > 0.95, na.rm = TRUE))]
head(data.new)

dim(data.new)
dim(pwy_df)

train_df <- cbind(data.new,
      traindf%>% dplyr::select(Participant_ID,
                               site_name,
                               diagnosis,
                               consent_age,
                               sex,
                               race,
                               Antibiotics)) %>%
  as.data.frame()

dim(train_df)

varlist_nocolin <- c(varlist[which(varlist %in% colnames(train_df))], 
                    "site_name", 
                    "sex", 
                    "race", 
                    "Antibiotics")
varstring_nocolin <- paste0(varlist_nocolin, collapse = " + ", sep = "")


numvariables <- c()
lambdavec <- seq(from = 10, to = 110, by = 5)
for(lambdy in lambdavec){
  lm1 <- glmmLasso(as.formula(paste0("diagnosis ~ ",varstring_nocolin)),
                   data = train_df,
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

lassoFeatures

lm1 <- glmmLasso(as.formula(paste0("diagnosis ~ ",varstring_nocolin)),
                 data = traindf, 
                 rnd = list(Participant_ID=~1),
                 lambda=25,
                 family = binomial(link = "logit"))
summary(lm1)
lassoFeatures <- names(lm1$coefficients[which(lm1$coefficients != 0)])
lassoFeatures <- lassoFeatures[grep("Participant_ID|site_name|diagnosis|consent_age|sex|race|Antibiotics|Intercept", lassoFeatures, invert = T)]
lassoFeatures <- unique(c(lassoFeatures, "Participant_ID", "site_name", "diagnosis", "consent_age", "sex", "race", "Antibiotics"))


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
mymod_ONLYcovar <- lme4::glmer(as.formula("diagnosis ~ site_name + consent_age + sex + race + Antibiotics + (1|Participant_ID)"), 
                               data = df_bestglm, 
                               family = binomial)

# predict based on only covariates
null_model_predictions <- as.data.frame(cbind(predictionDF$diagnosis, scale(predict(mymod_ONLYcovar, predictionDF, allow.new.levels = T))))
colnames(null_model_predictions) <- c("actual", "predicted")




# make a violin plot of the prediction
boxViolinPlot <- function(auc_df = avg_par_scores, covars = "", covars_only=F){
  
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
print("MODEL WITH ONLY FEATURES, NO COVARIATES")
PredPlot <- boxViolinPlot(auc_df = avg_par_scores, covars = "", covars_only=F)
PredPlot

ggsave("pred_features.png", width=2.5, height=2.5, units="in", dpi=320)


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

ggsave("pred_null.png", width=2.5, height=2.5, units="in", dpi=320)


pred_df
metatranscriptomics_pred_score_featuresonly <- tibble::rownames_to_column(pred_df, "External_ID")
write.table(metatranscriptomics_pred_score_featuresonly, "metatranscriptomics_features_scores.txt", sep="\t", col.names=T, row.names=F, quote=F)

# remove the intercept
featureplot_df <- mod_coef_df_nocovar[2:nrow(featureplot_df),]
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
      size = length(unique(df2$variable))*1/3,
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


# ## Mixed Effects Random Forests via MixRF ######################################################

myX <- as.data.frame(df_bestglm[,names(df_bestglm) %ni% c("diagnosis", "Participant_ID")])

myMixRF <- MixRF::MixRF(
  Y = as.numeric(pull(df_bestglm, diagnosis)),
  X = myX,
  random = "(1|Participant_ID)",
  MaxIterations = 10,
  data = df_bestglm
)

# ## RF Variable Importance ######################################################
caret::varImp(myMixRF$forest)
df_tmp <- caret::varImp(myMixRF$forest)
newdf <- df_tmp; colnames(newdf) <- c("imp")
# newdf <- data.frame(cbind(df_tmp$Overall)); colnames(newdf) <- c("imp"); rownames(newdf) <- rownames(df_tmp)
newdf$imp <- as.numeric(newdf$imp)
df2 <- newdf %>%
  tibble::rownames_to_column() %>%
  dplyr::rename("variable" = rowname) %>%
  dplyr::arrange(imp) %>%
  dplyr::mutate(variable = forcats::fct_inorder(variable))
plot_varimp2 <- ggplot2::ggplot(df2) +
  geom_segment(
    aes(
      x = variable,
      y = 0,
      xend = variable,
      yend = imp
    ),
    size = 1.5,
    alpha = 0.7
  ) +
  geom_point(aes(x = variable, y = imp, col = variable),
             size = 4,
             show.legend = F) +
  coord_flip() +
  labs(y = "Importance", x = NULL, title = "") +
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
      size = length(unique(df2$variable))*1/3,
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


levels(mergeytest[["sex"]])
levels(myX[["sex"]])
levels(mergeytest[["site_name"]])
levels(myX[["site_name"]])
levels(mergeytest[["race"]])
levels(myX[["race"]])

# ## RF Prediction ######################################################
# filter prediction dataframe to only complete cases
predictionDF <- mergeytest[complete.cases(mergeytest),]
predictionDF$consent_age <- as.factor(as.character(predictionDF$consent_age))
common <- names(myX)[names(myX) %in% names(predictionDF)]
predictionDF[common] <- lapply(common, function(x) {
  match.fun(paste0("as.", class(myX[[x]])))(predictionDF[[x]])
})
str(predictionDF[,common])
str(myX[,common])

# make the factors match between testing and training
myX$diagnosis <- as.numeric(pull(df_bestglm, diagnosis))
predictionDF <- predictionDF[,c(common,"diagnosis")]
predictionDF <- rbind(myX[1, ] , predictionDF)
predictionDF <- predictionDF[-1,]

# predictionDF$diagnosis <- as.factor(predictionDF$diagnosis)
postPred_RF <- as.data.frame(scale(predict(myMixRF$forest, newdata = predictionDF)))
# postPred_RF <- as.data.frame(predict(myglm, newdata = mergeytest)); postPred_RF <- ifelse(postPred_RF[,1] > 0.25, "1", "0")

comparePrediction <- cbind(as.character(predictionDF$diagnosis), postPred_RF)
comparePrediction <- as.data.frame(comparePrediction)
colnames(comparePrediction) <- c("actual", "predicted")
(comparePrediction$actual == comparePrediction$prediction)
summary((comparePrediction$actual == comparePrediction$prediction))
# View(comparePrediction)

PredPlotRF <- boxViolinPlot(pred_df = comparePrediction, predictionDF = predictionDF)
PredPlotRF




metatranscriptomics_glmer_scores <- pred_df
metatranscriptomics_glmer_scores <- tibble::rownames_to_column(metatranscriptomics_glmer_scores, "External_ID")
write.table(metatranscriptomics_glmer_scores, "metatranscriptomics_glmer_scores.txt", sep="\t", col.names=T, row.names=F, quote=F)
metatranscriptomics_MixRF_scores <- comparePrediction
metatranscriptomics_MixRF_scores <- tibble::rownames_to_column(metatranscriptomics_MixRF_scores, "External_ID")
write.table(metatranscriptomics_MixRF_scores, "metatranscriptomics_MixRF_scores.txt", sep="\t", col.names=T, row.names=F, quote=F)




# #---make a regression tree---#########################################################################################################################################################################
# 
# require(tree)
# mytree <- tree(diagnosis ~ ., data = training_ready_sub_vsurf_result)
# plot(mytree)
# text(mytree, pretty = 0, cex = .8)
