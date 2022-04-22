
setwd("/Users/chris/Documents/GRADSCHOOL/PolyOmicsRotation")

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
library(biomformat)
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
metadata <- subset(metadata, data_type == "metagenomics")
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
training_metadata <- fread("training_metadata.txt", sep = "\t")
training_metadata <- subset(training_metadata, data_type == "metagenomics")
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
testing_metadata <- fread("testing_metadata.txt", sep = "\t")
testing_metadata <- subset(testing_metadata, data_type == "metagenomics")
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
fname <- fread("taxonomic_profiles_3.tsv.gz", header=T)
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
mygrep <- grep("UNKNOWN", rownames(df_3), invert = F)
grep_species <- grep("\\|s__", rownames(df_3), invert = F)
mygrepni <- which(1:length(rownames(df_3)) %ni% mygrep)
if(length(commentRow) > 0){
  mygrepni <- mygrepni[which(mygrepni != commentRow)]
}

## Preprocessing, make data compositional for species --###################################
# grouped_df_3 <- df_3[mygrepni,]
# rownames(grouped_df_3) <- rownames(df_3)[mygrepni]
grouped_df_3 <- df_3[grep_species,]
rownames(grouped_df_3) <- rownames(df_3)[grep_species]
num_grouped_df_3 <- mutate_all(grouped_df_3, function(x) as.numeric(as.character(x)))

look <- as.data.frame(rownames(num_grouped_df_3)); colnames(look) <- c("look")
looksep <- look %>% separate(look,into=c("kingdom","phylum","class","order","family","genus","species"),convert=TRUE,sep="\\|")

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
namesiessep <- namesies %>% separate(namesies,into=c("namesies","junk"),convert=TRUE,sep="_profi")
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
mergeytest_colsd <- apply(mergeytest[complete.cases(mergeytest),grep("s__", names(mergeytest), invert = F)], 2, sd, na.rm=T)
mergey_colsd <- apply(mergey[complete.cases(mergey),grep("s__", names(mergey), invert = F)], 2, sd, na.rm=T)
# qthresh <- quantile(colsd, 0.05, na.rm=T)
hist(as.numeric(mergeytest_colsd))
hist(as.numeric(mergey_colsd))
toKeep <- c(intersect(c(names(mergeytest_colsd)[which(mergeytest_colsd > 0.1)]), c(names(mergey_colsd)[which(mergey_colsd > 0.1)])),
            c("External ID", "Participant ID", "site_name", "diagnosis", "consent_age", "sex", "race", "Antibiotics")
            )
length(toKeep)
length(mergey_colsd)
mergeytest <- mergeytest[,which(names(mergeytest) %in% toKeep)]
mergey <- mergey[,which(names(mergey) %in% toKeep)]

# rename the column names with just the species
cn <- as.data.frame(colnames(mergey)); colnames(cn) <- c("cn")
cn <- cn %>% separate(cn,into=c("junk","species"),convert=TRUE,sep="\\|s__")
cn <- cn$species
cn[is.na(cn)] <- c("Participant_ID", "site_name", "diagnosis", "consent_age", "sex", "race", "Antibiotics")
head(cn,20)
tail(cn,20)
# remove the square brackets
cn <- gsub(pattern = "\\[", replacement = "", cn)
cn <- gsub(pattern = "\\]", replacement = "", cn)
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

str(mergey[,1:30])
str(mergey[,(ncol(mergey)-30):ncol(mergey)])

traindf <- as.data.frame(mergey[complete.cases(mergey),])
dim(traindf)
dim(mergey)
varlist <- cn[which(cn %ni% c("diagnosis", "Participant_ID"))]
# varstring <- paste0(varlist[sample(1:length(varlist),200)], collapse = " + ", sep = "")
varstring <- paste0(varlist, collapse = " + ", sep = "")

## LASSO Lambda Search ######################################################

numvariables <- c()
lambdavec <- seq(from = 70, to = 110, by = 5)
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

lm1 <- glmmLasso(as.formula(paste0("diagnosis ~ ",varstring)),
                 data = traindf, 
                 rnd = list(Participant_ID=~1),
                 lambda=100,
                 family = binomial(link = "logit"))
summary(lm1)
lassoFeatures <- names(lm1$coefficients[which(lm1$coefficients != 0)])
# lassoFeatures <- lassoFeatures[lassoFeatures %ni% c("(Intercept)", "sexMale")]
lassoFeatures <- lassoFeatures[grep("Participant_ID|site_name|diagnosis|consent_age|sex|race|Antibiotics|Intercept", lassoFeatures, invert = T)]
# lassoFeatures <- unique(c(lassoFeatures, "Participant_ID", "site_name", "diagnosis", "consent_age", "sex", "race", "Antibiotics"))
lassoFeatures <- unique(c(lassoFeatures, "Participant_ID", "site_name", "diagnosis", "consent_age", "sex", "race"))

## Prediction Boxplot and AUC for LASSO ######################################################
df_bestglm <- as.data.frame(traindf[,c(lassoFeatures)])
df_bestglm$diagnosis <- as.factor(as.character(df_bestglm$diagnosis))
summary(df_bestglm$diagnosis)

varlist2 <- names(df_bestglm)[which(names(df_bestglm) %ni% c("diagnosis", "Participant_ID"))]
# varstring <- paste0(varlist[sample(1:length(varlist),200)], collapse = " + ", sep = "")
varstring2 <- paste0(varlist2, collapse = " + ", sep = "")

mymod <- lme4::glmer(as.formula(paste0("diagnosis ~ ",varstring2, " + (1|Participant_ID)")), 
             data = df_bestglm, 
             family = binomial)
mymodsum <- summary(mymod)
# prediction on reserved validation samples
# filter prediction dataframe to only complete cases
predictionDF <- mergeytest[complete.cases(mergeytest),]
pred_df <- as.data.frame(cbind(predictionDF$diagnosis, scale(predict(mymod, predictionDF, allow.new.levels = T)))); colnames(pred_df) <- c("actual", "predicted")
pred_df$actual <- as.factor(pred_df$actual)
# make a violin plot of the prediction
boxViolinPlot <- function(pred_df = pred_df, predictionDF = predictionDF){
  PredPlot <- ggplot(data = pred_df, aes(x = actual, y = predicted))+
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
  
  pred_df$actual <- as.factor(pred_df$actual)
  pred_df$predicted <- as.numeric(pred_df$predicted)
  caseControlGLM <- glm(as.formula(pred_df$actual ~ pred_df$predicted), data = pred_df, family = "binomial", na.action = na.omit)
  predpr <- predict(caseControlGLM, predictionDF, allow.new.levels = T, type = c("response"))
  caseControlroccurve <- pROC::roc(predictionDF$diagnosis ~ predpr, quiet=T, plot=F)
  caseControlroccurveCI <- pROC::roc(predictionDF$diagnosis ~ predpr, ci=T, quiet=T)
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
  minPRS <- min(pred_df$predicted)-(abs(max(pred_df$predicted)-min(pred_df$predicted)))*(label_disty*1.5*1.02)
  maxPRS <- max(pred_df$predicted) + (abs(max(pred_df$predicted)-min(pred_df$predicted)))*(label_disty*1.5*1.02)
  botLabLoc <- min(pred_df$predicted)-(abs(max(pred_df$predicted)-min(pred_df$predicted)))*(label_disty*1.5)
  topLabLoc <- max(pred_df$predicted) + (abs(max(pred_df$predicted)-min(pred_df$predicted)))*(label_disty*0.21)
  
  PredPlot <- PredPlot +
    geom_signif(textsize = 2.25, comparisons = list(c("0", "1")), annotations=caseControlAUCsummary, color="black", y_position = topLabLoc,tip_length=.03)+
    geom_signif(textsize = 2.25, comparisons = list(c("0", "1")), annotations=caseControlpthresh, color="black", y_position = botLabLoc,tip_length=-.03) +
    scale_y_continuous(breaks = seq(-100,100, by=2), limits = c(minPRS,maxPRS))
  return(PredPlot)
}

PredPlotLasso <- boxViolinPlot(pred_df = pred_df, predictionDF = predictionDF)
PredPlotLasso

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

# #---make a regression tree---#########################################################################################################################################################################
# 
# require(tree)
# mytree <- tree(diagnosis ~ ., data = training_ready_sub_vsurf_result)
# plot(mytree)
# text(mytree, pretty = 0, cex = .8)
