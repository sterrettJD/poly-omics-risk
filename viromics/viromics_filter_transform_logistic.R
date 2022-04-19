setwd("C:\\Users\\rgarr\\Documents\\poly-omics-risk")
#setwd("/Users/johnsterrett/Research-Projects/Team-rotation/poly-omics-scores/viromics/")


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
library(compositions)
library(lme4)
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

#VirMap
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

# add epsilon to all entries to set up for center log transform
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

# rename the columns for the merge with diagnosis
namesies <- as.data.frame(colnames(clr_num_grouped_df_3)); colnames(namesies) <- c("namesies")
head(namesies)
# namesiessep <- namesies %>% separate(namesies,into=c("namesies","junk"),convert=TRUE,sep="_profi")
namesiessep <- namesies
colnames(clr_num_grouped_df_3) <- namesiessep$namesies

intersecty <- intersect(c(as.character(namesiessep$namesies)), 
                        c(as.character(metadata1$`External ID`))
)

View(intersecty)
dim(namesiessep)
dim(metadata1)
length(intersecty)



# transpose transformed data
transp_clr_num_grouped_df_3 <- as.data.frame(t(clr_num_grouped_df_3))
transp_clr_num_grouped_df_3$`External ID` <- rownames(transp_clr_num_grouped_df_3)

## Merge with metadata --###################################
#read in training set and merge with viromics data
train_meta= fread("training_metadata.txt", sep='\t', header=FALSE)
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

#view(train_meta)


#read in testing metadata and merge with viromics
test_meta= fread("testing_metadata.txt", sep='\t', header=FALSE)
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

#view(test_meta)




# merge with diagnosis
# mergey <- merge(transp_clr_num_grouped_df_3, metadata1, by = "External ID")
# dim(transp_clr_num_grouped_df_3)
# dim(mergey)
# # make validation dataset that isn't used for training (the ids that are in the polyomic list)
# mergeytest <- mergey[which(as.character(mergey$`External ID`) %in% as.character(idlist)),]
# dim(mergeytest)
# # make training dataset (the ids that are NOT in the polyomic list)
# mergey <- subset(mergey, `External ID` %ni% idlist)
# dim(mergey)
# # make the ids the rownames for each dataframe, and then remove that column
# rownames(mergeytest) <- mergeytest$`External ID`
# mergeytest$`External ID` <- NULL
# rownames(mergey) <- mergey$`External ID`
# mergey$`External ID` <- NULL


# remove any columns that have no/little variation between samples...
mergeytest_colsd <- apply(mergeytest, 2, sd, na.rm=T)
mergey_colsd <- apply(mergey, 2, sd, na.rm=T)

# qthresh <- quantile(colsd, 0.05, na.rm=T)
hist(as.numeric(mergeytest_colsd))
hist(as.numeric(mergey_colsd))
toKeep <- intersect(c(names(mergeytest_colsd)[which(mergeytest_colsd > 0)]), c(names(mergey_colsd)[which(mergey_colsd > 0)]))
length(toKeep)
length(mergey_colsd)
mergeytest <- mergeytest[,which(names(mergeytest) %in% toKeep)]
mergey <- mergey[,which(names(mergey) %in% toKeep)]

# rename the column names with just the species
# cn <- as.data.frame(colnames(mergey)); colnames(cn) <- c("cn")
# cn <- cn %>% separate(cn,into=c("junk","species"),convert=TRUE,sep="species=")
# cn <- cn$species
# cn[length(cn)] <- "diagnosis"
# head(cn)
# tail(cn)
# colnames(mergeytest) <- cn
# colnames(mergey) <- cn

#combine mergys with metadata
mergey <- merge(transp_clr_num_grouped_df_3, train_meta, by = "External ID")
mergeytest <- merge(transp_clr_num_grouped_df_3, test_meta, by = "External ID")

## run logistic regression for each feature --###################################
featureNames <- c()
betas <- c()
pvals <- c()
for(i in 1:(ncol(mergey)-1)){
  # for(i in 1:1){
  # randName <- names(mergey)[sample(1:length(names(mergey)),1)]
  randName <- names(mergey)[i]
  mergeysub <- mergey[,c(randName, "Participant ID", "race", 
                         "consent_age", "site_name", 
                         "sex", "Antibiotics", "diagnosis")]
  
  colnames(mergeysub) <- c(randName, "Participant ID", "race", 
                           "consent_age", "site_name", 
                           "sex", "Antibiotics", "diagnosis")
  mergeysub$`Participant ID` <- as.factor(mergeysub$`Participant ID`)
  
  mymod <- glmer(as.formula(
    paste0("diagnosis ~ `", randName,"` + consent_age + sex + race + Antibiotics + (1|`Participant ID`) + (1|`site_name`)")), 
    data = mergeysub, family = "binomial", 
    #control = glmerControl(optimizer ="Nelder_Mead")
    control = glmerControl(optimizer ="bobyqa")
  )
  #mymod <- glm(as.formula(paste0("diagnosis ~ `", randName,"` + consent_age + sex + race + Antibiotics + (1|`Participant ID`) + (1|`site_name`)")), data = mergeysub, family = "binomial")
  #mymod <- glm(as.formula(paste0("diagnosis ~ `",randName,"`")), data = mergeysub, family = "binomial")
  mymodsum <- summary(mymod)
  featureNames <- c(featureNames, randName)
  betas <- c(betas, mymodsum$coefficients[2,1])
  pvals <- c(pvals, mymodsum$coefficients[2,4])
  # validate on testing set
  # if(mymodsum$coefficients[2,4] < 0.05){
  #   print(randName)
  #   print(cor(mergeytest$diagnosis, predict(mymod, mergeytest)))
  # }
  # Sys.sleep(1)
  # summary(mymod)
  # print(randName)
  # print(mymodsum$aic)
}




# multiple test p value for significance (bonferoni)
bonfsigthresh <- 0.05/(ncol(mergey)-1)
# make dataframe
tenResults <- as.data.frame(cbind(featureNames,betas,pvals))
tenResults$featureNames <- as.character(tenResults$featureNames)
tenResults$betas <- as.numeric(as.character(tenResults$betas))
tenResults$pvals <- as.numeric(as.character(tenResults$pvals))


# column for plot colors
tenResults$mycolors <- NA
tenResults$mycolors[tenResults$pvals < bonfsigthresh] <- "sig"
tenResults$mycolors[tenResults$pvals >= bonfsigthresh] <- "notsig"
tenResults$mycolors <- as.factor(tenResults$mycolors)


# column for plot labels
tenResults$delabels <- NA
tenResults$delabels[tenResults$pvals < bonfsigthresh] <- tenResults$featureNames[tenResults$pvals < bonfsigthresh]
tenResults$delabels[tenResults$pvals >= bonfsigthresh] <- NA
tenResults$delabels <- as.character(tenResults$delabels)

# make volcano plot
bplot <- ggplot(aes(x = betas, y = -log10(pvals), col=mycolors, label=delabels), data = tenResults) +
  # geom_bar(stat="identity", fill = "steelblue") + theme_minimal()
  geom_point() +
  theme_minimal() +
  geom_text()
bplot



# NO BONFERONNI
# column for plot colors
tenResults$significant <- NA
tenResults$significant[tenResults$pvals < 0.05] <- "sig"
tenResults$significant[tenResults$pvals >= 0.05] <- "notsig"
tenResults$significant <- as.factor(tenResults$significant)


# column for plot labels
tenResults$labels <- NA
tenResults$labels[tenResults$pvals < 0.05] <- tenResults$featureNames[tenResults$pvals < 0.05]
tenResults$labels[tenResults$pvals >= 0.05] <- NA
tenResults$labels <- as.character(tenResults$labels)

# make volcano plot
bplot <- ggplot(aes(x = betas, y = -log10(pvals), col=significant, label=labels), data = tenResults) +
  # geom_bar(stat="identity", fill = "steelblue") + theme_minimal()
  geom_point() +
  theme_minimal() +
  geom_text()
bplot



# extract the features that were significant and run a glm on the full model
sigFeatures <- tenResults$featureNames[tenResults$mycolors == "sig"]
df_bestglm <- mergey[,c(sigFeatures,"diagnosis")]
mymod <- glm(as.formula(paste0("diagnosis ~ .")), data = df_bestglm, family = "binomial")
mymodsum <- summary(mymod)
# prediction on reserved validation samples
pred_df <- as.data.frame(cbind(mergeytest$diagnosis, predict(mymod, mergeytest))); colnames(pred_df) <- c("actual", "predicted")
pred_df$actual <- as.factor(pred_df$actual)
# make a violin plot of the prediction
ggplot(data = pred_df, aes(x = actual, y = predicted))+
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
  ylab("Predicted Diagnosis")+
  xlab("Actual Diagnosis")


# bestglm
# function begins with a data frame containing explanatory variables and response variables. 
# The response variable should be in the last column. 
# Varieties of goodness-of-fit criteria can be specified in the IC argument.
# bestglm can only handle 15 variables.
subfeats <- sample(names(df_bestglm), 14)
subfeats <- subfeats[which(subfeats != "diagnosis")]; subfeats <- c(subfeats, "diagnosis")
bglm <- bestglm::bestglm(df_bestglm[,subfeats], family = binomial, IC = "BIC")
# prediction on reserved validation samples
pred_df <- as.data.frame(cbind(mergeytest$diagnosis, predict(bglm$BestModel, mergeytest))); colnames(pred_df) <- c("actual", "predicted")
pred_df$actual <- as.factor(pred_df$actual)
# make a violin plot of the prediction
ggplot(data = pred_df, aes(x = actual, y = predicted))+
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
  ylab("Predicted Diagnosis")+
  xlab("Actual Diagnosis")
