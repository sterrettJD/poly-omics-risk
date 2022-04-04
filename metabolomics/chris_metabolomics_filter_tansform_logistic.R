
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
`%ni%` <- Negate(`%in%`)

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
metadata <- subset(metadata, data_type == "metabolomics")
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

## Metabolomics HMP2_metabolomics --###################################
# fname <- fread("https://ibdmdb.org/tunnel/products/HMP2/WGS/1818/taxonomic_profiles_3.tsv.gz", header=T)
fname <- fread("https://ibdmdb.org/tunnel/products/HMP2/Metabolites/1723/HMP2_metabolomics.csv.gz", header=T)
# fname <- fread("pathabundances_3.tsv.gz", header=T)
fname <- as.data.frame(fname)
# data_glimpse(fname = fname, UNINTEGRATED = F, thresh = 0.8, normalize = T, myTitle = "Metatranscriptomes pathabundances_3")

df_3 <- fname
df_3 <- subset(df_3, Metabolite != "")
# df_3$Method <- NULL
# df_3$`Pooled QC sample CV` <- NULL
# df_3$`m/z` <- NULL
# df_3$RT <- NULL
# df_3$`HMDB (*Representative ID)` <- NULL
# df_3$Metabolite <- NULL
df_3 <- df_3 %>% dplyr::select(-Method, -`Pooled QC sample CV`, -`m/z`, -RT, -`HMDB (*Representative ID)`, -Metabolite)
rownames(df_3) <- df_3$Compound
df_3$Compound <- NULL
df_3[is.na(df_3) == TRUE] <- 0
dim(df_3)
hist(colSums(df_3))
hist(rowSums(df_3))

alldf <- as.data.frame(as.numeric(array(as.matrix(df_3)))); colnames(alldf) <- c("alldf")
hist(alldf$alldf[alldf$alldf])

num_grouped_df_3 <- df_3

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

# rename the columns for the merge with diagnosis
namesies <- as.data.frame(colnames(clr_num_grouped_df_3)); colnames(namesies) <- c("namesies")
head(namesies)
# namesies$namesies[grep("CS", namesies$namesies)] <- paste0(namesies$namesies[grep("CS", namesies$namesies)], "_P")
namesiessep <- namesies %>% separate(namesies,into=c("namesies","junk"),convert=TRUE,sep="\\_")
namesiessep <- namesies
colnames(clr_num_grouped_df_3) <- namesiessep$namesies

# View(cbind(c(as.character(namesiessep$namesies)), 
#                c(as.character(metadata1$`External ID`))))

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
mergey <- merge(transp_clr_num_grouped_df_3, metadata1, by = "External ID")
dim(transp_clr_num_grouped_df_3)
dim(mergey)
# make validation dataset that isn't used for training (the ids that are in the polyomic list)
mergeytest <- mergey[which(as.character(mergey$`External ID`) %in% as.character(idlist)),]
dim(mergeytest)
# make training dataset (the ids that are NOT in the polyomic list)
mergey <- subset(mergey, `External ID` %ni% idlist)
dim(mergey)
# make the ids the rownames for each dataframe, and then remove that column
rownames(mergeytest) <- mergeytest$`External ID`
mergeytest$`External ID` <- NULL
rownames(mergey) <- mergey$`External ID`
mergey$`External ID` <- NULL

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
# cn <- cn %>% separate(cn,into=c("junk","species"),convert=TRUE,sep="\\|s__")
# cn <- cn$species
# cn[length(cn)] <- "diagnosis"
# head(cn)
# tail(cn)
# colnames(mergeytest) <- cn
# colnames(mergey) <- cn

head(colnames(mergeytest))
head(colnames(mergey))

## run logistic regression for each feature --###################################
featureNames <- c()
betas <- c()
pvals <- c()
for(i in 1:(ncol(mergey)-1)){
# for(i in 1:1){
  # randName <- names(mergey)[sample(1:length(names(mergey)),1)]
  randName <- names(mergey)[i]
  mergeysub <- mergey[,c(randName, "diagnosis")]
  colnames(mergeysub) <- c(randName,"diagnosis")
  mymod <- glm(as.formula(paste0("diagnosis ~ `",randName,"`")), data = mergeysub, family = "binomial")
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
# bonfsigthresh <- 0.05/(ncol(mergey)-1)
bonfsigthresh <- 0.05

# make dataframe
tenResults <- as.data.frame(cbind(featureNames,betas,pvals))
tenResults$featureNames <- as.character(tenResults$featureNames)
tenResults$betas <- as.numeric(tenResults$betas)
tenResults$pvals <- as.numeric(tenResults$pvals)
# tenResults <- tenResults %>% separate(featureNames,into=c("kingdom","phylum","class","order","family","genus","species"),convert=TRUE,sep="\\|")

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


# stepAIC: Choose a model by AIC in a Stepwise Algorithm
# MASS::stepAIC(df_bestglm, scope = as.formula(paste0("diagnosis ~ .")), direction = c("both"))








## RANDOM FORESTS USING VSURF ######################################################
# User code starts here
stopifnot(require(VSURF))
stopifnot(require(tidyverse))
stopifnot(require(doParallel))
stopifnot(require(doFuture))
stopifnot(require(caret))
stopifnot(require(e1071))
stopifnot(require(randomForest))
stopifnot(require(magrittr))
stopifnot(require(dplyr))

unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
'%ni%' <- Negate('%in%')

#---load in things ---#########################################################################################################################################################################
data <- df_bestglm
data_testing <- mergeytest

data <- data[complete.cases(data),]

# R CMD BATCH --vanilla '--args fss' training.R output/training_diagnosis.Rout

#---RF setup ---#########################################################################################################################################################################

outcomeVariable <- "diagnosis"

num_cores = detectCores() - 1
print(
  paste0(
    'Number of cores being used = ',
    num_cores,
    ", of possible ",
    detectCores(),
    " cores"
  )
)
# registerDoParallel(num_cores)

# registerDoFuture()
# plan(multiprocess, workers = num_cores)

x <- data %>% dplyr::select(-diagnosis) #%>% dplyr::select(hosplos, gcs_use, admittoicudc1)
y <- data %>% dplyr::select(diagnosis)

dim(data)
x <- as.data.frame(x)
y <- y$diagnosis

number_trees = 2000
mtry_best <- max(floor(ncol(x) / 3), 1)
# mtry_best <- max(floor(ncol(x) / 2), 1)

# bestMtry <- randomForest::tuneRF(
#   x,
#   y,
#   stepFactor = 1,
#   improve = 1e-5,
#   ntree = number_trees,
#   na.action = na.omit
# )
# bestMtry <- as.data.frame(bestMtry)
# mtry_best <- bestMtry$mtry[which.min(bestMtry$OOBError)]

# # data(iris)
# # x <- iris[,1:4]
# # y <- iris[,5]
#
# # vozone <- VSURF(diagnosis ~ ., data = data, na.action = na.omit)
# # summary(vozone)
#
# #---THRESHOLDING ---#########################################################################################################################################################################
results.vsurf_thresh <- VSURF::VSURF_thres(
  x = x,
  y = y,
  na.action = na.omit,
  mtry = mtry_best,
  ntree = number_trees,
  parallel = TRUE,
  verbose = TRUE,
  ncores = num_cores,
  # nmj = 1,
  nfor.thres = 50,
  nmin = 1
)
length(results.vsurf_thresh$varselect.thres)
# save.image(
#   file = paste0("thresh_save",
#                 "_VSURFkeepers.Rdata")
# )
# load(paste0("./InputData/",
#             "thresh_save",
#             "_VSURFkeepers.Rdata"))
#---INTERPRETATION ---#########################################################################################################################################################################
results.vsurf_interp <- VSURF::VSURF_interp(
  x,
  y,
  na.action = na.omit,
  vars =  results.vsurf_thresh$varselect.thres,
  ntree = number_trees,
  parallel = TRUE,
  verbose = TRUE,
  ncores = num_cores,
  nfor.interp = 25, #25,
  nsd = 1
)
length(results.vsurf_interp$varselect.interp)
# save.image(
#   file = paste0("interp_save",
#                 "_VSURFkeepers.Rdata")
# )
#---PREDICTION ---#########################################################################################################################################################################
results.vsurf_pred <- VSURF::VSURF_pred(
  x,
  y,
  na.action = na.omit,
  err.interp = results.vsurf_interp$err.interp,
  varselect.interp = results.vsurf_interp$varselect.interp,
  ntree = number_trees,
  parallel = TRUE,
  verbose = TRUE,
  ncores = num_cores,
  nfor.pred = 25,
  nmj = 1 # increasing the value of nmj leads to selection of fewer variables:
)
length(results.vsurf_pred$varselect.pred)
# save.image(
#   file = paste0("pred_save",
#                 "_VSURFkeepers.Rdata")
# )

results.vsurf <- results.vsurf_pred
results.vsurf.OG <- results.vsurf

nVarInterp <-
  length(colnames(data[, results.vsurf_interp$varselect.interp]))
nVarPred <-
  length(colnames(data[, results.vsurf_pred$varselect.pred]))

# look at results of VSURF
summary(results.vsurf)
plot(results.vsurf_thresh)
plot(results.vsurf_interp)
plot(results.vsurf_pred)
results.vsurf_thresh$varselect.thres
results.vsurf_interp$varselect.interp
results.vsurf_pred$varselect.pred

# print the reduced number of variables that should be considered in model
colnames(x)[results.vsurf_thresh$varselect.thres]
colnames(x)[results.vsurf_interp$varselect.interp]
colnames(x)[results.vsurf_pred$varselect.pred] # The final list of variables to be included according to the VSURF methodology.
VSURF_thres_keepers <-
  colnames(x)[results.vsurf_thresh$varselect.thres]
VSURF_interp_keepers <-
  colnames(x)[results.vsurf_interp$varselect.interp]
VSURF_pred_keepers <-
  colnames(x)[results.vsurf_pred$varselect.pred]

# VSURF_thres_keepers <- c(VSURF_thres_keepers, "age", "female")
# VSURF_interp_keepers <- c(VSURF_interp_keepers, "age", "female")
# VSURF_pred_keepers <- c(VSURF_pred_keepers, "age", "female")

# # if one of the steps only produced 1 variable, use the previous VSURF list of names
# if (length(VSURF_pred_keepers) > 1) {
#   training_ready_sub_vsurf_result = dplyr::select(data,
#                                                   c(outcomeVariable, VSURF_pred_keepers))
#   training_ready_sub_vsurf_result_varImp = dplyr::select(data,
#                                                          c(outcomeVariable, VSURF_interp_keepers))
# } else if (length(VSURF_pred_keepers) <= 1 &&
#            length(VSURF_interp_keepers) > 1) {
#   training_ready_sub_vsurf_result = dplyr::select(data,
#                                                   c(outcomeVariable, VSURF_interp_keepers))
#   training_ready_sub_vsurf_result_varImp = dplyr::select(data,
#                                                          c(outcomeVariable, VSURF_interp_keepers))
# } else if (length(VSURF_pred_keepers) <= 1 &&
#            length(VSURF_interp_keepers) <= 1) {
#   training_ready_sub_vsurf_result = dplyr::select(data,
#                                                   c(outcomeVariable, VSURF_thres_keepers))
#   training_ready_sub_vsurf_result_varImp = dplyr::select(data,
#                                                          c(outcomeVariable, VSURF_thres_keepers))
# }

training_ready_sub_vsurf_result = dplyr::select(data,
                                                c(outcomeVariable, VSURF_pred_keepers))
training_ready_sub_vsurf_result_varImp = dplyr::select(data,
                                                       c(outcomeVariable, VSURF_interp_keepers))
training_ready_sub_vsurf_result_tresh = dplyr::select(data,
                                                      c(outcomeVariable, VSURF_thres_keepers))

glimpse(training_ready_sub_vsurf_result)
glimpse(training_ready_sub_vsurf_result_varImp)
glimpse(training_ready_sub_vsurf_result_tresh)



# #---RF All Variables---#########################################################################################################################################################################

x <- x
# x <- x %>% dplyr::select(-diagnosis) 
y <- y

if (length(data) > 2) {
  bestMtry <-
    randomForest::tuneRF(
      x,
      y,
      stepFactor = 1,
      improve = 1e-5,
      ntree = number_trees,
      na.action = nasaction
    )
  bestMtry <- as.data.frame(bestMtry)
  mtry_best <- bestMtry$mtry[which.min(bestMtry$OOBError)]
  # mtry_best <- max(floor(ncol(x)/3), 1)
} else{
  mtry_best <- 1
}
tunegrid <- expand.grid(
  .mtry = mtry_best,
  .splitrule = c('gini'),
  .min.node.size = c(5, 10, 20)
)
control <- caret::trainControl(method = "cv",
                               number = 3,
                               # repeats=3,
                               # verboseIter = T,
                               # classProbs = T,
                               allowParallel = TRUE,
                               verboseIter = TRUE)
mod_formula <- as.formula(paste(outcomeVariable, "~", "."))
unregister_dopar()
rf.mod.ALL <- caret::train(
  mod_formula,
  data = data,
  # data = training_ready_sub_vsurf_result_varImp,
  # data = training_ready_sub_vsurf_result,
  # method = 'ranger',
  method = 'rf',
  na.action = na.omit,
  keep.inbag = TRUE,
  replace = TRUE,
  # importance = "permutation", #***
  trControl = control,
  num.threads = num_cores
)
varImp(rf.mod.ALL)

# #---RF thresholding step ---#########################################################################################################################################################################

x <- training_ready_sub_vsurf_result_tresh
x <- x %>% dplyr::select(-diagnosis) 
y <- y

if (length(data) > 2) {
  bestMtry <-
    randomForest::tuneRF(
      x,
      y,
      stepFactor = 1,
      improve = 1e-5,
      ntree = number_trees,
      na.action = nasaction
    )
  bestMtry <- as.data.frame(bestMtry)
  mtry_best <- bestMtry$mtry[which.min(bestMtry$OOBError)]
  # mtry_best <- max(floor(ncol(x)/3), 1)
} else{
  mtry_best <- 1
}
tunegrid <- expand.grid(
  .mtry = mtry_best,
  .splitrule = c('gini'),
  .min.node.size = c(5, 10, 20)
)
control <- caret::trainControl(method = "cv",
                               number = 3,
                               # repeats=3,
                               # verboseIter = T,
                               # classProbs = T,
                               allowParallel = TRUE,
                               verboseIter = TRUE)
mod_formula <- as.formula(paste(outcomeVariable, "~", "."))
unregister_dopar()
rf.mod.THRESH <- caret::train(
  mod_formula,
  data = training_ready_sub_vsurf_result_tresh,
  # method = 'ranger',
  method = 'rf',
  na.action = na.omit,
  keep.inbag = TRUE,
  replace = TRUE,
  # importance = "permutation", #***
  trControl = control,
  num.threads = num_cores
)
varImp(rf.mod.THRESH)

#---RF Interpretation step---#########################################################################################################################################################################

x <- training_ready_sub_vsurf_result_varImp
x <- x %>% dplyr::select(-diagnosis) 
y <- y

if (length(data) > 2) {
  bestMtry <-
    randomForest::tuneRF(
      x,
      y,
      stepFactor = 1,
      improve = 1e-5,
      ntree = number_trees,
      na.action = nasaction
    )
  bestMtry <- as.data.frame(bestMtry)
  mtry_best <- bestMtry$mtry[which.min(bestMtry$OOBError)]
  # mtry_best <- max(floor(ncol(x)/3), 1)
} else{
  mtry_best <- 1
}
tunegrid <- expand.grid(
  .mtry = mtry_best,
  .splitrule = c('gini'),
  .min.node.size = c(5, 10, 20)
)
control <- caret::trainControl(method = "cv",
                               number = 3,
                               # repeats=3,
                               # verboseIter = T,
                               # classProbs = T,
                               allowParallel = TRUE,
                               verboseIter = TRUE)
mod_formula <- as.formula(paste(outcomeVariable, "~", "."))
unregister_dopar()
rf.mod.INTERP <- caret::train(
  mod_formula,
  data = training_ready_sub_vsurf_result_varImp,
  # method = 'ranger',
  method = 'rf',
  na.action = na.omit,
  keep.inbag = TRUE,
  replace = TRUE,
  # importance = "permutation", #***
  trControl = control,
  num.threads = num_cores
)
varImp(rf.mod.INTERP)

df_tmp <- varImp(rf.mod.INTERP)
df_tmp2 <- as.data.frame(df_tmp$importance); colnames(df_tmp2) <- c("imp")
df2 <- df_tmp2 %>%
  tibble::rownames_to_column() %>%
  dplyr::rename("variable" = rowname) %>%
  dplyr::arrange(imp) %>%
  dplyr::mutate(variable = forcats::fct_inorder(variable))
plot_varimp <- ggplot2::ggplot(df2) +
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
      size = length(unique(df2$variable))*1/6,
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
plot_varimp

#---RF Prediction step---#########################################################################################################################################################################

x <- training_ready_sub_vsurf_result
x <- x %>% dplyr::select(-diagnosis) 
y <- y

if (length(data) > 2) {
  bestMtry <-
    randomForest::tuneRF(
      x,
      y,
      stepFactor = 1,
      improve = 1e-5,
      ntree = number_trees,
      na.action = nasaction
    )
  bestMtry <- as.data.frame(bestMtry)
  mtry_best <- bestMtry$mtry[which.min(bestMtry$OOBError)]
  # mtry_best <- max(floor(ncol(x)/3), 1)
} else{
  mtry_best <- 1
}
tunegrid <- expand.grid(
  .mtry = mtry_best,
  .splitrule = c('gini'),
  .min.node.size = c(5, 10, 20)
)
control <- caret::trainControl(method = "cv",
                               number = 3,
                               # repeats=3,
                               # verboseIter = T,
                               # classProbs = T,
                               allowParallel = TRUE,
                               verboseIter = TRUE)
mod_formula <- as.formula(paste(outcomeVariable, "~", "."))
unregister_dopar()
rf.mod.PRED <- caret::train(
  mod_formula,
  data = training_ready_sub_vsurf_result,
  # method = 'ranger',
  method = 'rf',
  na.action = na.omit,
  keep.inbag = TRUE,
  replace = TRUE,
  # importance = "permutation", #***
  trControl = control,
  num.threads = num_cores
)
varImp(rf.mod.PRED)

df_tmp <- varImp(rf.mod.PRED)
df_tmp2 <- as.data.frame(df_tmp$importance); colnames(df_tmp2) <- c("imp")
df2 <- df_tmp2 %>%
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

#---simple glm---#########################################################################################################################################################################

print(paste0("diagnosis ~ ", paste(VSURF_interp_keepers, sep = "", collapse = " + ")))
myglm <- glm(as.formula(paste0("diagnosis ~ ", paste(VSURF_interp_keepers, sep = "", collapse = " + "))),
             data = data, na.action = na.omit, family = binomial()
)

#---evaluate prediction---#########################################################################################################################################################################

# postPred_RF <- as.data.frame(predict(rf.mod.ALL, newdata = mergeytest))
# postPred_RF <- as.data.frame(predict(rf.mod.THRESH, newdata = mergeytest))
# postPred_RF <- as.data.frame(predict(rf.mod.INTERP, newdata = mergeytest))
postPred_RF <- as.data.frame(predict(rf.mod.PRED, newdata = mergeytest))
# postPred_RF <- as.data.frame(predict(myglm, newdata = mergeytest)); postPred_RF <- ifelse(postPred_RF[,1] > 0.25, "1", "0")

comparePrediction <- cbind(as.character(mergeytest$diagnosis), postPred_RF)
comparePrediction <- as.data.frame(comparePrediction)
colnames(comparePrediction) <- c("actual", "predicted")
(comparePrediction$actual == comparePrediction$prediction)
summary((comparePrediction$actual == comparePrediction$prediction))
# View(comparePrediction)

ggplot(data = comparePrediction, aes(x = actual, y = predicted))+
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

print(VSURF_thres_keepers)
print(VSURF_interp_keepers)
print(VSURF_pred_keepers)

# rtn <- rf.mod.PRED
rtn <- rf.mod.INTERP #the interpretation set of selected parameters seems to do the best 
# rtn <- rf.mod.ALL
# rtn <- myglm

#---Look at effect size directions---#########################################################################################################################################################################

sub_tenResults <- subset(tenResults, featureNames %in% VSURF_pred_keepers)

# make volcano plot
bplot2 <- ggplot(aes(x = betas, y = -log10(pvals), col=mycolors, label=delabels), data = sub_tenResults) +
  # geom_bar(stat="identity", fill = "steelblue") + theme_minimal()
  geom_point() +
  theme_minimal() +
  geom_text() +
  geom_vline(xintercept = 0)
bplot2

#---make a regression tree---#########################################################################################################################################################################

require(tree)
mytree <- tree(diagnosis ~ ., data = training_ready_sub_vsurf_result)
plot(mytree)
text(mytree, pretty = 0, cex = .8)
