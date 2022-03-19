
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
library(data.table)
library(ggplot2)
library(R.utils)
library(tidyverse)
library(UpSetR)
library(cowplot)
library(biomformat)
library(compositions)
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
metadata <- subset(metadata, data_type == "metagenomics")
# metadata <- subset(metadata, `External ID` %in% idlist)

# fill missinig enries with NA and view thte columns that are mostly non-NA
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
# separate dataframe out by rows which were annotatted vs rows which weren't
commentRow <- grep("# ", rownames(df_3), invert = F)
mygrep <- grep("UNKNOWN", rownames(df_3), invert = F)
grep_species <- grep("\\|s__", rownames(df_3), invert = F)
mygrepni <- which(1:length(rownames(df_3)) %ni% mygrep)
if(length(commentRow) > 0){
  mygrepni <- mygrepni[which(mygrepni != commentRow)]
}

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


# transform by center log transform across rows (features)
clr_num_grouped_df_3 <- compositions::clr(num_grouped_df_3,)
# clr_num_grouped_df_3 <- as.data.frame(apply(num_grouped_df_3, 1, clr))


alldf <- as.data.frame(as.numeric(array(as.matrix(clr_num_grouped_df_3)))); colnames(alldf) <- c("alldf")
hist(alldf$alldf)

# rename the columns for the merge with diagnosis
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
# rownames(transp_clr_num_grouped_df_3) <- colnames(clr_num_grouped_df_3)
# colnames(transp_clr_num_grouped_df_3) <- rownames(clr_num_grouped_df_3)
transp_clr_num_grouped_df_3$`External ID` <- rownames(transp_clr_num_grouped_df_3)

# merge with diagnosis
mergey <- merge(transp_clr_num_grouped_df_3, metadata1, by = "External ID")
dim(transp_clr_num_grouped_df_3)
dim(mergey)
mergeytest <- mergey[which(as.character(mergey$`External ID`) %in% as.character(idlist)),]
dim(mergeytest)
mergey <- subset(mergey, `External ID` %ni% idlist)
dim(mergey)
rownames(mergey) <- mergey$`External ID`
mergey$`External ID` <- NULL

# run logistic regression
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
  if(mymodsum$coefficients[2,4] < 0.05){
    print(randName)
    print(cor(mergeytest$diagnosis, predict(mymod, mergeytest)))
  }
  Sys.sleep(1)
  # summary(mymod)
  # print(randName)
  # print(mymodsum$aic)
}
tenResults <- as.data.frame(cbind(featureNames,betas,pvals))
tenResults$featureNames <- as.character(tenResults$featureNames)
tenResults$betas <- as.numeric(tenResults$betas)
tenResults$pvals <- as.numeric(tenResults$pvals)
tenResults <- tenResults %>% separate(featureNames,into=c("kingdom","phylum","class","order","family","genus","species"),convert=TRUE,sep="\\|")

bonfsigthresh <- 0.05/(ncol(mergey)-1)

tenResults$mycolors <- NA
tenResults$mycolors[tenResults$pvals < bonfsigthresh] <- "sig"
tenResults$mycolors[tenResults$pvals >= bonfsigthresh] <- "notsig"
tenResults$mycolors <- as.factor(tenResults$mycolors)

tenResults$delabels <- NA
tenResults$delabels[tenResults$pvals < bonfsigthresh] <- tenResults$species[tenResults$pvals < bonfsigthresh]
tenResults$delabels[tenResults$pvals >= bonfsigthresh] <- NA
tenResults$delabels <- as.character(tenResults$delabels)

bplot <- ggplot(aes(x = betas, y = -log10(pvals), col=mycolors, label=delabels), data = tenResults) +
  # geom_bar(stat="identity", fill = "steelblue") + theme_minimal()
  geom_point() +
  theme_minimal() +
  geom_text()
bplot



