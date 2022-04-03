setwd("C:\\Users\\rgarr\\Documents\\poly-omics-risk")

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
`%ni%` <- Negate(`%in%`)

## metadata --###################################
metadata <- fread("https://ibdmdb.org/tunnel/products/HMP2/Metadata/hmp2_metadata.csv", header=T, stringsAsFactors=T)
str(metadata)


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
metadata <- subset(metadata, data_type == "viromics")
dim(metadata)
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

## Viromics taxonomic_profiles --###################################
# #MetaPhlAn2
# fname <- fread("https://ibdmdb.org/tunnel/products/HMP2/Viromics/1732/taxonomic_profiles.tsv.gz", header=T)
# fname <- fread("taxonomic_profiles_3.tsv.gz", header=T)
# fname <- as.data.frame(fname)

#VirMap
fname <- fread("https://ibdmdb.org/tunnel/products/HMP2/Viromics/1732/virome_virmap_analysis.tsv.gz", header=T)
fname <- fread("virome_virmap_analysis.tsv.gz", header=T)
fname <- as.data.frame(fname)

df_3 <- fname
#View(df_3)
print("data dimensions:")
print(dim(df_3))
print("data glimpse:")
str(df_3[1:5,1:4])
rownames(df_3) <- df_3$`Virus`
df_3$`Virus` <- NULL
View(as.data.frame(rownames(df_3)))
# separate dataframe out by rows which were annotated vs rows which weren't
commentRow <- grep("# ", rownames(df_3), invert = F)
mygrep <- grep("UNKNOWN", rownames(df_3), invert = F)
grep_species <- grep("\\|s__", rownames(df_3), invert = F)
mygrepni <- which(1:length(rownames(df_3)) %ni% mygrep)
if(length(commentRow) > 0){
  mygrepni <- mygrepni[which(mygrepni != commentRow)]
}


#######NEED TO FIX######
## Preprocessing, make data compositional for species --###################################
# grouped_df_3 <- df_3[mygrepni,]
# rownames(grouped_df_3) <- rownames(df_3)[mygrepni]
# grouped_df_3 <- df_3[grep_species,]
# rownames(grouped_df_3) <- rownames(df_3)[grep_species]
# num_grouped_df_3 <- mutate_all(grouped_df_3, function(x) as.numeric(as.character(x)))
# View(num_grouped_df_3)
# 
# look <- as.data.frame(rownames(num_grouped_df_3)); colnames(look) <- c("look")
# looksep <- look %>% separate(look,into=c("kingdom","phylum","class","order","family","genus","species"),convert=TRUE,sep="\\|")



# is it compositional?
summary(colSums(num_grouped_df_3))
hist(colSums(num_grouped_df_3))

# add epsilon to all entries to set up for center log transform
pepsi <- 1E-06
num_grouped_df_3 <- num_grouped_df_3 + pepsi
View(num_grouped_df_3)


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


