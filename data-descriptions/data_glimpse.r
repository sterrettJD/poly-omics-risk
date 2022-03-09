setwd("/Users/chris/Documents/GRADSCHOOL/PolyOmicsRotation")

list.files()

# install.packages("data.table")
# install.packages("ggplot2")
# install.packages("R.utils")
# install.packages("tidyverse")
library(data.table)
library(ggplot2)
library(R.utils)
library(tidyverse)

c_make_hist <- function(df, whichCol){
  p1 <- ggplot(df, aes_string(x=whichCol)) + 
    geom_histogram(color="black", fill="salmon", alpha = 0.5, bins = 30) +
    geom_vline(aes(xintercept=mean(.data[[whichCol]])),
               color="blue", linetype="dashed", size=1)
  return(p1)
}

c_make_violin <- function(df, whichCol){
  p1 <- ggplot(data = df,aes_string(x = NA, y = whichCol))+
    scale_fill_viridis_d( option = "D")+
    theme_dark(base_size = 14)+
    geom_violin(fill = "gray70",alpha=0.4, position = position_dodge(width = .5),size=1,color="gray22",width=.5,lwd=.2) +
    geom_boxplot(fill = "gray95",notch = F, shape=21, outlier.size = -1, color="gray32",lwd=.5, alpha = .75)+
    theme(plot.title = element_text(size=8))+theme(plot.title = element_text(hjust = 0.5))+
    theme(axis.title.x = element_text(size=14))+
    theme(axis.text.x = element_text(colour = "black",size=14))+
    theme(axis.text.y = element_text(colour = "black",size=14))+
    theme( axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"))+
    theme(panel.grid.major.x = element_blank())+
    theme(panel.background = element_rect(fill = 'white'), panel.grid = element_line(color='gray80'))+
    # labs(title = addToTitle)+
    ylab("Row Means")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  return(p1)
}

## species_counts_table --###################################
# https://ibdmdb.org/tunnel/dataset_summary/HMP2/WGS/1818/summary/summary.html
species_counts_table <- fread("species_counts_table.tsv", header=T)
head(species_counts_table)
dim(species_counts_table)
summary(species_counts_table)
# it looks like we have 1301 individuals

# https://ibdmdb.org/tunnel/public/HMP2/WGS/1818/products
# merged tables
## ecs_3 --###################################
# ecs_3 <- fread("ecs_3.tsv.gz", header=T)
ecs_3 <- fread("https://ibdmdb.org/tunnel/products/HMP2/WGS/1818/ecs_3.tsv.gz", header=T)
ecs_3 <- as.data.frame(ecs_3)
dim(ecs_3)
str(ecs_3[1:5,1:4])
rownames(ecs_3) <- ecs_3$`Feature\\Sample`
ecs_3$`Feature\\Sample` <- NULL
# View(as.data.frame(rownames(ecs_3)))
ungrouped_ecs_3 <- ecs_3[grep("UNINTEGRATED",rownames(pathabundances_3), invert = F),]
grouped_ecs_3 <- ecs_3[grep("UNINTEGRATED",rownames(pathabundances_3), invert = T),]
num_ungrouped_ecs_3 <- mutate_all(ungrouped_ecs_3, function(x) as.numeric(as.character(x)))
num_grouped_ecs_3 <- mutate_all(grouped_ecs_3, function(x) as.numeric(as.character(x)))
# ROW MEANS
# rmeans <- as.data.frame(rowMeans(num_ungrouped_ecs_3)); colnames(rmeans) <- c("rmeans")
rmeans <- as.data.frame(rowMeans(num_grouped_ecs_3)); colnames(rmeans) <- c("rmeans")
summary(rmeans$rmeans)
c_make_hist(df = rmeans, whichCol = "rmeans")
c_make_violin(df = rmeans, whichCol = "rmeans")
# ROW STANDARD DEVIATIONS
rsds <- as.data.frame(apply(num_ungrouped_ecs_3, 1, sd)); colnames(rsds) <- c("rsds")
# rsds <- as.data.frame(apply(num_grouped_ecs_3, 1, sd)); colnames(rsds) <- c("rsds")
summary(rsds$rsds)
c_make_hist(df = rsds, whichCol = "rsds")
c_make_violin(df = rsds, whichCol = "rsds")

## pathabundances_3 --###################################
# pathabundances_3 <- fread("pathabundances_3.tsv.gz", header=T)
pathabundances_3 <- fread("https://ibdmdb.org/tunnel/products/HMP2/WGS/1818/pathabundances_3.tsv.gz", header=T)
pathabundances_3 <- as.data.frame(pathabundances_3)
dim(pathabundances_3)
str(pathabundances_3[1:5,1:4])
rownames(pathabundances_3) <- pathabundances_3$`Feature\\Sample`
pathabundances_3$`Feature\\Sample` <- NULL
# View(as.data.frame(rownames(pathabundances_3)))
ungrouped_pathabundances_3 <- pathabundances_3[grep("UNINTEGRATED",rownames(pathabundances_3), invert = F),]
grouped_pathabundances_3 <- pathabundances_3[grep("UNINTEGRATED",rownames(pathabundances_3), invert = T),]
num_ungrouped_pathabundances_3 <- mutate_all(ungrouped_pathabundances_3, function(x) as.numeric(as.character(x)))
num_grouped_pathabundances_3 <- mutate_all(grouped_pathabundances_3, function(x) as.numeric(as.character(x)))
# ROW MEANS
# rmeans <- as.data.frame(rowMeans(num_ungrouped_pathabundances_3)); colnames(rmeans) <- c("rmeans")
rmeans <- as.data.frame(rowMeans(num_grouped_pathabundances_3)); colnames(rmeans) <- c("rmeans")
summary(rmeans$rmeans)
c_make_hist(df = rmeans, whichCol = "rmeans")
c_make_violin(df = rmeans, whichCol = "rmeans")
# ROW STANDARD DEVIATIONS
rsds <- as.data.frame(apply(num_ungrouped_pathabundances_3, 1, sd)); colnames(rsds) <- c("rsds")
# rsds <- as.data.frame(apply(num_grouped_pathabundances_3, 1, sd)); colnames(rsds) <- c("rsds")
summary(rsds$rsds)
c_make_hist(df = rsds, whichCol = "rsds")
c_make_violin(df = rsds, whichCol = "rsds")

## taxonomic_profiles_3 --###################################
# taxonomic_profiles_3 <- fread("taxonomic_profiles_3.tsv.gz", header=T)
taxonomic_profiles_3 <- fread("https://ibdmdb.org/tunnel/products/HMP2/WGS/1818/taxonomic_profiles_3.tsv.gz", header=T)
taxonomic_profiles_3 <- as.data.frame(taxonomic_profiles_3)
dim(taxonomic_profiles_3)
str(taxonomic_profiles_3[1:5,1:4])
rownames(taxonomic_profiles_3) <- taxonomic_profiles_3$`Feature\\Sample`
taxonomic_profiles_3$`Feature\\Sample` <- NULL
# View(as.data.frame(rownames(taxonomic_profiles_3)))
ungrouped_taxonomic_profiles_3 <- taxonomic_profiles_3[grep("UNINTEGRATED",rownames(pathabundances_3), invert = F),]
grouped_taxonomic_profiles_3 <- taxonomic_profiles_3[grep("UNINTEGRATED",rownames(pathabundances_3), invert = T),]
num_ungrouped_taxonomic_profiles_3 <- mutate_all(ungrouped_taxonomic_profiles_3, function(x) as.numeric(as.character(x)))
num_grouped_taxonomic_profiles_3 <- mutate_all(grouped_taxonomic_profiles_3, function(x) as.numeric(as.character(x)))
# ROW MEANS
# rmeans <- as.data.frame(rowMeans(num_ungrouped_taxonomic_profiles_3)); colnames(rmeans) <- c("rmeans")
rmeans <- as.data.frame(rowMeans(num_grouped_taxonomic_profiles_3)); colnames(rmeans) <- c("rmeans")
summary(rmeans$rmeans)
c_make_hist(df = rmeans, whichCol = "rmeans")
c_make_violin(df = rmeans, whichCol = "rmeans")
# ROW STANDARD DEVIATIONS
rsds <- as.data.frame(apply(num_ungrouped_taxonomic_profiles_3, 1, sd)); colnames(rsds) <- c("rsds")
# rsds <- as.data.frame(apply(num_grouped_taxonomic_profiles_3, 1, sd)); colnames(rsds) <- c("rsds")
summary(rsds$rsds)
c_make_hist(df = rsds, whichCol = "rsds")
c_make_violin(df = rsds, whichCol = "rsds")

