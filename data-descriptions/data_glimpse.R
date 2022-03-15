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
`%ni%` <- Negate(`%in%`)

pdf("metagenomics_data_glimpse.pdf", width = 8.5, height = 11)

quick_hist = function(values_vec, breaks=50) {
  res = hist(values_vec, plot=FALSE, breaks=breaks)
  
  dat = data.frame(xmin=head(res$breaks, -1L),
                   xmax=tail(res$breaks, -1L),
                   ymin=0.0,
                   ymax=res$counts)
  
  ggplot(dat, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)) +
    geom_rect(size=0.5, colour="grey30", fill="grey80")
}

c_make_hist <- function(df, whichCol){
  p1 <- ggplot(df, aes_string(x=whichCol)) + 
    geom_histogram(color="black", fill="salmon", alpha = 0.5, bins = 30) +
    geom_vline(aes(xintercept=mean(.data[[whichCol]])),
               color="blue", linetype="dashed", size=1) +
    theme(axis.title.x=element_blank())
  return(p1)
}

c_make_violin <- function(df, whichCol){
  p1 <- ggplot(data = df,aes_string(x = NA, y = whichCol))+
    scale_fill_viridis_d( option = "D")+
    theme_dark()+
    geom_violin(fill = "gray70",alpha=0.4, position = position_dodge(width = .5),size=1,color="gray22",width=.5,lwd=.2) +
    geom_boxplot(fill = "gray95",notch = F, shape=21, outlier.size = -1, color="gray32",lwd=.5, alpha = .75)+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(axis.title.x = element_text(size=14))+
    theme(axis.text.x = element_text(colour = "black"))+
    theme(axis.text.y = element_text(colour = "black"))+
    theme( axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"))+
    theme(panel.grid.major.x = element_blank())+
    theme(panel.background = element_rect(fill = 'white'), panel.grid = element_line(color='gray80'))+
    # labs(title = addToTitle)+
    # ylab("Row Means")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  return(p1)
}




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

# subset to jsut one datatype to explore diagnosis counts
metadata <- subset(metadata, data_type == "metagenomics")
metadata <- subset(metadata, `External ID` %in% idlist)

# fill missinig enries with NA and view thte columns that are mostly non-NA
metadata[metadata==""] <- NA
isna <- sapply(metadata, function(x) sum(is.na(x)))
isna[isna < 100]
summary(metadata$diagnosis)



## species_counts_table --###################################
# https://ibdmdb.org/tunnel/dataset_summary/HMP2/WGS/1818/summary/summary.html
species_counts_table <- fread("species_counts_table.tsv", header=T)
head(species_counts_table)
dim(species_counts_table)
summary(species_counts_table)
# it looks like we have 1301 individuals


data_glimpse <- function(fname, UNINTEGRATED = F, thresh = NA, normalize = T, myTitle=""){
  ecs_3 <- fname
  print("data dimensions:")
  print(dim(ecs_3))
  print("data glimpse:")
  str(ecs_3[1:5,1:4])
  rownames(ecs_3) <- ecs_3$`Feature\\Sample`
  ecs_3$`Feature\\Sample` <- NULL
  # View(as.data.frame(rownames(ecs_3)))
  # separate dataframe out by rows which were annotatted vs rows which weren't
  commentRow <- grep("# ", rownames(ecs_3), invert = F)
  mygrep <- grep("UNINTEGRATED|UNGROUPED|UNMAPPED", rownames(ecs_3), invert = F)
  mygrepni <- which(1:length(rownames(ecs_3)) %ni% mygrep)
  if(length(commentRow) > 0){
    mygrepni <- mygrepni[which(mygrepni != commentRow)]
  }
  if(UNINTEGRATED == T){
    ungrouped_ecs_3 <- ecs_3[mygrep,]
    rownames(ungrouped_ecs_3) <- rownames(ecs_3)[mygrep]
    num_ungrouped_ecs_3 <- mutate_all(ungrouped_ecs_3, function(x) as.numeric(as.character(x)))
  }else if(UNINTEGRATED == F){
    grouped_ecs_3 <- ecs_3[mygrepni,]
    rownames(grouped_ecs_3) <- rownames(ecs_3)[mygrepni]
    num_grouped_ecs_3 <- mutate_all(grouped_ecs_3, function(x) as.numeric(as.character(x)))
  }
  # normalize by sample sum
  if(normalize == T){
    if(UNINTEGRATED == T){
      num_ungrouped_ecs_3 <- as.data.frame(scale(num_ungrouped_ecs_3, center=FALSE, scale=colSums(num_ungrouped_ecs_3)))
    }else if(UNINTEGRATED == F){
      num_grouped_ecs_3 <- as.data.frame(scale(num_grouped_ecs_3, center=FALSE, scale=colSums(num_grouped_ecs_3)))
    }
  }
  # ROW MEANS
  if(UNINTEGRATED == T){
    rmeans <- as.data.frame(rowMeans(num_ungrouped_ecs_3, na.rm=T)); colnames(rmeans) <- c("rmeans")
  }else if(UNINTEGRATED == F){
    rmeans <- as.data.frame(rowMeans(num_grouped_ecs_3, na.rm=T)); colnames(rmeans) <- c("rmeans")
  }
  q1 <- quantile(rmeans, 0.8, na.rm=T)
  rmeans <- c(rmeans[rmeans < q1])
  # summary(rmeans$rmeans)
  # p1 <- c_make_hist(df = rmeans, whichCol = "rmeans") + ggtitle("Row Means")
  p1 <- quick_hist(rmeans) + theme_bw() + ggtitle("Row Means")
  p2 <- c_make_violin(df = as.data.frame(rmeans), whichCol = "rmeans") + ggtitle("Row Means") + ylab("Row Means")
  # ROW STANDARD DEVIATIONS
  if(UNINTEGRATED == T){
    rsds <- as.data.frame(apply(num_ungrouped_ecs_3, 1, sd, na.rm=T)); colnames(rsds) <- c("rsds")
  }else if(UNINTEGRATED == F){
    rsds <- as.data.frame(apply(num_grouped_ecs_3, 1, sd, na.rm=T)); colnames(rsds) <- c("rsds")
  }
  q2 <- quantile(rsds, 0.8, na.rm=T)
  rsds <- c(rsds[rsds < q2])
  # summary(rsds$rsds)
  # p3 <- c_make_hist(df = rsds, whichCol = "rsds") + ggtitle("Row SDs")
  p3 <- quick_hist(rsds) + theme_bw() + ggtitle("Row SDs")
  p4 <- c_make_violin(df = as.data.frame(rsds), whichCol = "rsds") + ggtitle("Row SDs") + ylab("Row SDs")
  # ALL ENTRIES
  if(UNINTEGRATED == T){
    alldf <- as.data.frame(as.numeric(array(as.matrix(num_ungrouped_ecs_3)))); colnames(alldf) <- c("alldf")
  }else if(UNINTEGRATED == F){
    alldf <- as.data.frame(as.numeric(array(as.matrix(num_grouped_ecs_3)))); colnames(alldf) <- c("alldf")
  }
  # summary(alldf$alldf) #this is super useful for choosing the thresh
  if(is.na(thresh)){
    myquantthresh <- quantile(alldf$alldf[alldf$alldf>0 & alldf$alldf < 0.0005], 0.85, na.rm=T)
  }else{
    myquantthresh <- quantile(alldf$alldf[alldf$alldf>0 & alldf$alldf < 0.0005], thresh, na.rm=T)
  }
  
  p5 <- quick_hist(alldf$alldf[alldf$alldf>0 & alldf$alldf < myquantthresh]) + theme_bw() + ggtitle(paste0("All Non-Zero Entries Below ", myquantthresh))
  pg1 <- plot_grid(p2, p1, p4, p3,  labels = c('A', 'B', 'C', 'D'), label_size = 12, ncol = 2)
  pg2 <- plot_grid(p5, labels = c('E'), label_size = 12, ncol = 1)
  pg3 <- plot_grid(pg1, pg2, ncol = 1)
  title <- ggdraw() + draw_label(myTitle, fontface='bold')
  pg <- plot_grid(title, pg3, ncol=1, rel_heights=c(0.1, 1)) # rel_heights values control title margins
  return(pg)
}



# Metagenomes
# https://ibdmdb.org/tunnel/public/HMP2/WGS/1818/products
# merged tables
## Metagenomes ecs_3 --###################################
fname <- fread("https://ibdmdb.org/tunnel/products/HMP2/WGS/1818/ecs_3.tsv.gz", header=T)
# fname <- fread("ecs_3.tsv.gz", header=T)
fname <- as.data.frame(fname)
data_glimpse(fname = fname, UNINTEGRATED = F, thresh = 0.8, normalize = T, myTitle = "Metagenomes ecs_3")

## Metagenomes pathabundances_3 --###################################
fname <- fread("https://ibdmdb.org/tunnel/products/HMP2/WGS/1818/pathabundances_3.tsv.gz", header=T)
# fname <- fread("pathabundances_3.tsv.gz", header=T)
fname <- as.data.frame(fname)
data_glimpse(fname = fname, UNINTEGRATED = F, thresh = 0.8, normalize = T, myTitle = "Metagenomes pathabundances_3")

## Metagenomes taxonomic_profiles_3 --###################################
fname <- fread("https://ibdmdb.org/tunnel/products/HMP2/WGS/1818/taxonomic_profiles_3.tsv.gz", header=T)
# fname <- fread("taxonomic_profiles_3.tsv.gz", header=T)
fname <- as.data.frame(fname)
data_glimpse(fname = fname, UNINTEGRATED = F, thresh = 0.8, normalize = T, myTitle = "Metagenomes taxonomic_profiles_3")



# Metatranscriptomes
# https://ibdmdb.org/tunnel/public/HMP2/MTX/1750/products
# merged tables
## Metatranscriptomes ecs_3 --###################################
fname <- fread("https://ibdmdb.org/tunnel/products/HMP2/MTX/1750/ecs_3.tsv.gz", header=T)
# fname <- fread("ecs_3.tsv.gz", header=T)
fname <- as.data.frame(fname)
data_glimpse(fname = fname, UNINTEGRATED = F, thresh = 0.8, normalize = T, myTitle = "Metatranscriptomes ecs_3")

## Metatranscriptomes pathabundances_3 --###################################
fname <- fread("https://ibdmdb.org/tunnel/products/HMP2/MTX/1750/pathabundances_3.tsv.gz", header=T)
# fname <- fread("pathabundances_3.tsv.gz", header=T)
fname <- as.data.frame(fname)
data_glimpse(fname = fname, UNINTEGRATED = F, thresh = 0.8, normalize = T, myTitle = "Metatranscriptomes pathabundances_3")





# # https://ibdmdb.org/tunnel/public/HMP2/Metabolites/1723/products
# # merged tables
# ## metabolomics ecs_3 --###################################
# fname <- fread("https://ibdmdb.org/tunnel/products/HMP2/Metabolites/1723/HMP2_metabolomics.csv.gz", header=T)
# # fname <- fread("ecs_3.tsv.gz", header=T)
# fname <- as.data.frame(fname)
# data_glimpse(fname = fname, UNINTEGRATED = F, thresh = 0.8, normalize = T, myTitle = "metabolomics ecs_3")


# # https://ibdmdb.org/tunnel/public/HMP2/Viromics/1732/products
# # merged tables
# ## viromics ecs_3 --###################################
# # fname <- fread("https://ibdmdb.org/tunnel/products/HMP2/Viromics/1732/taxonomic_profiles.biom.gz", header=T)
# fname <- read_biom("taxonomic_profiles.biom.gz")
# # fname <- fread("ecs_3.tsv.gz", header=T)
# fname <- as.data.frame(fname)
# data_glimpse(fname = fname, UNINTEGRATED = F, thresh = 0.8, normalize = T, myTitle = "viromics ecs_3")
# 
# ## viromics ecs_3 --###################################
# fname <- fread("https://ibdmdb.org/tunnel/products/HMP2/Viromics/1732/virome_virmap_analysis.tsv.gz", header=T)
# # fname <- fread("ecs_3.tsv.gz", header=T)
# fname <- as.data.frame(fname)
# data_glimpse(fname = fname, UNINTEGRATED = F, thresh = 0.8, normalize = T, myTitle = "viromics ecs_3")


dev.off()
