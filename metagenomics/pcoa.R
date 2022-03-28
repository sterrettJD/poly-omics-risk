library(tidyverse)
library(data.table)
library(vegan) #for distance matrix
library(ape) # for pcoa
`%ni%` <- Negate(`%in%`)

setwd("/Users/johnsterrett/Research-Projects/Team-rotation/poly-omics-scores/metagenomics")

metadata <- fread("https://ibdmdb.org/tunnel/products/HMP2/Metadata/hmp2_metadata.csv", header=T, stringsAsFactors=T)

taxprof <- fread("taxonomic_profiles_3.tsv.gz", header = T) %>% as.data.frame()

# set the rownames
rownames(taxprof) <- taxprof$`Feature\\Sample`
taxprof$`Feature\\Sample` <- NULL


# separate dataframe out by rows which were annotated vs rows which weren't
commentRow <- grep("# ", rownames(taxprof), invert = F)
mygrep <- grep("UNKNOWN", rownames(taxprof), invert = F)
grep_species <- grep("\\|s__", rownames(taxprof), invert = F)
mygrepni <- which(1:length(rownames(taxprof)) %ni% mygrep)
if(length(commentRow) > 0){
    mygrepni <- mygrepni[which(mygrepni != commentRow)]
}

# grab species level annotation
grouped_taxprof <- taxprof[grep_species,]
rownames(grouped_taxprof) <- rownames(taxprof)[grep_species]
num_grouped_taxprof <- mutate_all(grouped_taxprof, function(x) as.numeric(as.character(x)))

look <- as.data.frame(rownames(num_grouped_df_3)); colnames(look) <- c("look")
looksep <- look %>% separate(look,into=c("kingdom","phylum","class","order","family","genus","species"),convert=TRUE,sep="\\|")


# remove individuals with colsums that don't add up to 100
num_grouped_taxprof <- num_grouped_taxprof[,which(colSums(num_grouped_taxprof) > 5)]

# make a bray curtis distance matrix
metagenome_bc_distmat <- vegdist(t(num_grouped_taxprof), method="bray")

metagenome_pcoa <- pcoa(metagenome_bc_distmat)



