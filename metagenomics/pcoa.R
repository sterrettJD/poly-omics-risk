library(tidyverse)
library(data.table)
setwd("/Users/johnsterrett/Research-Projects/Team-rotation/poly-omics-scores/metagenomics")

metadata <- fread("https://ibdmdb.org/tunnel/products/HMP2/Metadata/hmp2_metadata.csv", header=T, stringsAsFactors=T)

taxprof <- fread("taxonomic_profiles_3.tsv.gz", header = T) %>% as.data.frame()

# set the rownames
rownames(taxprof) <- taxprof$`Feature\\Sample`
taxprof$`Feature\\Sample` <- NULL


