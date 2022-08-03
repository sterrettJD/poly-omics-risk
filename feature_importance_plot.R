setwd("/Users/johnsterrett/Research-Projects/Team-rotation/poly-omics-scores/")

mbl <- fread("metabolomics/variable_importance.txt", sep="\t", drop=1)
mgn <- fread("metagenomics/variable_importance.txt", sep="\t", drop=1)
mts <- fread("metatranscriptomics/variable_importance.txt", sep="\t", drop=1)
vir <- fread("viromics/variable_importance.txt", sep="\t", drop=1)

# Add a column for omic
mbl$omic <- "MBL"; mgn$omic <- "MGN"; mts$omic <- "MTS"; vir$omic <- "VIR"

# mbl add a column for column
mbl$Column <- mbl$Feature %>% 
    as.character() %>% 
    sapply(function(x) substr(x, 1, 3))





