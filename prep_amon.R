setwd("/Users/johnsterrett/Research-Projects/poly-omics-scores/")
#setwd("~/OneDrive - UCB-O365/Documents/GRADSCHOOL/PolyOmicsRotation/poly-omics-risk2")
library(tidyverse)
library(data.table)

mgn <- fread("metagenomics/variable_importance.txt", sep="\t", drop=1)
# read ec table
mgn.func <- fread("metagenomics/ecs_3.tsv")
# filter weird row
mgn.func <- mgn.func[mgn.func$`Feature\\Sample`!="# Gene Family: NO_NAME",]
rownames(mgn.func) <- mgn.func$`Feature\\Sample`

list.of.ec.org <- mgn.func$`Feature\\Sample`

# get all ecs belonging to bugs that we selected for the models
bugs <- lapply(list.of.ec.org, 
                FUN=function(x) str_split(string=x, pattern=".s__")[[1]][2])

index.important <- which(bugs %in% mgn$variable)
ecs.org.important <- list.of.ec.org[index.important]

# filter the ec df
mgn.func.important <- mgn.func %>% as.data.frame() %>%
    filter( rownames(mgn.func) %in% ecs.org.important,
            !grepl( "UNMAPPED", x=rownames(mgn.func) ),
            !grepl( "UNGROUPED", x=rownames(mgn.func) )
            )
rownames(mgn.func.important) <- mgn.func.important$`Feature\\Sample`
mgn.func.important$`Feature\\Sample` <- NULL

mgn.func.important <- mgn.func.important %>% 
    mutate_all(as.numeric)

# filter out ECs that are very low abundance
mgn.func.important <- mgn.func.important[rowMeans(mgn.func.important) >= 0.001,]

# export the ecs
ecs <- rownames(mgn.func.important)
ecs.only <- lapply(ecs, 
                    FUN=function(x) str_split(x, pattern=":")[[1]][1])
unique.ecs <- unique(ecs.only)

fwrite(unique.ecs, file="selected_taxa_ECs.txt", sep="\t")


mbl <- fread("metabolomics/variable_importance.txt", sep="\t", drop=1)
metabolites <- mbl$variable

fwrite(as.list(metabolites), file="selected_metabolites.txt", sep="\t")
