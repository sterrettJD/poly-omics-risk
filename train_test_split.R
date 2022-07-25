#setwd("/Users/chris/Documents/GRADSCHOOL/PolyOmicsRotation")

setwd("/Users/johnsterrett/Research-Projects/Team-rotation/poly-omics-scores")
list.files()

library(data.table)
library(ggplot2)
library(R.utils)
library(tidyverse)
library(UpSetR)
library(cowplot)
#library(biomformat)
library(compositions)
library(bestglm)
library(MASS)
library(tree)
`%ni%` <- Negate(`%in%`)

## metadata --###################################
metadata <- fread("https://ibdmdb.org/tunnel/products/HMP2/Metadata/hmp2_metadata.csv", header=T, stringsAsFactors=T)
str(metadata)



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

subset_metadata <- metadata %>%
    subset(data_type %in% c("metabolomics", 
                            "metagenomics", 
                            "metatranscriptomics", 
                            "viromics"))


samples_with_all_omics <- metadata[which(metadata$`External ID` %in% idlist),]

samples_with_all_omics$ibd <- samples_with_all_omics$diagnosis!="nonIBD"

samples_with_all_omics %>% 
    count(ibd,`Participant ID`) %>% 
    group_by(ibd) %>%
    arrange(desc(n), by_group=ibd)

samples_with_all_omics[sample(1:nrow(samples_with_all_omics), 30),] %>% count(ibd,`Participant ID`)


dims <- c()


for(par in unique(metadata$`Participant ID`)){
    parRows <- nrow(metadata[metadata$`Participant ID`==par])
    print(nrow(metadata[metadata$`Participant ID`==par]))
    dims <- c(dims,parRows)
}
mean(dims)


grouped_meta <- group_by(metadata, col=`Participant ID`, col=data_type)

counts_df <- metadata %>% count(as.character(diagnosis)=="nonIBD", `Participant ID`, data_type) 

colnames(counts_df) <- c("nonIBD", "Participant ID", "data_type", "n")

sum(counts_df$n) == nrow(metadata)

counts_df <- subset(counts_df, 
                    data_type %in% c("metabolomics", 
                                     "metagenomics", 
                                     "metatranscriptomics", 
                                     "viromics"))

par_with_all_omics <- c()
par_without_all_omics <- c()
for(par in unique(counts_df$`Participant ID`)){
    participant_df <-subset(counts_df, `Participant ID` == par)
    nonibd <- participant_df$nonIBD[1]
    if(is.na(nonibd)){
        print(paste0("diagnosis empty for ", par))
    }else if(nonibd==FALSE){
        diagnosis <- "IBD"
    }else{
        diagnosis <- "nonIBD"
    }
    
    
    if(nrow(participant_df)!=4){
        print("LESS THAN 4 omics")
        nsamples <- sum(participant_df$n)
        nomics <- sum(participant_df$n>0) 
        print(paste0(diagnosis, " sample ", par, " has ", nsamples, " samples over ", nomics, " omics layers" ))
        par_without_all_omics <- c(par_without_all_omics, par)
    }else{
        nsamples <- sum(participant_df$n)
        nomics <- sum(participant_df$n>0) 
        print(paste0(diagnosis, " sample ", par, " has ", nsamples, " samples over ", nomics, " omics layers" ))
        par_with_all_omics <- c(par_with_all_omics, par)
        
        
    }
}

length(par_with_all_omics)
length(par_without_all_omics)
length(unique(counts_df$`Participant ID`)) == (length(par_with_all_omics) +length(par_without_all_omics))

# only 1 participant doesn't have any of the 4
# 26 don't have 4 omics
# 104 have all 4 omics

counts_df[counts_df$`Participant ID` %in% length(par_without_all_omics)]

# sample 30 participants with all omics layers
sampled <- par_with_all_omics[sample(1:length(par_with_all_omics), 30)]


# grab those from the dataframe
sampled_df <- metadata[metadata$`Participant ID` %in% sampled]
sampled_df %<>% 
    count(as.character(diagnosis)=="nonIBD", `Participant ID`, data_type) %>%
    subset(data_type %in% c("metabolomics", 
                            "metagenomics", 
                            "metatranscriptomics", 
                            "viromics"))

summary(sampled_df)

length(unique(subset_metadata$`External ID`))

total_samples_per_participant <- c()
for(par in unique(sampled_df$`Participant ID`)){
    participant_df <-subset(metadata, `Participant ID` == par)
    print(paste0(par, " has ", nrow(participant_df), " total samples"))
    total_samples_per_participant <- c(total_samples_per_participant, nrow(participant_df))
}

summary(total_samples_per_participant)

dim(metadata[metadata$`Participant ID` %in% sampled])
dim(metadata)


parlist <- c()
quadn <- c()
nonquadn <- c()
for(par in unique(subset_metadata$`Participant ID`)){
    pars_df <- subset(subset_metadata, `Participant ID` == par)
    pars_df_quad <- subset(pars_df, `External ID` %in% idlist)
    pars_df_nonquad <- subset(pars_df, `External ID` %ni% idlist)
    parsquad_nsamples <- nrow(pars_df_quad)
    parsnonquad_nsamples <- nrow(pars_df_nonquad)
    
    parlist <- c(parlist, par)
    quadn <- c(quadn, parsquad_nsamples)
    nonquadn <- c(nonquadn, parsnonquad_nsamples)
}

pars_quad_df <- data.frame(parlist, quadn, nonquadn)
pars_quad_df$sum <- pars_quad_df$quadn + pars_quad_df$nonquadn 
pars_quad_df$percent <- pars_quad_df$quadn/pars_quad_df$sum 

pars_quad_df[order(-pars_quad_df$percent, -pars_quad_df$sum),]

pars_quad_df[order(pars_quad_df$nonquadn),] %>% 
    subset(quadn>mean(quadn)) %>% arrange(desc(by=percent))

plotvecquad <- c()
plotvecnonquad <- c()
for(i in 1:50){
    top_par_quad <- pars_quad_df[order(pars_quad_df$nonquadn),] %>% 
        subset(quadn>mean(quadn)) %>% 
        arrange(desc(by=percent)) 
    
    top_par_quad <- top_par_quad[1:i,]
    #top_par_quad <- top_par_quad[sample(1:nrow(top_par_quad), 26),]
    
    colSums(top_par_quad[,c("quadn", "nonquadn")])
    
    plotvecquad <- c(plotvecquad, sum(top_par_quad[,c("quadn")]))
    plotvecnonquad <- c(plotvecnonquad, sum(top_par_quad[,c("nonquadn")]))
}

plotdf <- data.frame(1:50,plotvecquad,plotvecnonquad)
colnames(plotdf) <- c("index", "quad", "nonquad")
ggplot(data = plotdf, aes(x = index)) + 
    geom_point(aes(y = quad), color = "red") + 
    geom_point(aes(y = nonquad), color = "blue")



top_par_quad <- pars_quad_df[order(pars_quad_df$nonquadn),] %>% 
    subset(quadn>mean(quadn)) %>% 
    arrange(desc(by=percent)) 

top_par_quad <- top_par_quad[1:30,]
#top_par_quad <- top_par_quad[sample(1:nrow(top_par_quad), 26),]

sum(colSums(top_par_quad[,c("quadn", "nonquadn")]))

top_par_quad$uniquequadsamplen <- top_par_quad$quadn/4
top_par_quad


testing_metadata <- subset(subset_metadata, `Participant ID` %in% top_par_quad$parlist)

testing_metadata <- testing_metadata[,c("External ID",
                                        "Participant ID",
                                        "race",
                                        "data_type",
                                        "consent_age",
                                        "site_name",
                                        "diagnosis",
                                        "sex",
                                        "Antibiotics")]



training_metadata <- subset(subset_metadata, `Participant ID` %ni% top_par_quad$parlist)

nrow(training_metadata) + nrow(testing_metadata) == nrow(subset_metadata)



training_metadata <- training_metadata[,c("External ID",
                                          "Participant ID",
                                          "race",
                                          "data_type",
                                          "consent_age",
                                          "site_name",
                                          "diagnosis",
                                          "sex",
                                          "Antibiotics")]

summary(training_metadata$diagnosis)
summary(testing_metadata$diagnosis)

dim(training_metadata[!duplicated(training_metadata$`Participant ID`),])
dim(testing_metadata[!duplicated(testing_metadata$`Participant ID`),])

summary(training_metadata[!duplicated(training_metadata$`Participant ID`),diagnosis])
summary(testing_metadata[!duplicated(testing_metadata$`Participant ID`),diagnosis])



# ONLY RUN THIS IF YOU WANT TO DEAL WITH THE CONSEQUENCES
testing_metadata %>%
    write.table("name.txt", sep = "\t", row.names=F, col.names=F, quote=F)
training_metadata %>%
    write.table("name.txt", sep = "\t", row.names=F, col.names=F, quote=F)