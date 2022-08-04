library(tidyverse)
library(data.table)
library(vegan) #for distance matrix
library(ape) # for pcoa
library(ggplot2)
library(stringr)
`%ni%` <- Negate(`%in%`)

setwd("/Users/johnsterrett/Research-Projects/Team-rotation/poly-omics-scores/metagenomics")

metadata <- fread("https://ibdmdb.org/tunnel/products/HMP2/Metadata/hmp2_metadata.csv", header=T, stringsAsFactors=T)
metadata <- subset(metadata, data_type == "metagenomics")

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

look <- as.data.frame(rownames(num_grouped_taxprof)) 
colnames(look) <- c("look")
looksep <- look %>% separate(look,into=c("kingdom","phylum","class","order","family","genus","species"),convert=TRUE,sep="\\|")


# remove individuals with colsums that don't add up to 100
num_grouped_taxprof <- num_grouped_taxprof[,which(colSums(num_grouped_taxprof) > 5)]

# make a bray curtis distance matrix
metagenome_bc_distmat <- vegdist(t(num_grouped_taxprof), method="bray")

# Do PCoA
metagenome_pcoa <- pcoa(metagenome_bc_distmat)
metagenome_pcoa_coords <- as.data.frame(metagenome_pcoa$vectors)

#plot
ggplot(data = metagenome_pcoa_coords, mapping = aes(x=Axis.1, y=Axis.2, )) +
    geom_point()

# merge datasets
rownames(metagenome_pcoa_coords)
metadata$"External ID" 
rownames(metagenome_pcoa_coords) <- rownames(metagenome_pcoa_coords) %>% str_replace("_profile", "")
metagenome_pcoa_coords$`External ID` <- rownames(metagenome_pcoa_coords)
metagenome_pcoa_n_metadata <- merge(metagenome_pcoa_coords, metadata, by = "External ID")


ggplot(data = metagenome_pcoa_n_metadata, mapping = aes(x=Axis.1, y=Axis.6, color=diagnosis)) +
    geom_point()

# Scree plot
plot(metagenome_pcoa$values$Relative_eig[1:10])
# looks like the elbow is at 7

axis.cols <- metagenome_pcoa_n_metadata[,grep("Axis", 
                                              colnames(metagenome_pcoa_n_metadata), 
                                              fixed=T, invert=F)]

# Get cor with IBD  and then use the first 7
axis.cor<-cor(x=axis.cols,
              y=metagenome_pcoa_n_metadata$diagnosis!="nonIBD") %>%
    data.frame()
axis.cor.of.interest <- data.frame(axis.cor[1:7,])

# plot cor
ggplot(data = NULL, mapping = aes(x=rownames(axis.cor.of.interest), y=as.numeric(axis.cor.of.interest[,1]))) +
    geom_bar(stat = "identity") +
    xlab("Axis (of PCoA)") +
    ylab("Pearson correlation coef with IBD")

# Looks like axis 4 has the strongest cor with IBD
ggplot(data = metagenome_pcoa_n_metadata, mapping = aes(x=Axis.4, y=Axis.6, color=diagnosis)) +
    geom_point()

###### Is beta diversity significantly associated with IBD group?
new.dist.labels <- labels(metagenome_bc_distmat) %>% str_replace("_profile", "")
metagenome_bc_distmat<-usedist::dist_setNames(d=metagenome_bc_distmat, new.dist.labels)

diag <- as.factor(metadata$diagnosis)

diag.adonis.results <- adonis2(formula = metagenome_bc_distmat ~ diagnosis, 
                               data = metagenome_pcoa_n_metadata,  method=NULL,)


diag.adonis.results # SIGNIFICANT


# Dispoersion test
dispersion<-betadisper(metagenome_bc_distmat, group=metagenome_pcoa_n_metadata$diagnosis)
permutest(dispersion)

plot(dispersion)

