setwd("/Users/johnsterrett/Research-Projects/Team-rotation/poly-omics-scores/")

library(tidyverse)
library(ggplot2)

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
mbl$variable <- mbl$variable %>% as.character()
mbl$new_variable <- mbl$variable


for(m in unique(mbl$variable)){
    if(nrow(mbl %>% filter(variable==m)) > 1){
        #mbl[mbl$variable==m, "variable"] <- paste(mbl[mbl$variable==m, "variable"],
        #                                           mbl[mbl$variable==m, "Column"]) 
        
        mbl[mbl$variable==m,] <- mbl[mbl$variable==m,] %>% 
            unite("variable", c(variable, Column),
                  remove=F)
    }
}





mbl <- mbl %>% dplyr::arrange(Estimate)

ggplot(mbl) +
    geom_segment(
        aes(
            x = variable,
            y = 0,
            xend = variable,
            yend = Estimate
        ),
        size = 1.5,
        alpha = 0.7
    ) +
    geom_point(aes(x = variable, y = Estimate, col = variable),
               size = 2,
               show.legend = F) +
    coord_flip() +
    labs(y = "Weight", x = NULL, title = "") +
    theme_bw() +
    theme(legend.title = element_text(size = 14)) +
    theme(
        axis.text.x = element_text(
            color = "black",
            size = 8,
            angle = 0,
            hjust = .5,
            vjust = .5
        ),
        axis.text.y = element_text(
            color = "black",
            size = length(unique(mbl$variable))*1/5,
            angle = 0
        ),
        axis.title.x = element_text(
            color = "black",
            size = 8,
            angle = 0
        ),
        axis.title.y = element_text(
            color = "black",
            size = 13,
            angle = 90
        )
    )



