# setwd("/Users/johnsterrett/Research-Projects/Team-rotation/poly-omics-scores/")
setwd("~/OneDrive - UCB-O365/Documents/GRADSCHOOL/PolyOmicsRotation/poly-omics-risk2")
library(tidyverse)
library(ggplot2)
library(data.table)
library(ggpubr)

mbl <- fread("metabolomics/variable_importance.txt", sep="\t", drop=1)
mgn <- fread("metagenomics/variable_importance.txt", sep="\t", drop=1)
mts <- fread("metatranscriptomics/variable_importance.txt", sep="\t", drop=1)
vir <- fread("viromics/variable_importance.txt", sep="\t", drop=1)

# Add a column for omic
mbl$omic <- "MBL"; mgn$omic <- "MGN"; mts$omic <- "MTS"; vir$omic <- "VIR"

# mbl add a column for column
mbl$Column <- mbl$Feature %>% 
    as.character() %>% 
    sapply(function(x) substr(x, 1, 4)) 

mbl$Column <- mbl$Column %>% 
    sapply(function(x) ifelse(x=="C18n", yes="C18 ",
                              no=ifelse(x=="C8p_", yes="C8",
                                        no=ifelse(x=="HILn", yes="HIL (neg)",
                                                  no="HIL (pos)"))))


mbl$variable <- mbl$variable %>% as.character()
mbl$new_variable <- mbl$variable


for(m in unique(mbl$variable)){
    if(nrow(mbl %>% filter(variable==m)) > 1){
        mbl[mbl$variable==m,] <- mbl[mbl$variable==m,] %>% 
            unite("variable", c(variable, Column), sep = " - ",
                  remove=F)
    }
}

# update MTS pathway names
mts$variable <- mts$variable %>% sapply(function(x) str_split(x, ":", simplify = T)[1,2])

# remove _
mgn$variable <- mgn$variable %>% sapply(function(x) gsub("_", " ", x, fixed=TRUE))
vir$variable <- vir$variable %>% sapply(function(x) gsub("_", " ", x, fixed=TRUE))


# plotting
make_plot <- function(data, feature_text_scaling) {
    
    data$variable <- reorder(data$variable, data$Estimate)
    
    ggplot(data) +
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
                size = feature_text_scaling,
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
}



mgn_plot <- make_plot(mgn, feature_text_scaling=8)
mts_plot <- make_plot(mts, feature_text_scaling=6)
vir_plot <- make_plot(vir, feature_text_scaling=8)
mbl_plot <- make_plot(mbl, feature_text_scaling=8)

p <- ggarrange(mgn_plot, mts_plot, vir_plot, mbl_plot,
               labels = c("MGN", "MTS", "VRM", "MBL"), align = "h", 
               label.x = c(0.5, 0.5, 0.5, 0.5) 
               )

p
ggsave("variable_importance_quad.png", height = 8)



