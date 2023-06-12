library(tidyverse)
library(data.table)
library(VennDiagram)
library(plotrix)

all.v.human <- fread("amon_out_all_v_human/origin_table.tsv") %>% 
    filter(detected==T)

prod.microbe <- sum(all.v.human$all_taxa)
prod.host <- sum(all.v.human$human)
prod.both <- sum(all.v.human$all_taxa & all.v.human$human)

# move to new plotting page
grid.newpage()
draw.pairwise.venn(area1=prod.microbe, area2=prod.host, cross.area=prod.both, 
                   cex = 3, cat.cex = 2, cat.dist = -0.08,
                   category=c("Produced by microbes",
                              "Produced by host"),
                   fill=c("Yellow","Blue"))

all.v.human[all.v.human$all_taxa & !all.v.human$human,]


selected.v.human <- fread("amon_out_selected_v_human/origin_table.tsv") %>% 
    filter(detected==T)

prod.microbe <- sum(selected.v.human$selected_taxa)
prod.host <- sum(selected.v.human$human)
prod.both <- sum(selected.v.human$selected_taxa & selected.v.human$human)

# move to new plotting page
grid.newpage()
draw.pairwise.venn(area1=prod.microbe, area2=prod.host, cross.area=prod.both, 
                   cex = 3, cat.cex = 2, cat.dist = -0.08,
                   category=c("Produced by selected microbes",
                              "Produced by host"),
                   fill=c("Yellow","Blue"))


selected.v.all <- fread("amon_out_selected_v_all/origin_table.tsv") %>% 
    filter(detected==T)

prod.all <- sum(selected.v.all$all_taxa)
prod.selected <- sum(selected.v.all$selected_taxa)
prod.both <- sum(selected.v.all$all_taxa & selected.v.all$selected_taxa)

# move to new plotting page
grid.newpage()
draw.pairwise.venn(area1=prod.selected, area2=prod.all, cross.area=prod.both, 
                   cex = 3, cat.cex = 2, cat.dist = -0.08,
                   category=c("Produced by selected microbes",
                              "Produced by any microbes"),
                   fill=c("Yellow","Blue"))

selected.v.all

