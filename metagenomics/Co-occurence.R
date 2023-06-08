#Set working direction
setwd("~/Documents/PolyOmic_IBD")

library(dplyr)

metagenomics <- read.csv("MGN_preLASSO_training_DF.csv")
df2 <- read.csv("MGN_postLASSO_training_DF.csv")
meta_wod <- metagenomics[, -c(238, 239, 240, 241, 242, 243, 244, 245, 246, 247)]
meta_wod
features_wod <- df2[, -c(15, 16, 17, 18, 19, 20, 21, 22, 23, 24)]
features_wod

#Save the txt files just in case.
write.table(meta_wod, file = "meta_wod.txt", sep = "\t", row.names = FALSE)
write.table(meta_wod, file = "features_wod.txt", sep = "\t", row.names = FALSE)

#build the pearson correlation index matrix 
cor_matrix <- cor(meta_wod, method = "pearson")
cor_matrix_Df2 <- cor(features_wod, method = "pearson")

#save the matrixs just in case.
write.table(cor_matrix, file = "cor_matrix.txt", sep = "\t", row.names = FALSE)
write.table(cor_matrix_Df2, file = "cor_matrix_Df2.txt", sep = "\t", row.names = FALSE)

# Set heatmap colors, choose any of them. I choose the second one.
# heatmap_colors <- colorRampPalette(c("black", "red"))(2)
heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)

# Create the heatmap using pheatmap and save it.

#For total metagenome
p <- pheatmap(cor_matrix, 
         color = heatmap_colors,
         main = "Co-occurrence heatmap",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         scale = "none",
         fontsize_row = 3, # set row label font size
         fontsize_col = 3) # set column label font size
ggsave("heatmap_all.png", plot = p, dpi = 300, width = 8, height = 6, units = "in")

#For selected features
j <-pheatmap(cor_matrix_Df2, 
         color = heatmap_colors,
         main = "Co-occurrence heatmap",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         scale = "none")
ggsave("heatmap_14features.png", plot = j, dpi = 300, width = 8, height = 6, units = "in")

#Let's see how many pairs have >0.5 Pearson correlation index, because it is difficult just observing the complex heatmap for the total metagenome

#How many >0.5
high_cor_indices <- which(cor_matrix > 0.5 & upper.tri(cor_matrix), arr.ind = TRUE)
num_high_cor_pairs <- nrow(high_cor_indices)
cat("There are", num_high_cor_pairs, "pairs with correlation coefficient > 0.5.\n")

#Pairs names:
high_cor_indices <- which(cor_matrix_Df2 > 0.5 & upper.tri(cor_matrix), arr.ind = TRUE)
# Get the variable names for each index
high_cor_pairs <- rownames(cor_matrix)[high_cor_indices[,1]] + " - " + rownames(cor_matrix)[high_cor_indices[,2]]
cat("The pairs with correlation coefficient > 0.5 are:\n", high_cor_pairs, "\n")
