library(ggplot2)

counts <- read.csv("C:/Users/hesam/Desktop/University of Leicester/Sem2/Steered Project/Old/100/Gene_Count_matrix_all_genes_100_with_names.csv",
                   stringsAsFactors = FALSE, check.names = FALSE)

counts$gene_name <- gsub('"', '', counts$gene_name)
counts$gene_name <- gsub("'", "", counts$gene_name)

genes_of_interest <- c("HLA-A", "HLA-B", "HLA-C", "B2M", "TAPBP", "CALR", "ERAP1", "HSPA5", 
                       "PDIA3", "TAP2", "SEC61A1", "SEC61A2", "SEC61B", "SEC61G")

filtered_counts <- counts[counts$gene_name %in% genes_of_interest, ]

library(dplyr)
filtered_agg <- filtered_counts %>%
  group_by(gene_name) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE))

filtered_matrix <- as.data.frame(filtered_agg)
rownames(filtered_matrix) <- filtered_matrix$gene_name
filtered_matrix$gene_name <- NULL

gene_means <- rowMeans(filtered_matrix)

barplot(gene_means,
        las = 2,
        col = "steelblue",
        main = "Average Expression of MHCI Genes",
        ylab = "Mean Expression (Counts)")
