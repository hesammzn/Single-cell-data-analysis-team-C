library(Seurat)
library(ggplot2)
library(Matrix)
library(dplyr)
library(reshape2)


counts <- read.csv("Alt_with_gene_names.csv", stringsAsFactors = FALSE, check.names = FALSE)

# Sum duplicate gene names
counts_agg <- counts %>%
  group_by(Gene = counts[[1]]) %>%
  summarise(across(where(is.numeric), sum))

# Convert to matrix with gene names as rownames
counts_matrix <- as.data.frame(counts_agg)
rownames(counts_matrix) <- counts_matrix$Gene
counts_matrix$Gene <- NULL

# Define genes of interest
genes_of_interest <- c("HLA-A", "HLA-B", "HLA-C", "B2M", "TAPBP", "CALR", "ERAP1", "HSPA5", 
                       "PDIA3", "TAP2", "SEC61A1", "SEC61A2", "SEC61B", "SEC61G")

# Filter for genes of interest
filtered_counts <- counts_matrix[rownames(counts_matrix) %in% genes_of_interest, ]

# Save result
write.csv(filtered_counts, "MHCI_expression_per_cell.csv")


gene_means <- rowMeans(filtered_counts)

barplot(gene_means,
        las = 2,
        col = "steelblue",
        main = "Average Expression of MHCI Genes",
        ylab = "Mean Expression (Counts)")
