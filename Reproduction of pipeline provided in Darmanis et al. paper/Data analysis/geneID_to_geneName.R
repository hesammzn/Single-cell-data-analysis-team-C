gene_map <- read.table("C:/Users/hesam/Desktop/University of Leicester/Sem2/Steered Project/CleanCOunts/gene_id_gene_name.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(gene_map) <- c("gene_id", "gene_name")

expr_table <- read.csv("C:/Users/hesam/Desktop/University of Leicester/Sem2/Steered Project/CleanCOunts/Top20_Genes_Per_Cluster_CleanCounts.csv", stringsAsFactors = FALSE)

expr_table$gene_short <- sub("\\..*", "", expr_table$Gene)
gene_map$gene_short <- sub("\\..*", "", gene_map$gene_id)

merged <- merge(expr_table, gene_map, by = "gene_short", all.x = TRUE)

merged$Gene <- merged$gene_name
final_table <- merged[, c("cluster", "Gene", "AvgExpression")]

write.csv(final_table, "top20_genes_named_CleanCounts.csv", row.names = FALSE)
