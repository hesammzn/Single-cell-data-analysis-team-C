gene_map <- read.table("/scratch/alice/h/hm435/Alternative/trimmed_with_fastp_adapter/all_quants/gene_id_gene_name.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(gene_map) <- c("gene_id", "gene_name")

expr_table <- read.csv("/home/h/hm435/gene_counts_matrix.csv", stringsAsFactors = FALSE)

colnames(expr_table)[1] <- "Gene"
expr_table$gene_short <- sub("\\..*", "", expr_table$Gene)

gene_map$gene_short <- sub("\\..*", "", gene_map$gene_id)

merged <- merge(expr_table, gene_map, by = "gene_short", all.x = TRUE)

merged$Gene <- merged$gene_name
merged <- subset(merged, select = -c(gene_short, gene_id, gene_name))

write.csv(merged, "Alt_with_gene_names.csv", row.names = FALSE)


