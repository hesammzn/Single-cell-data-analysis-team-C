library(readr)

csv_file <- "C:/Users/hesam/Desktop/University of Leicester/Sem2/Steered Project/Matrices/GeneCountMatrix.csv"
tsv_file <- "C:/Users/hesam/Desktop/University of Leicester/Sem2/Steered Project/Matrices/SRR1974543_ReadsPerGene.out.tab"

csv_data <- read.csv(csv_file, stringsAsFactors = FALSE)
tsv_data <- read.table(tsv_file, sep = "\t", header = FALSE, skip = 4)

new_first_column <- tsv_data[, 1]
csv_data[, 1] <- new_first_column

write.csv(csv_data, "modified_gene_count_matrix.csv", row.names = FALSE)
