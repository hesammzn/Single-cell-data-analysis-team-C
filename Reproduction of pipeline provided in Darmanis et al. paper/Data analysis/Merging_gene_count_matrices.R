path <- "C:/Users/hesam/Desktop/University of Leicester/Sem2/Steered Project/Matrices"

files <- list.files(path = path, pattern = "ReadsPerGene.out.tab", full.names = TRUE)

count_list <- list()


for (file in files) {
  df <- read.table(file, header = FALSE, skip = 4)
  df <- df[, 2, drop = FALSE]
  
  colnames(df) <- gsub("_ReadsPerGene.out.tab", "", basename(file))
  
  df$Gene_ID <- rownames(df)
  rownames(df) <- NULL
  
  count_list[[basename(file)]] <- df
}

merged_counts <- Reduce(function(x, y) full_join(x, y, by = "Gene_ID"), count_list)

write.csv(merged_counts, file.path(path, "GeneCountMatrix.csv"), row.names = TRUE)


