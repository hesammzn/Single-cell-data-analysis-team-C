library(tximport)
library(readr)
library(GenomicFeatures)


#Listing all files ending in .sf
files <- list.files("/scratch/alice/h/hm435/Alternative/trimmed_with_fastp_adapter/all_quants", pattern = "_quant\\.sf$", full.names = TRUE)

#Extracting sample names
names(files) <- gsub("_quant\\.sf$", "", basename(files))
head(files)

txdb <- makeTxDbFromGFF("/scratch/alice/h/hm435/Alternative/trimmed_with_fastp_adapter/all_quants/gencode.v47.annotation.gtf", format = "gtf")

k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
head(tx2gene)

#Cleaning the transcript names by removing the tx version manually
tx2gene$TXNAME <- sub("\\..*", "", tx2gene$TXNAME)

txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, ignoreAfterBar = TRUE)
dim(txi$counts)
head(txi$counts[, 1:3])
write.csv(txi$counts, "gene_counts_matrix.csv")