library(scde)
library(Matrix)
library(parallel)
library(tsne)
library(mclust)
library(scatterplot3d)
library(FactoMineR)
library(dplyr)

counts <- read.csv("/home/h/hm435/Desktop/Steered_project/R/Modified_GeneCountMatrix.csv", 
                   row.names = 1, header = TRUE)
n.cores <- 1
counts <- as.matrix(counts)
counts <- counts[rowSums(counts) >= 100, ]

scde.fitted.model <- readRDS("/home/h/hm435/Desktop/Steered_project/R/second_error_with100/SCDE_ErrorModels.rds")

scde.prior <- scde.expression.prior(models = scde.fitted.model, counts = counts)

jp <- scde.posteriors(models = scde.fitted.model, counts, scde.prior, 
                      return.individual.posterior.modes = TRUE, n.cores = n.cores)

jp$jp.modes <- log(as.numeric(colnames(jp$jp)))[max.col(jp$jp)]
p.mode.fail <- scde.failure.probability(models = scde.fitted.model, magnitudes = jp$jp.modes)
p.self.fail <- scde.failure.probability(models = scde.fitted.model, counts = counts)

# Weight and magnitude matrices
matw <- 1 - sqrt(p.self.fail * sqrt(p.self.fail * p.mode.fail))
mat <- log10(exp(jp$modes) + 1)

#PCA with FactoMineR on transposed matrix
expr_data <- t(mat) 
expr_data <- as.data.frame(expr_data)
pca_res <- PCA(expr_data, scale.unit = TRUE, ncp = 30, graph = FALSE)
pca_matrix <- pca_res$ind$coord

#t-SNE on top PCs
tsne_res <- tsne(pca_matrix, k = 3)

#GMM clustering
gmm_res <- Mclust(tsne_res, G = 10)
clusters <- gmm_res$classification
colors <- rainbow(length(unique(clusters)))[clusters]

#3D plot
scatterplot3d(tsne_res[,1], tsne_res[,2], tsne_res[,3],
              color = colors,
              pch = 20,
              xlab = "t-SNE 1", ylab = "t-SNE 2", zlab = "t-SNE 3",
              xlim = c(-60, 60), ylim = c(-60, 60), zlim = c(-60, 60),
              main = "Clusters obtained with pipeline in paper")

# Convert to vector
clusters <- as.matrix(clusters)
cluster_vec <- as.vector(clusters)
names(cluster_vec) <- colnames(counts)

# Get top 20 genes per cluster
top_genes_by_cluster <- list()

for (cl in sort(unique(cluster_vec))) {
  cl_cells <- names(cluster_vec)[cluster_vec == cl]
  cl_counts <- counts[, cl_cells, drop = FALSE]
  gene_means <- rowMeans(cl_counts)
  top_genes <- sort(gene_means, decreasing = TRUE)[1:20]
  top_genes_by_cluster[[paste0("Cluster_", cl)]] <- top_genes
}

#Formatting top genes into a data frame
top_genes_df <- do.call(rbind, lapply(names(top_genes_by_cluster), function(cluster_name) {
  data.frame(
    cluster = cluster_name,
    Gene = names(top_genes_by_cluster[[cluster_name]]),
    AvgExpression = as.numeric(top_genes_by_cluster[[cluster_name]])
  )
}))


write.csv(top_genes_df, "Top20_Genes_Per_Cluster_100.csv", row.names = FALSE)
