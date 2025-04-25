library(Seurat)
library(dplyr)
library(plotly)
library(RColorBrewer)
library(ggplot2)

counts <- read.csv("C:/Users/hesam/Desktop/University of Leicester/Sem2/Steered Project/Alt/Alt_with_gene_names.csv", stringsAsFactors = FALSE)

#checking for number of duplicates
sum(duplicated(counts[,1]))

#Using dplyr library to sum the duplicate gene names
counts_agg <- counts %>%
  group_by(Gene = counts[[1]]) %>%
  summarise(across(where(is.numeric), sum))
counts_matrix <- as.data.frame(counts_agg)
rownames(counts_matrix) <- counts_matrix$Gene
counts_matrix$Gene <- NULL

#Checking again to see if we still have duplicates
sum(duplicated(rownames(counts_matrix)))


seurat_obj <- CreateSeuratObject(counts = counts_matrix, project = "MyProject", min.cells = 3, min.features = 200)

#Normalize
seurat_obj <- NormalizeData(seurat_obj)

#Variable features
seurat_obj <- FindVariableFeatures(seurat_obj)

#Scale
seurat_obj <- ScaleData(seurat_obj)

#PCA
seurat_obj <- RunPCA(seurat_obj)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes


#UMAP and Clustering
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 1.2)
DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 1)


#Identifying top expressed genes in each cluster
cluster_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#Extracting top 20 genes in each cluster
top20 <- cluster_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 20)

write.csv(top20, "top_markers_per_cluster_ALT.csv", row.names = FALSE)

#For checking our results
FeaturePlot(seurat_obj, features = c("VIP", "CCK", "CALB2", "GRIN3A"))

#Violin plot for cells in different cycles
VlnPlot(seurat_obj, features = c("S.Score", "G2M.Score"), group.by = "seurat_clusters", pt.size = 0.1)

table <- table(seurat_obj$seurat_clusters)
