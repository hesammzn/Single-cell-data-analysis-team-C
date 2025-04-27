library(Seurat)
library(dplyr)
library(plotly)
library(RColorBrewer)
library(ggplot2)
library(scatterplot3d)

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
rownames(counts_matrix) <- toupper(rownames(counts_matrix))

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

write.csv(top20, "top_markers_per_cluster_ALT222222222.csv", row.names = FALSE)

#For checking our results
FeaturePlot(seurat_obj, features = c("AQP4", "TMEM144", "AIF1", "TNR", "SATB2", "TAC1", "CCK", "B2M", "TOP2A"))

seurat_obj <- CellCycleScoring(
  object = seurat_obj,
  s.features = s.genes,
  g2m.features = g2m.genes,
  set.ident = TRUE
)
#Violin plot for cells in different cycles
VlnPlot(seurat_obj, features = c("S.Score", "G2M.Score"), group.by = "seurat_clusters", pt.size = 0.1)

table <- table(seurat_obj$seurat_clusters)



cluster_ids <- c(
  "0" = "Fetal_Quiescent\nnewly born neurons",
  "1" = "Neuroendocrines",
  "2" = "Astrocytes",
  "3" = "GABAergic neurons",
  "4" = "Oligodendrocytes",
  "5" = "OPCs",
  "6" = "Microglia",
  "7" = "Endothelial cells",
  "8" = "Fetal_Replicating neuronal\nprogenitors"
)


seurat_obj$celltype <- plyr::mapvalues(
  x = seurat_obj$seurat_clusters,
  from = names(cluster_ids),
  to = cluster_ids
)
celltype_counts <- table(seurat_obj$celltype)

new_labels <- paste0(names(celltype_counts), " (n=", celltype_counts, ")")
names(new_labels) <- names(celltype_counts)

celltype_with_counts_vector <- new_labels[seurat_obj$celltype]

names(celltype_with_counts_vector) <- colnames(seurat_obj)

seurat_obj <- AddMetaData(seurat_obj, metadata = celltype_with_counts_vector, col.name = "celltype_with_counts")


#Creating figure with clusters
DimPlot(seurat_obj, reduction = "umap", group.by = "celltype_with_counts", label = TRUE, repel = TRUE, pt.size = 1, label.size = 3) +
  ggtitle("Clusters with number of cells") +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 8),
    plot.caption = element_text(size = 8)
  ) +
  theme(legend.position = "right") +
  labs(x = "UMAP_1", y = "UMAP_2")



#3D plot
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10, n.components = 3)

umap_coords <- Embeddings(seurat_obj, "umap") %>% as.data.frame()
umap_coords$celltype_with_counts <- seurat_obj$celltype_with_counts

celltypes <- unique(umap_coords$celltype_with_counts)
colors <- RColorBrewer::brewer.pal(n = length(celltypes), name = "Set3")
color_map <- setNames(colors, celltypes)
point_colors <- color_map[as.character(umap_coords$celltype_with_counts)]

point_colors <- as.character(point_colors)  
scatterplot3d(
  x = umap_coords$umap_1,
  y = umap_coords$umap_2,
  z = umap_coords$umap_3,
  color = point_colors,
  pch = 16,
  angle = 55,
  scale.y = 0.7,
  main = "3D UMAP Clustering with Cell Type Labels",
  xlab = "UMAP 1",
  ylab = "UMAP 2",
  zlab = "UMAP 3"
)

#Analyzing MHC pathway type I in different cells.
mhci_genes <- c("HLA-A", "HLA-B", "HLA-C", "B2M", "TAPBP", "CALR", "ERAP1", "HSPA5", 
                "PDIA3", "TAP2", "SEC61A1", "SEC61A2", "SEC61B", "SEC61G")
mhci_genes_found <- mhci_genes[mhci_genes %in% rownames(seurat_obj)]
mhci_data <- GetAssayData(seurat_obj, slot = "data")[mhci_genes_found, ]

celltypes <- seurat_obj$celltype

mhci_means_by_celltype <- as.data.frame(t(mhci_data)) %>%
  mutate(CellType = celltypes) %>%
  group_by(CellType) %>%
  summarise(across(all_of(mhci_genes_found), mean))

mhci_melted <- reshape2::melt(mhci_means_by_celltype, id.vars = "CellType", variable.name = "Gene", value.name = "AverageExpression")

#Plotting the expression pattern
ggplot(mhci_melted, aes(x = Gene, y = AverageExpression, fill = CellType)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  labs(title = "MHCI Gene Expression by Cell Type", y = "Average Expression", x = "Gene") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(size = 6)
  ) + scale_fill_brewer(palette = "Set3")



cluster_ids <- c(
  "0" = "Fetal_Quiescent\nnewly born neurons",
  "1" = "Adult cells",
  "2" = "Adult cells",
  "3" = "Adult cells",
  "4" = "Adult cells",
  "5" = "Adult cells",
  "6" = "Adult cells",
  "7" = "Adult cells",
  "8" = "Fetal_Replicating neuronal\nprogenitors"
)


seurat_obj$celltype <- plyr::mapvalues(
  x = seurat_obj$seurat_clusters,
  from = names(cluster_ids),
  to = cluster_ids
)
celltype_counts <- table(seurat_obj$celltype)

new_labels <- paste0(names(celltype_counts), " (n=", celltype_counts, ")")
names(new_labels) <- names(celltype_counts)

celltype_with_counts_vector <- new_labels[seurat_obj$celltype]

names(celltype_with_counts_vector) <- colnames(seurat_obj)

seurat_obj <- AddMetaData(seurat_obj, metadata = celltype_with_counts_vector, col.name = "celltype_with_counts")


#Creating figure with clusters for adult and fetal cells
DimPlot(seurat_obj, reduction = "umap", group.by = "celltype_with_counts", label = TRUE, repel = TRUE, pt.size = 1, label.size = 3) +
  ggtitle("Clusters with number of cells") +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 8),
    plot.caption = element_text(size = 8)
  ) +
  theme(legend.position = "right") +
  labs(x = "UMAP_1", y = "UMAP_2")


#Validating manual annotation with PanglaoDB
glia_markers <- read.delim("C:/Users/hesam/Desktop/University of Leicester/PanglaoDB_markers_27_Mar_2020.tsv")


colnames(glia_markers)


glia_markers$official.gene.symbol <- toupper(glia_markers$official.gene.symbol)

astrocyte_markers <- glia_markers$official.gene.symbol[glia_markers$cell.type == "Astrocytes"]
oligodendrocyte_markers <- glia_markers$official.gene.symbol[glia_markers$cell.type == "Oligodendrocytes"]
opc_markers <- glia_markers$official.gene.symbol[glia_markers$cell.type == "Oligodendrocyte progenitor cells"]
microglia_markers <- glia_markers$official.gene.symbol[glia_markers$cell.type == "Microglia"]
neuron_markers <- glia_markers$official.gene.symbol[glia_markers$cell.type == "Neurons"]
vascular_markers <- glia_markers$official.gene.symbol[glia_markers$cell.type == "Endothelial cells"]
Neuroendocrine_cells <- glia_markers$official.gene.symbol[glia_markers$cell.type == "Neuroendocrine cells"]
GABAergic_neurons <- glia_markers$official.gene.symbol[glia_markers$cell.type == "GABAergic neurons"]
Radial_glia_cells <- glia_markers$official.gene.symbol[glia_markers$cell.type == "Radial glia cells"]


seurat_obj <- AddModuleScore(seurat_obj, features = list(astrocyte_markers), name = "Astrocytes")
seurat_obj <- AddModuleScore(seurat_obj, features = list(oligodendrocyte_markers), name = "Oligodendrocytes")
seurat_obj <- AddModuleScore(seurat_obj, features = list(opc_markers), name = "Oligodendrocyte progenitor cells")
seurat_obj <- AddModuleScore(seurat_obj, features = list(microglia_markers), name = "Microglia")
seurat_obj <- AddModuleScore(seurat_obj, features = list(neuron_markers), name = "Neuron")
seurat_obj <- AddModuleScore(seurat_obj, features = list(vascular_markers), name = "Vascular")
seurat_obj <- AddModuleScore(seurat_obj, features = list(Neuroendocrine_cells), name = "Neuroendocrine cells")
seurat_obj <- AddModuleScore(seurat_obj, features = list(GABAergic_neurons), name = "GABAergic neurons")
seurat_obj <- AddModuleScore(seurat_obj, features = list(Radial_glia_cells), name = "Radial glia cells")


VlnPlot(seurat_obj, 
        features = c("Astrocytes1", "Oligodendrocytes1", "Oligodendrocyte progenitor cells1", "Microglia1", "Neuron1", "Vascular1", "Neuroendocrine cells1","GABAergic neurons1","Radial glia cells1"), 
        group.by = "seurat_clusters", 
        pt.size = 0.1, 
        ncol = 3) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 10)
  )
