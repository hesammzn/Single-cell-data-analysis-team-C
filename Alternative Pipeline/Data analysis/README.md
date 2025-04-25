The files in this folder were executed in this order:

1- Merging_gene_counts.R:
	We used this script to import transcript-level quantification results and aggregate them into gene-level count data using tximport (version 1.28.0).

2- Extracting_geneNames_and_geneIDs_from_annotation_file.sh:
	This code was executed in terminal to get gene names with associated gene IDs from annotation GTF file.

3- geneID_to_geneName.R:
	We changed the gene IDs to gene names using this script.

4- Seurat.R:
	Script for dimension reduction, clustering and different plots with Seurat (version 5.2.1).

5- MHC_gene_expression.R:
	Script for analysing the existence of MHC pathway in the dataset.
	



Other files in this folder:

gene_id_gene_name.tsv:
	Tsv file for gene IDs together with gene names extracted from GTF annotation file.

Alt_with_gene_names.csv:
	Table for all genes expressed in the analysis with their name.

top_markers_per_cluster_ALT.csv:
	20 most expressed genes among each cluster.