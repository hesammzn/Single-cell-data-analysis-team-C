The scripts in this folder were executed in this order:

1- Merging_gene_count_matrices:
	This script is used to merge all the gene matrices produced with STAR.

2- Adding_gene_IDs:
	Since merging the gene matrices do not have gene_ID column, we had to add it manually from one of the gene matrices.

3- scde_script:
	The script for calculating the pairwise distance, cell clustering, extracting top 20 expressed genes and creating 3D plots.

4- geneID_to_geneName:
	Since we had gene IDs instead of gene names, we used this script to change the IDs to names.

5- MHCI_expression:
	This scipt was used to create barplots to check MHCI pathway existence in the dataset.



Other files in the folder:

Modified_GeneCountMatrix:
	Gene counts across all cells with their IDs.

SLURM_SCDE:
	Bash script for running error model as SLURM since it takes so long.

SCDE_ErrorModels_100.rds:
	The error model created with scde package.

top20_genes_names:
	Top 20 genes expressed in different clusters with the gene names instead of IDs.




