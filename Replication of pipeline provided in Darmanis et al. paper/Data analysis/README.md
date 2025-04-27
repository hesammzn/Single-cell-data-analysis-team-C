### Script Execution Order

The scripts in this folder were executed in this order:

1. `Merging_gene_count_matrices.R`: This script is used to merge all the gene matrices produced with STAR (version 2.7.11a).

2. `Adding_gene_IDs.R`: Since merging the gene matrices do not have gene_ID column, we had to add it manually from one of the gene matrices.

3. `scde_script.R`: The script for calculating the pairwise distance, cell clustering, extracting top 20 expressed genes and creating a 3D plot with scde (version 2.27.1).

4. `Extracting_geneNames_and_geneIDs_from_annotation_file.sh`: This code was executed in terminal to get gene names with associated gene IDs from annotation GTF file.

5. `geneID_to_geneName.R`: Since we had gene IDs instead of gene names, we used this script to change the IDs to names.

---

### Other Files in the Folder

- `Modified_GeneCountMatrix.csv`: Gene counts across all cells with their IDs.

- `SLURM_SCDE.sh`: Bash script for running error model as SLURM since it takes so long.

- `SCDE_ErrorModels_100.rds`: The error model created with scde package.

- `top20_genes_names.csv`: Top 20 genes expressed in different clusters with the gene names instead of IDs.

- `gene_id_gene_name.tsv`: All gene IDs together with gene names obtained from hg19 GTF annotation file.
