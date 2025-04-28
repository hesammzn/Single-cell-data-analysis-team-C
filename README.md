### Team C

- **Hesam Moazzen**  
  `hm435@student.le.ac.uk`

- **Abirajh Arulrajah**  
  `aa1152@student.le.ac.uk`

- **Sanskriti Sanskriti**  
  `s82@student.le.ac.uk`

- **Priyadarshini Venkataramesh**  
  `pv97@student.le.ac.uk`

- **Alice Haskell**  
  `aeh35@student.le.ac.uk`

---

### Project Overview

This repository contains scripts and tools used for single-cell RNA-seq data analysis using two pipelines:
- One based on **Darmanis et al. (2015)**
- Another using **tools published after 2016**

---

### `bash_for_downloading_files.sh`

This bash script will download the FASTQ files using the accession file (cell numbers) provided in the data section of the paper. It is a SLURM job which uses SRA toolkit version 3.0.0.

---

### `accession.txt`

Used for downloading the files.

---

### Repository Structure

There are 3 folders in this repository:

- **Replication of pipeline provided in Darmanis et al. paper**  
  Contains all of the scripts for the re-analysis of the pipeline provided in supplementary information.

- **Alternative pipeline**  
  Contains scripts for a pipeline produced with tools published after 2016.

- **Website**  
  Includes the code for creating the website in the LAMP server.

```
│   accession.txt
│   bash_for_downloading_files.sh
│   README.md
│
├───Alternative Pipeline
│   ├───Data analysis
│   │       Alt_with_gene_names.csv
│   │       Extracting_geneNames_and_geneIDs_from_annotation_file.sh
│   │       geneID_to_geneName.R
│   │       gene_id_gene_name.tsv
│   │       Merging_gene_counts.R
│   │       README.md
│   │       Seurat.R
│   │       top_markers_per_cluster_ALT.csv
│   │
│   └───Processing fastq files
│           accessions.txt
│           fastp_bash.sh
│           README.md
│           salmon_index.sh
│           salmon_map.sh
│
├───Replication of pipeline provided in Darmanis et al. paper
│   ├───Data analysis
│   │       Adding_gene_IDs.R
│   │       Extracting_geneNames_and_geneIDs_from_annotation_file.sh
│   │       geneID_to_geneName.R
│   │       gene_id_gene_name.tsv
│   │       Merging_gene_count_matrices.R
│   │       Modified_GeneCountMatrix.csv
│   │       README.md
│   │       SCDE_ErrorModels_100.rds
│   │       scde_script.R
│   │       SLURM_SCDE.sh
│   │       top20_genes_named.csv
│   │
│   └───Processing fastq files
│           index_star.sh
│           prinseq.sh
│           README.md
│           STAR_MAP.sh
│           trimgalore.sh
│
└───Website on the lamp server
        Website.R
 ```
