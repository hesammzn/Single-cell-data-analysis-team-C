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
- One based on **Darmanis et al.**
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

- **Reproduction of pipeline provided in Darmanis et al. paper**  
  Contains all of the scripts for the re-analysis of the pipeline provided in supplementary information.

- **Alternative pipeline**  
  Contains scripts for a pipeline produced with tools published after 2016.

- **Website**  
  Includes the code for creating the website in the LAMP server.


 ```├── bash_for_downloading_files.sh # SLURM script to download FASTQ files ├── accession/ # Accession files for downloading data ├── darmanis_pipeline/ # Scripts based on Darmanis et al. paper ├── alternative_pipeline/ # Custom pipeline with tools post-2016 ├── website/ # Code for building the website (LAMP server) └── README.md # Project overview ```
