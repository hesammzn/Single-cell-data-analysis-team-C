### Script Descriptions

---

### `prinseq.sh`

This script is for trimming the fastq files using prinseq with processing 6 files at a time which means 6 CPUs are needed. The perl file of prinseq is needed to run this script (version 0.20.4).

---

### `Trim_Galore.sh`

This script is for removing the nextera adapters from sequences. Trim Galore perl file (0.6.10) is needed to run this script.

---

### `index_star.sh`

This script is used for creating gene indices with STAR tool (version 2.7.11a). hg19 was used in this script and the annotation version is GRCh37.p13. This script should be submitted as a SLURM job.

---

### `STAR_MAP.sh`

This script is for mapping the reads to hg19 and creating gene count matrix at the end. This script should be submitted as a SLURM job.
