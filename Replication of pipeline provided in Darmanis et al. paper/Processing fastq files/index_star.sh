#!/bin/bash
#SBATCH --job-name=simple_job
#SBATCH --cpus-per-task=16       
#SBATCH --mem=64G                
#SBATCH --time=23:59:59          
#SBATCH --export=NONE            


#This script is for creating genome indexes for hg19 using STAR.

module load star


STAR --runThreadN 16 \
     --runMode genomeGenerate \
     --genomeDir STAR_Index \
     --genomeFastaFiles GRCh37.p13.genome.fa \
     --sjdbGTFfile gencode.v19.chr_patch_hapl_scaff.annotation.gtf \
     --sjdbOverhang 100 \
     --limitGenomeGenerateRAM 64000000000  

