#!/bin/bash
#SBATCH --job-name=simple_job
#SBATCH --cpus-per-task=16       
#SBATCH --mem=64G                
#SBATCH --time=23:59:59          
#SBATCH --export=NONE            


#This script is for mapping the reads to hg19 and creating a gene matrix.

module load star

for acc in $(cat accessions.txt); do
    # Define input file names
    read1="${acc}_cleaned_1_val_1.fq"
    read2="${acc}_cleaned_2_val_2.fq"
    
    
    echo "Processing $read1 and $read2 with STAR..."
    
    
    
    STAR --readFilesIn "$read1" "$read2" \
    	 --runThreadN 16 \
   	 --genomeDir /lustre/alice3/scratch/alice/h/hm435/fastq_files/FASTQ2/cleaned/cleaned_with_trimgalore/STAR/STAR_Index \
    	 --outFilterType BySJout \
    	 --outFileNamePrefix "${acc}_" \
    	 --outFilterMultimapNmax 20 \
    	 --alignSJoverhangMin 8 \
    	 --alignSJDBoverhangMin 1 \
    	 --outFilterMismatchNmax 999 \
    	 --outFilterMismatchNoverLmax 0.04 \
    	 --alignIntronMin 20 \
    	 --alignIntronMax 1000000 \
    	 --alignMatesGapMax 1000000 \
   	 --outSAMstrandField intronMotif \
   	 --quantMode GeneCounts

done

