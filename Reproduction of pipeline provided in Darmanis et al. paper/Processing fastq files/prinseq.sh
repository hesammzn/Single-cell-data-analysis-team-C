#!/bin/bash


module load fastqc

#This script is for trimming data with prinseq and checking the quality of them using fastqc.

process_sample() {
    acc="$1"
    echo "Processing $acc..."
    
    
    read1="${acc}_1.fastq"
    read2="${acc}_2.fastq"

    
    out_prefix="${acc}_cleaned"
    discarded_prefix="${acc}_discarded"
    
    
    perl prinseq-lite.pl -fastq "$read1" -fastq2 "$read2" \
        -out_format 3 \
        -out_good "$out_prefix" -out_bad "$discarded_prefix" \
        -min_len 30 -trim_left 10 -trim_qual_right 25 \
        -lc_method entropy -lc_threshold 65

    echo "Finished processing $acc"
}

export -f process_sample  


cat accessions.txt | parallel -j 6 process_sample

