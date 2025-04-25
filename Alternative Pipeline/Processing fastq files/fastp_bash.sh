#!/bin/bash


export PATH=$PATH:$(pwd)

while read acc; do
    echo "Processing $acc..."

    read1="${acc}_1.fastq"
    read2="${acc}_2.fastq"
    
    out1="${acc}_1_trimmed.fastq"
    out2="${acc}_2_trimmed.fastq"
    html_report="${acc}_report.html"

    # Run fastp with Nextera trimming and improved matching to prinseq
    fastp -i "$read1" -I "$read2" \
      -o "$out1" -O "$out2" \
      --adapter_sequence CTGTCTCTTATACACATCT \
      --adapter_sequence_r2 CTGTCTCTTATACACATCT \
      --trim_front1 10 --trim_front2 10 \
      --cut_tail --cut_mean_quality 25 \
      --length_required 30 \
      --thread 4 --html "$html_report"


    echo "Finished processing $acc"
    echo "------------------------"
done < accessions.txt
