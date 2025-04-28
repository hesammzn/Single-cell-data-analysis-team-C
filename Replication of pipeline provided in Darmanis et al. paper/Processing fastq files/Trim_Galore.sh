#!/bin/bash

#This script is used for removing adapters from sequences. 


for acc in $(cat accessions.txt); do
    read1="${acc}_cleaned_1.fastq"
    read2="${acc}_cleaned_2.fastq"
    
    echo "Processing $read1 and $read2 with Trim Galore..."


    perl trim_galore --paired --stringency 1 --nextera "$read1" "$read2"

    echo "Finished processing $read1 and $read2"

done

