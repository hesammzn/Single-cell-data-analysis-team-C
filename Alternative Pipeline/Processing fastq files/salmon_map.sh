#!/bin/bash

export PATH=$PATH:$(pwd)

while read acc; do
  echo "Quantifying $acc..."
  salmon quant -i salmon_index \
               -l A \
               -1 ${acc}_1_trimmed.fastq \
               -2 ${acc}_2_trimmed.fastq \
               -o ${acc}_quant \
               -p 4
  echo "Done with $acc"
done < accessions.txt
