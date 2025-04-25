#!/bin/bash

#Accession hm435

#SBATCH --job-name=simple_job
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=23:59:59
#SBATCH --export=NONE


#This script is for downloading files using SRA toolkit.

module load sratoolkit/3.0.0-5fetwpi

while read -r accession; do
    echo "Downloading $accession..."
    fasterq-dump "$accession" --split-files --skip-technical
    if [[ $? -eq 0 ]]; then
        echo "$accession download completed."
    else
        echo "Failed to download $accession."
    fi
done < accession.txt

echo "All downloads completed."
