#!/bin/sh
#SBATCH --job-name=hm435
#SBATCH --time=23:59:59
#SBATCH --cpus-per-task=1
#SBATCH --mem=32g
#SBATCH --export=NONE


module load R/4.3.1

cd $SLURM_SUBMIT_DIR

R --vanilla -s -f /home/h/hm435/Desktop/Steered_project/R/scde.R