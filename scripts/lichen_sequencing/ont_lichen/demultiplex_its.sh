#!/bin/bash

#SBATCH --array=1-8
#SBATCH --mem-per-cpu=32G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/lichen_sequencing/ont_lichen//demultiplex_its_%A_%a.out
#SBATCH --error=log/lichen_sequencing/ont_lichen/demultiplex_its_%A_%a.err
#SBATCH --partition=scavenger

# Load dependencies
source $(conda info --base)/etc/profile.d/conda.sh
conda activate demultiplex

# Make a copy of the reads to unzip
cp analyses/reads/ont_lichen/pool${SLURM_ARRAY_TASK_ID}_all.fastq.gz \
    analyses/lichen_sequencing/ont_lichen/demultiplex/its/pool${SLURM_ARRAY_TASK_ID}_all.fastq.gz
gunzip analyses/lichen_sequencing/ont_lichen/demultiplex/its/pool${SLURM_ARRAY_TASK_ID}_all.fastq.gz

# Demultiplex ITS reads with barcodes in plus strand
demultiplex match -m 3 -p analyses/lichen_sequencing/ont_lichen/demultiplex/its \
    analyses/lichen_sequencing/ont_lichen/demultiplex/pool${SLURM_ARRAY_TASK_ID}_its_barcoded_plus.tsv \
    analyses/lichen_sequencing/ont_lichen/demultiplex/its/pool${SLURM_ARRAY_TASK_ID}_all.fastq

# Remove the copy of the reads
rm analyses/lichen_sequencing/ont_lichen/demultiplex/its/pool${SLURM_ARRAY_TASK_ID}_all.fastq