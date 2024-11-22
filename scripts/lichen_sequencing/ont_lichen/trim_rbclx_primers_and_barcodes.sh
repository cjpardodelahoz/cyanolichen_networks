#!/bin/bash

#SBATCH --array=1-2638
#SBATCH --mem-per-cpu=4G  # adjust as needed
#SBATCH -c 4 # number of threads per process
#SBATCH --output=log/lichen_sequencing/ont_lichen/trim_rbclx_primers_and_barcodes_%A_%a.out
#SBATCH --error=log/lichen_sequencing/ont_lichen/trim_rbclx_primers_and_barcodes_%A_%a.err
#SBATCH --partition=scavenger

# Load dependencies
source $(conda info --base)/etc/profile.d/conda.sh
conda activate cutadapt
# Sample variable
sample=$(cat analyses/lichen_sequencing/ont_lichen/demultiplex/all_rbclx_barcodes_plus.tsv | cut -f 1 | sed -n ${SLURM_ARRAY_TASK_ID}p)
# Variables with paths
reads_path="analyses/lichen_sequencing/ont_lichen/demultiplex/rbclx"
barcode_path="analyses/lichen_sequencing/ont_lichen/demultiplex/"
# Variables with forward and reverse primer-barcode sequences
forward=$(grep "${sample}" ${barcode_path}/all_rbclx_barcodes_plus.tsv  | cut -f 2)
reverse=$(grep "${sample}" ${barcode_path}/all_rbclx_barcodes_plus.tsv  | cut -f 3)
# Trim rbclx reads
cutadapt -g ${forward}...${reverse} \
 --revcomp \
 --discard-untrimmed \
 --minimum-length 600 \
 --maximum-length 1000 \
 --cores 1 \
 -o ${reads_path}/${sample%%_plus}.fastq \
 ${reads_path}/pool*_all_${sample}.fastq