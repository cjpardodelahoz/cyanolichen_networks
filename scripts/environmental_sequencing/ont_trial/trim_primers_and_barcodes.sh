#!/bin/bash

#SBATCH --array=1-117%50
#SBATCH --mem-per-cpu=4G  # adjust as needed
#SBATCH -c 4 # number of threads per process
#SBATCH --output=log/environmental_sequencing/ont_trial/trim_primers_and_barcodes_%A_%a.out
#SBATCH --error=log/environmental_sequencing/ont_trial/trim_primers_and_barcodes_%A_%a.err
#SBATCH --partition=scavenger

# Load dependencies
source $(conda info --base)/etc/profile.d/conda.sh
conda activate cutadapt
# Sample variable
sample=$(cat misc_files/environmental_sequencing/ont_trial/all_trial_samples.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)
# Variables with paths
reads_path="analyses/environmental_sequencing/ont_trial/demultiplex/reads"
barcode_path="analyses/environmental_sequencing/ont_trial/demultiplex"
# Variables with forward and reverse primer-barcode sequences
forward=$(grep "${sample}" ${barcode_path}/barcoded_primers.tsv | cut -f 2)
reverse=$(grep "${sample}" ${barcode_path}/barcoded_primers.tsv | cut -f 3)
# Run cutadapt
cutadapt -g ${forward}...${reverse} \
 --revcomp \
 --discard-untrimmed \
 --minimum-length 1300 \
 --maximum-length 1600 \
 --cores 4 \
 -o ${reads_path}/${sample}.fastq \
 ${reads_path}/${sample}_barcoded.fastq