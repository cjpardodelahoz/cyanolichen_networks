#!/bin/bash

#SBATCH --array=1-2
#SBATCH --mem-per-cpu=32G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/environmental_sequencing/ont_trial/demultiplex_raw_%A_%a.out
#SBATCH --error=log/environmental_sequencing/ont_trial/demultiplex_raw_%A_%a.err
#SBATCH --partition=scavenger

# Load dependencies
source $(conda info --base)/etc/profile.d/conda.sh
conda activate demultiplex
# Variable with pool names
pool=$(cat misc_files/environmental_sequencing/ont_trial/ont_pools.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)
# Demultiplex searching for forward and reverse primers separately
demultiplex match -m 2 -p analyses/environmental_sequencing/ont_trial/demultiplex/reads \
 analyses/environmental_sequencing/ont_trial/demultiplex/barcodes_${pool}_f_f.tsv \
 analyses/reads/ont_trial/${pool}_all.fastq
demultiplex match -m 2 -p analyses/environmental_sequencing/ont_trial/demultiplex/reads \
 analyses/environmental_sequencing/ont_trial/demultiplex/barcodes_${pool}_r_r.tsv \
 analyses/environmental_sequencing/ont_trial/demultiplex/reads/${pool}_all_UNKNOWN.fastq