#!/bin/bash

#SBATCH --array=1-2
#SBATCH --mem-per-cpu=16G  # adjust as needed
#SBATCH -c 4 # number of threads per process
#SBATCH --output=log/environmental_sequencing/ont_trial/fastqc_%A_%a.out
#SBATCH --error=log/environmental_sequencing/ont_trial/fastqc_%A_%a.err
#SBATCH --partition=scavenger

# Load fastqc
module load FastQC/0.11.7
# Pool variable
pool=$(cat misc_files/environmental_sequencing/ont_trial/ont_pools.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)
# Run Fastqc
fastqc --outdir analyses/environmental_sequencing/ont_trial/qc -t 4 analyses/reads/ont_trial/${pool}_all.fastq