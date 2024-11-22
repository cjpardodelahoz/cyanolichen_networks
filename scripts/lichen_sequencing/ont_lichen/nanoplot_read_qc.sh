#!/bin/bash

#SBATCH --array=1-8
#SBATCH --mem-per-cpu=16G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/lichen_sequencing/ont_lichen/read_qc_%A_%a.out
#SBATCH --error=log/lichen_sequencing/ont_lichen/read_qc_%A_%a.err
#SBATCH --partition=scavenger

# Load dependencies
source $(conda info --base)/etc/profile.d/conda.sh
conda activate nanoplot

# Directory for output
mkdir -p analyses/lichen_sequencing/ont_lichen/qc
# Get QC report for each pool
NanoPlot --fastq analyses/reads/ont_lichen/pool${SLURM_ARRAY_TASK_ID}_all.fastq.gz \
 --outdir analyses/lichen_sequencing/ont_lichen/qc \
 -p pool${SLURM_ARRAY_TASK_ID}_all_nanoplot
