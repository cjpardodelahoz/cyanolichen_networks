#!/bin/bash

#SBATCH --array=1-8
#SBATCH --mem-per-cpu=16G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/lichen_sequencing/ont_lichen/cat_reads_%A_%a.out
#SBATCH --error=log/lichen_sequencing/ont_lichen/cat_reads_%A_%a.err
#SBATCH --partition=scavenger


# Merge reads by pool and make working copy
cat data/reads/10298_CP_Pool1/10298_CPP1/20240828_1258_X1_FAZ71248_ff4b8dff/fastq_pass/barcode0${SLURM_ARRAY_TASK_ID}/*.fastq.gz > \
  analyses/reads/ont_lichen/pool${SLURM_ARRAY_TASK_ID}_all.fastq.gz