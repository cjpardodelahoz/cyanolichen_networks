#!/bin/bash

#SBATCH --array=2-385
#SBATCH --mem-per-cpu=16G  # adjust as needed
#SBATCH -c 4 # number of threads per process
#SBATCH --output=log/environmental_sequencing/pacbio_env/fastqc_batch2_%A_%a.out
#SBATCH --error=log/environmental_sequencing/pacbio_env/fastqc_batch2_%A_%a.err
#SBATCH --partition=scavenger

# Load FastQC module
module load FastQC/0.11.7
# # Variable witht he sample name
sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p misc_files/environmental_sequencing/pacbio_env/revio_order_10413_barcode.csv | cut -d, -f2)
# Run FastQC
fastqc --outdir analyses/environmental_sequencing/pacbio_env/qc/fastqc \
 -t 4 analyses/reads/pacbio_env/${sample}.fastq.gz