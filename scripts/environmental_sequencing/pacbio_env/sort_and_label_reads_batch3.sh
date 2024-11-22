#!/bin/bash

#SBATCH --array=2-322
#SBATCH --mem-per-cpu=32G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/environmental_sequencing/pacbio_env/sort_and_label_reads_batch3.out
#SBATCH --error=log/environmental_sequencing/pacbio_env/sort_and_label_reads_batch3.err
#SBATCH --partition=scavenger

# Variable with the barcode combination
barcode=$(sed -n ${SLURM_ARRAY_TASK_ID}p misc_files/environmental_sequencing/pacbio_env/revio_order_10483_barcode.csv | cut -d, -f1)
# Variable witht he sample name
sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p misc_files/environmental_sequencing/pacbio_env/revio_order_10483_barcode.csv | cut -d, -f2)
# Make a directory for the copy of the reads
mkdir -p analyses/reads/pacbio_env
# Make a working copy of the reads with the sample name
cp data/reads/0000000406/outputs/fastx_files/*${barcode}*fastq.gz analyses/reads/pacbio_env/${sample}.fastq.gz
