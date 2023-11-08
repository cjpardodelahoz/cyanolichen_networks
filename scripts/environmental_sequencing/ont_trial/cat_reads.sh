#!/bin/bash

#SBATCH --mem-per-cpu=32G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/environmental_sequencing/ont_trial/cat_reads.out
#SBATCH --error=log/environmental_sequencing/ont_trial/cat_reads.err
#SBATCH --partition=scavenger

# Folder for copy of reads
mkdir -p analyses/reads/ont_trial
# Copy reads
cp data/reads/8729_Delivery/fastq_pass/barcode01_8729_Pool1/*.fastq.gz \
 analyses/reads/ont_trial
cp data/reads/8729_Delivery/fastq_pass/barcode02_8729_Pool2/*.fastq.gz \
 analyses/reads/ont_trial
# Unzip reads
gunzip analyses/reads/ont_trial/*fastq.gz
# Cat reads
cat analyses/reads/ont_trial/*barcode01* > analyses/reads/ont_trial/pool1_all.fastq
cat analyses/reads/ont_trial/*barcode02* > analyses/reads/ont_trial/pool2_all.fastq
# Remove individual read files
rm analyses/reads/ont_trial/*barcode*
