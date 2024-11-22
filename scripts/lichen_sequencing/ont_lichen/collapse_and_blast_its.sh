#!/bin/bash

#SBATCH --mem-per-cpu=4G  # adjust as needed
#SBATCH -c 4 # number of threads per process
#SBATCH --output=log/lichen_sequencing/ont_lichen/collapse_and_blast_its.out
#SBATCH --error=log/lichen_sequencing/ont_lichen/collapse_and_blast_its.err
#SBATCH --partition=scavenger


# Script to identify ITS haplotypes from the lichen sequencing data and blast them against the ITS sequences from Pardo-de la Hoz et al. 2024

# Load conda dependencies
source $(conda info --base)/etc/profile.d/conda.sh
conda activate blast

# Create output directory
mkdir -p analyses/lichen_sequencing/ont_lichen/id/its

# Path variables
vsearch_path="/hpc/group/bio1/carlos/apps/vsearch-2.28.1-linux-x86_64/bin"

# Collapse sequences into haplotypes and print clustering report
${vsearch_path}/vsearch --cluster_fast analyses/lichen_sequencing/ont_lichen/consensus/lichen_its.fasta \
    --id 1 \
    --centroids analyses/lichen_sequencing/ont_lichen/id/its/its_haplotypes.fasta \
    --uc analyses/lichen_sequencing/ont_lichen/id/its/its_haplotypes.tsv \
    --qmask none \
    --threads 1

# Blast haplotypes against the ITS sequences from Pardo-de la Hoz et al. 2024 and print results as a table with qseqid, ssqid, pident, length, mismatch
blastn -query analyses/lichen_sequencing/ont_lichen/id/its/its_haplotypes.fasta \
    -db analyses/lichen_sequencing/ont_lichen/error_correction/its/its_abmi_pardodelahoz2024 \
    -outfmt "6 qseqid sseqid pident length mismatch" \
    -max_target_seqs 1 \
    -out analyses/lichen_sequencing/ont_lichen/id/its/its_haplotypes_blast.tsv \
    -num_threads 1