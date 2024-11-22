#!/bin/bash

#SBATCH --mem-per-cpu=4G  # adjust as needed
#SBATCH -c 2 # number of threads per process
#SBATCH --output=log/lichen_sequencing/ont_lichen/collapse_multiscale_its.out
#SBATCH --error=log/lichen_sequencing/ont_lichen/collapse_multiscale_its.err
#SBATCH --partition=scavenger


# Script to collapse ITS and rbcLX sequences into haplotypes, including ONT sequences and sequences from Pardo-de la Hoz et al. 2024

# Output directories
mkdir -p analyses/lichen_sequencing/multiscale/haplotypes/its
mkdir -p analyses/lichen_sequencing/multiscale/haplotypes/rbclx

# Path variables
vsearch_path="/hpc/group/bio1/carlos/apps/vsearch-2.28.1-linux-x86_64/bin"
seqkit_path="/hpc/group/bio1/carlos/apps"
its_out_path="analyses/lichen_sequencing/multiscale/haplotypes/its"
rbclx_out_path="analyses/lichen_sequencing/multiscale/haplotypes/rbclx"

# Merge ONT and ABMI sequences
cat analyses/lichen_sequencing/ont_lichen/consensus/lichen_its_trimmed.fasta \
    data/sequences/its_abmi_pardodelahoz2024.fasta > \
    ${its_out_path}/its_multiscale.fasta
cat analyses/lichen_sequencing/ont_lichen/consensus/lichen_rbclx_trimmed.fasta \
    data/sequences/rbclx_abmi_pardodelahoz2024.fasta > \
    ${rbclx_out_path}/rbclx_multiscale.fasta  

# Remove rbclx sequences with length < 700 bp and its sequences with length < 500 bp
${seqkit_path}/seqkit seq -m 700 ${rbclx_out_path}/rbclx_multiscale.fasta > ${rbclx_out_path}/rbclx_multiscale_filtered.fasta
${seqkit_path}/seqkit seq -m 450 ${its_out_path}/its_multiscale.fasta > ${its_out_path}/its_multiscale_filtered.fasta 

# Collapse ITS sequences into haplotypes and print clustering report
${vsearch_path}/vsearch --cluster_fast ${its_out_path}/its_multiscale_filtered.fasta \
    --id 1 \
    --centroids ${its_out_path}/its_multiscale_haplotypes.fasta \
    --uc ${its_out_path}/its_multiscale_haplotypes.tsv \
    --qmask none \
    --threads 2

# Collapse rbcLX sequences into haplotypes and print clustering report
${vsearch_path}/vsearch --cluster_fast ${rbclx_out_path}/rbclx_multiscale_filtered.fasta \
    --id 1 \
    --centroids ${rbclx_out_path}/rbclx_multiscale_haplotypes.fasta \
    --uc ${rbclx_out_path}/rbclx_multiscale_haplotypes.tsv \
    --qmask none \
    --threads 2