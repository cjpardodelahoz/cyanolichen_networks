#!/bin/bash

#SBATCH --mem-per-cpu=32G  # adjust as needed
#SBATCH -c 16 # number of threads per process
#SBATCH --output=log/environmental_sequencing/pacbio_env/unoise_asv_call.out
#SBATCH --error=log/environmental_sequencing/pacbio_env/unoise_asv_call.err
#SBATCH --partition=scavenger

# Unoise output directory
mkdir -p analyses/environmental_sequencing/pacbio_env/unoise/seqs
mkdir -p analyses/environmental_sequencing/pacbio_env/unoise/counts
# Path variables
reads_path="analyses/reads/pacbio_env"
vsearch_path="/hpc/group/bio1/carlos/apps/vsearch-2.28.1-linux-x86_64/bin"
seqs_path="analyses/environmental_sequencing/pacbio_env/unoise/seqs"
counts_path="analyses/environmental_sequencing/pacbio_env/unoise/counts"

# ASV calling pipeline

# Pool all the reads from batch 1
cat ${reads_path}/*_trimmed.fastq.gz > ${reads_path}/batch123_pooled.fastq.gz
# Remove sequences < 1300 bp and > 1600 bp and with expected errors > 1
${vsearch_path}/vsearch --fastq_filter ${reads_path}/batch123_pooled.fastq.gz \
 --fastaout ${seqs_path}/batch123_filtered.fasta \
 --fastq_maxee 1 \
 --fastq_minlen 1300 \
 --fastq_maxlen 1600
# Dereplicate the reads and remove seqs with abundance < 4
${vsearch_path}/vsearch --derep_fulllength ${seqs_path}/batch123_filtered.fasta \
 --output ${seqs_path}/batch123_derep.fasta \
 --sizeout \
 --minuniquesize 10 \
 --relabel asv
# Call ASVs
${vsearch_path}/vsearch --cluster_unoise ${seqs_path}/batch123_derep.fasta \
 --id 1 \
 --centroids ${seqs_path}/batch123_asvs.fasta \
 --unoise_alpha 2.0 \
 --minsize 8 \
 --threads 16
# Filter chimeras
${vsearch_path}/vsearch --uchime_denovo ${seqs_path}/batch123_asvs.fasta \
 --nonchimeras ${seqs_path}/batch123_asvs_nonchimeras.fasta
# Map reads to ASVs
${vsearch_path}/vsearch --usearch_global ${reads_path}/batch123_pooled.fastq.gz \
 --db ${seqs_path}/batch123_asvs_nonchimeras.fasta \
 --id 0.99 \
 --threads 16 \
 --otutabout ${counts_path}/batch123_otutab.txt
