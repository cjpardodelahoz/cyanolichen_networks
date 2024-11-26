#!/bin/bash

#SBATCH --mem-per-cpu=16G  # adjust as needed
#SBATCH -c 4 # number of threads per process
#SBATCH --output=log/environmental_sequencing/pacbio_env/blast_nostoc_ref_batch123.out
#SBATCH --error=log/environmental_sequencing/pacbio_env/blast_nostoc_ref_batch123.err
#SBATCH --partition=scavenger

# Load dependencies
source $(conda info --base)/etc/profile.d/conda.sh
conda activate blast
# Outout directories
mkdir -p analyses/environmental_sequencing/pacbio_env/unoise/blast/dbs
mkdir -p analyses/environmental_sequencing/pacbio_env/unoise/blast/out
# Variables with paths and sample
sample="batch123_asvs"
out_path="analyses/environmental_sequencing/pacbio_env/unoise/blast/out"
db_path="analyses/environmental_sequencing/pacbio_env/unoise/blast/dbs"
seqs_path="analyses/environmental_sequencing/pacbio_env/unoise/seqs"

# Make BLAST database
makeblastdb -in data/placement_refs/16s_nostoc_ref.fasta \
    -dbtype nucl \
    -parse_seqids \
    -out ${db_path}/nostoc_ref_16s

# All-to-all BLASTn search
blastn -query ${seqs_path}/${sample}_nostocales.fasta \
    -db ${db_path}/nostoc_ref_16s \
    -outfmt '6 qseqid sseqid length pident nident mismatch gaps' \
    -out ${out_path}/${sample}_nostocales_to_nostoc_ref.txt

