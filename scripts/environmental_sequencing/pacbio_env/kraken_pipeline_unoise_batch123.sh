#!/bin/bash

#SBATCH --mem-per-cpu=8G  # adjust as needed
#SBATCH -c 4 # number of threads per process
#SBATCH --output=log/environmental_sequencing/pacbio_env/kraken_pipeline_unoise_batch123.out
#SBATCH --error=log/environmental_sequencing/pacbio_env/kraken_pipeline_unoise_batch123.err
#SBATCH --partition=scavenger

# Load dependencies
source $(conda info --base)/etc/profile.d/conda.sh
conda activate kraken2
# Outout directory
mkdir -p analyses/environmental_sequencing/pacbio_env/unoise/taxonomy/kraken
# Variables with paths and sample
sample="batch123_asvs"
out_path="analyses/environmental_sequencing/pacbio_env/unoise/taxonomy/kraken"
db_path="/hpc/group/bio1/carlos/apps/kraken2/greengenes/greengenes"
seqs_path="analyses/environmental_sequencing/pacbio_env/unoise/seqs"

# Assign taxonomy to ASVs

# Run kraken2
kraken2 --d ${db_path} \
 --threads 4 \
 ${seqs_path}/${sample}_nonchimeras.fasta \
 --report ${out_path}/${sample}.k2report > \
 ${out_path}/${sample}.kraken2

# Extract Nostocales ASVs

# Get the sequences
extract_kraken_reads.py -k ${out_path}/${sample}.kraken2 \
  -s ${seqs_path}/${sample}_nonchimeras.fasta \
  -o ${seqs_path}/${sample}_nostocales.fasta \
  -r ${out_path}/${sample}.k2report  \
  --include-children \
  -t 503
# Get the list of Nostocales ASV headers
grep ">" ${seqs_path}/${sample}_nostocales.fasta | sed 's/>//g' > ${out_path}/${sample}_nostocales.txt