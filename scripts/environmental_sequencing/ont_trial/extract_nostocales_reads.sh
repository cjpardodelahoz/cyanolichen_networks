#!/bin/bash

#SBATCH --array=1-16,18-23,25-79,81-89,91-117%25
#SBATCH --mem-per-cpu=8G  # adjust as needed
#SBATCH -c 4 # number of threads per process
#SBATCH --output=log/environmental_sequencing/ont_trial/extract_nostocales_reads_%A_%a.out
#SBATCH --error=log/environmental_sequencing/ont_trial/extract_nostocales_reads_%A_%a.err
#SBATCH --partition=scavenger

# Load dependencies
source $(conda info --base)/etc/profile.d/conda.sh
conda activate kraken2
seqkit_path="/hpc/group/bio1/carlos/apps"
fastx_path="/hpc/group/bio1/carlos/apps/fastx/bin"
# Sample variable
sample=$(cat misc_files/environmental_sequencing/ont_trial/all_trial_samples.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)
# Variables with paths
out_path="analyses/environmental_sequencing/ont_trial/kraken/greengenes"
db_path="/hpc/group/bio1/carlos/apps/kraken2/greengenes/greengenes"
reads_path="analyses/environmental_sequencing/ont_trial/demultiplex/reads"
# Extract Nostocales reads in fastq
extract_kraken_reads.py -k ${out_path}/${sample}/${sample}.kraken2 \
  -s ${reads_path}/${sample}.fastq \
  -o ${out_path}/${sample}/${sample}_nostocales.fastq \
  -r ${out_path}/${sample}/${sample}.k2report  \
  --include-children \
  -t 503 \
  --fastq-output
# Extract Nostocales reads in fasta
extract_kraken_reads.py -k ${out_path}/${sample}/${sample}.kraken2 \
  -s ${reads_path}/${sample}.fastq \
  -o ${out_path}/${sample}/${sample}_tmp.fasta \
  -r ${out_path}/${sample}/${sample}.k2report  \
  --include-children \
  -t 503
# Relabel nostocales sequences to remove read crap
cat ${out_path}/${sample}/${sample}_tmp.fasta | \
 ${seqkit_path}/seqkit replace -p .+ -r "${sample}_nostocales_{nr}" > \
 ${out_path}/${sample}/${sample}_nostocales.fasta
# Collapse identical sequences and keep read counts
${fastx_path}/fasta_formatter -i ${out_path}/${sample}/${sample}_nostocales.fasta | \
 ${fastx_path}/fastx_collapser -o ${out_path}/${sample}/${sample}_nostocales_tmp.fasta
# Relabel the fasta sequences by adding sample name as prefix
cat ${out_path}/${sample}/${sample}_nostocales_tmp.fasta | \
 ${seqkit_path}/seqkit replace -p ^ -r "${sample}_nostocales_" > \
 ${out_path}/${sample}/${sample}_nostocales_unique.fasta
# Remove the tmp files
rm ${out_path}/${sample}/${sample}_tmp.fasta
rm ${out_path}/${sample}/${sample}_nostocales_tmp.fasta

