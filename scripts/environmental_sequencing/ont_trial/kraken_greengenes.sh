#!/bin/bash

#SBATCH --array=1-117%25
#SBATCH --mem-per-cpu=8G  # adjust as needed
#SBATCH -c 4 # number of threads per process
#SBATCH --output=log/environmental_sequencing/ont_trial/kraken_greengenes_%A_%a.out
#SBATCH --error=log/environmental_sequencing/ont_trial/kraken_greengenes_%A_%a.err
#SBATCH --partition=scavenger

# Load dependencies
source $(conda info --base)/etc/profile.d/conda.sh
conda activate kraken2
# Sample variable
sample=$(cat misc_files/environmental_sequencing/ont_trial/all_trial_samples.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)
# Variables with paths
out_path="analyses/environmental_sequencing/ont_trial/kraken/greengenes"
db_path="/hpc/group/bio1/carlos/apps/kraken2/greengenes/greengenes"
reads_path="analyses/environmental_sequencing/ont_trial/demultiplex/reads"
# Make sample-specific directory for output
mkdir -p ${out_path}/${sample}
# Run kraken2
kraken2 --d ${db_path} \
 --threads 4 \
 ${reads_path}/${sample}.fastq \
 --report ${out_path}/${sample}/${sample}.k2report > \
 ${out_path}/${sample}/${sample}.kraken2
# Run Bracken on Kraken2 results
bracken -d ${db_path} \
 -i ${out_path}/${sample}/${sample}.k2report \
 -o ${out_path}/${sample}/${sample}.bracken \
 -r 1500
# Convert bracken report to mpa
kreport2mpa.py -r ${out_path}/${sample}/${sample}_bracken_species.k2report \
 -o ${out_path}/${sample}/${sample}_bracken.mpa
 
 