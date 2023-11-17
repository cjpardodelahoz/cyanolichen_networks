#!/bin/bash

#SBATCH --array=1-17%10
#SBATCH --mem-per-cpu=8G  # adjust as needed
#SBATCH -c 16 # number of threads per process
#SBATCH --output=log/environmental_sequencing/ont_trial/nanoclust_b_%A_%a.out
#SBATCH --error=log/environmental_sequencing/ont_trial/nanoclust_b_%A_%a.err
#SBATCH --partition=common

# Load dependencies
source $(conda info --base)/etc/profile.d/conda.sh
conda activate nanoclust
module load Java/11.0.8 # Requires Java >8
module load NCBI-BLAST/2.12.0-rhel8
# Paths to nextflow and nanoclust pipeline
nextflow_path="/hpc/group/bio1/carlos/apps/nextflow_22.10.6"
nanoclust_path="/hpc/group/bio1/carlos/apps/NanoCLUST"
base_dir="/hpc/group/bio1/carlos/nostoc_communities"
# Sample variable
sample=$(cat misc_files/environmental_sequencing/ont_trial/nanoclust_b_samples.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)
# Make sample-specific directory for output
mkdir -p ${base_dir}/analyses/environmental_sequencing/ont_trial/nanoclust/b/${sample}
# Run Nanoclust pipeline
${nextflow_path}/nextflow run -work-dir ${base_dir}/analyses/environmental_sequencing/ont_trial/nanoclust/b/${sample} \
 ${nanoclust_path}/main.nf \
 --reads "${base_dir}/analyses/environmental_sequencing/ont_trial/demultiplex/reads/${sample}.fastq" \
 --db "${nanoclust_path}/db/16S_ribosomal_RNA" \
 --tax "${nanoclust_path}/db/taxdb" \
 --min_read_length 1200 \
 --max_read_length 1600 \
 --outdir "${base_dir}/analyses/environmental_sequencing/ont_trial/nanoclust/b/${sample}"
 
 