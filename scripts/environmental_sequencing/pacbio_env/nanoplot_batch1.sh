#!/bin/bash

#SBATCH --array=2-385
#SBATCH --mem-per-cpu=8G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/environmental_sequencing/pacbio_env/nanoplot_batch1_%A_%a.out
#SBATCH --error=log/environmental_sequencing/pacbio_env/nanoplot_batch1_%A_%a.err
#SBATCH --partition=scavenger

# Load NanoPlot conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate nanoplot
# # Variable witht he sample name
sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p misc_files/environmental_sequencing/pacbio_env/revio_order_10070_barcode.csv | cut -d, -f2)
# Run FastQC
NanoPlot --fastq analyses/reads/pacbio_env/${sample}.fastq.gz \
 --outdir analyses/environmental_sequencing/pacbio_env/qc/nanoplot \
 -p ${sample}_nanoplot_