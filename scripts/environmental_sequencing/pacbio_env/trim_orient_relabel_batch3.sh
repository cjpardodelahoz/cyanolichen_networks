#!/bin/bash

#SBATCH --array=2-322
#SBATCH --mem-per-cpu=32G  # adjust as needed
#SBATCH -c 1 # number of threads per process
#SBATCH --output=log/environmental_sequencing/pacbio_env/trim_orient_relabel_batch3_%A_%a.out
#SBATCH --error=log/environmental_sequencing/pacbio_env/trim_orient_relabel_batch3_%A_%a.err
#SBATCH --partition=scavenger

# Load dependencies
source $(conda info --base)/etc/profile.d/conda.sh
conda activate cutadapt

# Variable witht he sample name
sample=$(sed -n ${SLURM_ARRAY_TASK_ID}p misc_files/environmental_sequencing/pacbio_env/revio_order_10483_barcode.csv | cut -d, -f2)
# Path variables
vsearch_path="/hpc/group/bio1/carlos/apps/vsearch-2.28.1-linux-x86_64/bin"
reads_path="analyses/reads/pacbio_env"

# Orient the reads
${vsearch_path}/vsearch --orient analyses/reads/pacbio_env/${sample}.fastq.gz \
 --db ${vsearch_path}/gg_12_8_2000.fasta \
 --relabel "sample="${sample}";" \
 --fastqout ${reads_path}/${sample}_oriented.fastq

# Trim primers with cutadapt
cutadapt -g AGTTTGATCCTGGCTCAG...AAGTCGTAACAAGGTAACC \
 --revcomp \
 --discard-untrimmed \
 --minimum-length 1300 \
 --maximum-length 1600 \
 --cores 1 \
 -o ${reads_path}/${sample}_trimmed.fastq \
 ${reads_path}/${sample}_oriented.fastq

# Rempve the oriented reads
rm ${reads_path}/${sample}_oriented.fastq
# Compress the trimmed reads
gzip ${reads_path}/${sample}_trimmed.fastq

