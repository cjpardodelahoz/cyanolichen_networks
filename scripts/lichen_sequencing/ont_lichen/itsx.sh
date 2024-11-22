#!/bin/bash

#SBATCH --mem-per-cpu=4G  # adjust as needed
#SBATCH -c 4 # number of threads per process
#SBATCH --output=log/lichen_sequencing/ont_lichen/itsx.out
#SBATCH --error=log/lichen_sequencing/ont_lichen/itsx.err
#SBATCH --partition=scavenger


# Script to to extract the ITS region (i.e., ITS1, 5.8S, and ITS2) from the lichen ONT sequences

# Load conda dependencies
source $(conda info --base)/etc/profile.d/conda.sh
conda activate itsx

# Output directory
mkdir -p analyses/lichen_sequencing/ont_lichen/consensus/itsx

# Path variables
out_path="analyses/lichen_sequencing/ont_lichen/consensus/itsx"
vsearch_path="/hpc/group/bio1/carlos/apps/vsearch-2.28.1-linux-x86_64/bin"
amas_path="/hpc/group/bio1/carlos/apps/AMAS/amas"

# Extract ITS region with ITSx
ITSx -i analyses/lichen_sequencing/ont_lichen/consensus/lichen_its.fasta \
    -o ${out_path}/itsx \
    --save_regions all \
    -t F \
    --cpu 4
# Remove crap from ITSx
for out in $(ls ${out_path}/*.fasta) ; do
	sed -i "s|\|F.*||" ${out}
done
# Concatenate ITS1, 5.8S, and ITS2 with AMAS.py
${amas_path}/AMAS.py concat \
    -i ${out_path}/*ITS1.fasta \
    ${out_path}/*5_8S.fasta \
    ${out_path}/*ITS2.fasta \
    -f fasta \
    -d dna \
    -t ${out_path%%/itsx}/lichen_its_trimmed.fasta
rm partitions.txt