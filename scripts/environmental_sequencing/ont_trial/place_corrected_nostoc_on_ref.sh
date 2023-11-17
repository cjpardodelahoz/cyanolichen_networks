#!/bin/bash

#SBATCH --mem-per-cpu=4G  # adjust as needed
#SBATCH -c 8 # number of threads per process
#SBATCH --output=log/environmental_sequencing/ont_trial/place_corrected_nostoc_on_ref.out
#SBATCH --error=log/environmental_sequencing/ont_trial/place_corrected_nostoc_on_ref.err
#SBATCH --partition=common

# Load dependencies
module load RAxML/8.2.12
mafft_path="/hpc/home/cjp47/mafft-7.475-with-extensions/bin"
trimal_path="/hpc/group/bio1/carlos/apps/trimAl/source"
amas_path="/hpc/group/bio1/carlos/apps/AMAS/amas"
# Variables with paths
# Path to wd. Replace with the path where you have the nostoc project dir
# This is because RAxML doesn't like relative paths
wd="/hpc/group/bio1/carlos/nostoc_communities"
placement_path="analyses/environmental_sequencing/ont_trial/placement"
placement_ref_path="analyses/placement_refs/"
# Merge ref with query 16S seqs
cat ${placement_ref_path}/16s_nostoc_ref.fasta \
 ${placement_path}/nostoc/sequences/all_nostoc_consensus.fasta > \
 ${placement_path}/nostoc/sequences/corrected_plus_ref_nostoc.fasta
# Align ref with query rbclx seqs
${mafft_path}/mafft --retree 1 --maxiterate 0 --adjustdirection \
 ${placement_path}/nostoc/sequences/corrected_plus_ref_nostoc.fasta > \
 ${placement_path}/nostoc/alignments/corrected_plus_ref_nostoc_aln.fasta
# Remove gaps from alignment
${trimal_path}/trimal -in ${placement_path}/nostoc/alignments/corrected_plus_ref_nostoc_aln.fasta \
 -out ${placement_path}/nostoc/alignments/corrected_plus_ref_nostoc_aln_ng.fasta \
 -fasta -nogaps
# Fix headers from trimal output
sed -i "s| .*||" ${placement_path}/nostoc/alignments/corrected_plus_ref_nostoc_aln_ng.fasta
# Convert to phylip format
${amas_path}/AMAS.py convert -i ${placement_path}/nostoc/alignments/corrected_plus_ref_nostoc_aln_ng.fasta \
 -f fasta \
 -d dna \
 -u phylip
mv corrected_plus_ref_nostoc_aln_ng.fasta-out.phy ${placement_path}/nostoc/alignments/corrected_plus_ref_nostoc_aln_ng.phy
# Run placement 
raxmlHPC-PTHREADS-SSE3 -f v \
 -s ${placement_path}/nostoc/alignments/corrected_plus_ref_nostoc_aln_ng.phy \
 -t analyses/placement_refs/16s_nostoc_ref.tree \
 -w ${wd}/${placement_path}/nostoc/trees/ \
 -m GTRCAT -n corrected_nostoc_epa_result -T 4
