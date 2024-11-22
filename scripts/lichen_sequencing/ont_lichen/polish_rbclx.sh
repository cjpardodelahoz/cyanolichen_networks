#!/bin/bash

#SBATCH --array=1-2638
#SBATCH --mem-per-cpu=4G  # adjust as needed
#SBATCH -c 4 # number of threads per process
#SBATCH --output=log/lichen_sequencing/ont_lichen/polish_rbclx_%A_%a.out
#SBATCH --error=log/lichen_sequencing/ont_lichen/polish_rbclx_%A_%a.err
#SBATCH --partition=scavenger

# Load conda dependencies
source $(conda info --base)/etc/profile.d/conda.sh
conda activate ssumi
# Path variables
vsearch_path="/hpc/group/bio1/carlos/apps/vsearch-2.28.1-linux-x86_64/bin"
seqkit_path="/hpc/group/bio1/carlos/apps"
reads_path="analyses/lichen_sequencing/ont_lichen/demultiplex/rbclx"
out_path="analyses/lichen_sequencing/ont_lichen/error_correction/rbclx"

# Sample variable
sample=$(cat analyses/lichen_sequencing/ont_lichen/demultiplex/all_rbclx_barcodes_plus.tsv | cut -f 1 | sed -n ${SLURM_ARRAY_TASK_ID}p)
#sample=$(cat analyses/lichen_sequencing/ont_lichen/demultiplex/all_rbclx_barcodes_plus.tsv | cut -f 1 | sed -n 1p)

# READ FILTERING

# Convert fastq to fasta using seqkit and remove anything after the first space in the headers
${seqkit_path}/seqkit fq2fa ${reads_path}/${sample%%_plus}.fastq | sed 's/ .*$//' > ${out_path}/${sample%%_plus}.fasta
# Run a blastn search against the Peltigera rbclx sequences from Pardo-de la Hoz et al. 2024 Then, 
blastn -query ${out_path}/${sample%%_plus}.fasta \
    -db analyses/lichen_sequencing/ont_lichen/error_correction/rbclx/rbclx_abmi_pardodelahoz2024 \
    -outfmt "6 qseqid pident length" \
    -max_target_seqs 1 \
    -num_threads 1 > ${out_path}/${sample%%_plus}_blast.tsv
# Remove the reads from the fastq file with either pident < 97 or length < 500
awk '$2 >= 97 && $3 >= 500 {print $1}' ${out_path}/${sample%%_plus}_blast.tsv > ${out_path}/${sample%%_plus}_blast_ids.txt
${seqkit_path}/seqkit grep -f ${out_path}/${sample%%_plus}_blast_ids.txt ${reads_path}/${sample%%_plus}.fastq -o ${out_path}/${sample%%_plus}_pelt.fastq
# Filter reads by expected error rate
${vsearch_path}/vsearch --fastq_filter ${out_path}/${sample%%_plus}_pelt.fastq \
    --fastq_qmin 1 \
    --fastq_qmax 50 \
    --fastq_ascii 33 \
    --fastq_maxee_rate 0.03 \
    --fastqout ${out_path}/${sample%%_plus}_filtered.fastq
# Get centroid sequences
${vsearch_path}/vsearch --cluster_fast ${out_path}/${sample%%_plus}_filtered.fastq \
    --id 0.95 \
    --sizeout \
    --centroids ${out_path}/${sample%%_plus}_centroids.fasta
# Get seed for Racon
${vsearch_path}/vsearch --sortbysize ${out_path}/${sample%%_plus}_centroids.fasta \
    --topn 1 \
    --relabel seed \
    --output ${out_path}/${sample%%_plus}_seed.fasta
# Copy of the seed for looping
cp ${out_path}/${sample%%_plus}_seed.fasta ${out_path}/${sample%%_plus}_racon3x.fasta

# POLISHING

# Polish with Racon 3 times
for i in `seq 1 3`; do
    minimap2 -t 1 \
        -x map-ont \
        ${out_path}/${sample%%_plus}_racon3x.fasta \
        ${out_path}/${sample%%_plus}_filtered.fastq > \
        ${out_path}/${sample%%_plus}_ovlp.paf
    
    racon -t 1 \
        -m 8 \
        -x -6 \
        -g -8 \
        -w 500 \
        -e 0.04 \
        --no-trimming \
        ${out_path}/${sample%%_plus}_filtered.fastq \
        ${out_path}/${sample%%_plus}_ovlp.paf \
        ${out_path}/${sample%%_plus}_racon3x.fasta > \
        ${out_path}/${sample%%_plus}_tmp.fasta

    mv ${out_path}/${sample%%_plus}_tmp.fasta ${out_path}/${sample%%_plus}_racon3x.fasta

    rm ${out_path}/${sample%%_plus}_ovlp.paf
done

# Copy of the seed for looping
cp ${out_path}/${sample%%_plus}_racon3x.fasta ${out_path}/${sample%%_plus}_polished.fasta
# Polish with 2 rounds of Medaka
for i in `seq 1 2`; do
    medaka_consensus -i ${out_path}/${sample%%_plus}_filtered.fastq \
        -d ${out_path}/${sample%%_plus}_polished.fasta \
        -m r1041_e82_400bps_hac_v4.3.0 \
        -o ${out_path}/medaka_rbclx_${sample}_${i} \
        -t 1
    rm ${out_path}/${sample%%_plus}_polished.fasta.*
    mv ${out_path}/medaka_rbclx_${sample}_${i}/consensus.fasta ${out_path}/${sample%%_plus}_polished.fasta
done
# Remove Medaka temporary files
rm -r ${out_path}/medaka_rbclx_${sample}_*

# Replace header with sample name
sed -i "s/>.*/>${sample%%_rbclx_plus}/" ${out_path}/${sample%%_plus}_polished.fasta

# Print read count used for polishing
n_reads=$(grep -c "^+$" ${out_path}/${sample%%_plus}_filtered.fastq)
echo "${sample%%_rbclx_plus},${n_reads}" > ${out_path}/${sample%%_plus}_read_count.csv