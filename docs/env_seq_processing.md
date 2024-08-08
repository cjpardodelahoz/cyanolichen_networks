# Environmental sequencing trial with 16S and ONT reads

## Read QC, demultiplexing and trimming of primers and barcodes

Make a working copy and concatenate the reads:

```sh
sbatch scripts/environmental_sequencing/ont_trial/cat_reads.sh
```

Get QC report for ONT reads:

```sh
conda activate nanoplot
mkdir -p analyses/environmental_sequencing/ont_trial/qc
NanoPlot --fastq analyses/reads/ont_trial/pool1_all.fastq \
 --outdir analyses/environmental_sequencing/ont_trial/qc \
 -p pool1_all_nanoplot
NanoPlot --fastq analyses/reads/ont_trial/pool2_all.fastq \
 --outdir analyses/environmental_sequencing/ont_trial/qc \
 -p pool2_all_nanoplot
```

Prepare the barcode files for demultiplexing the ONT reads by environmental sample:

```sh
Rscript scripts/environmental_sequencing/ont_trial/sort_barcodes.R
```

Demultiplex the reads searching separately for the forward and reverse barcodes:

```sh
mkdir -p analyses/environmental_sequencing/ont_trial/demultiplex/reads
sbatch scripts/environmental_sequencing/ont_trial/demultiplex_raw.sh
```

I tested this on one of the read files from pool1 which was 107 MB in size. When I used the barcode combinations that included the sequence of the 16S primer, there were about 90Mb (85%) of unbinned reads. When I removed the 16S primers from the barcode combinations, there were 74 MB (70%) of unbinned reads. When demultiplexed using only the barcodes (without primers) and searching separately for the forward and reverse barcodes, only 25 MB (24%) of the reads were unbinned.

Now, we will concatenate the reads from the forward and reverse searches in a sample-wise manner:

```sh
for sample in $(cat misc_files/environmental_sequencing/ont_trial/pool1_samples.txt) ; do
 cat analyses/environmental_sequencing/ont_trial/demultiplex/reads/*${sample}_f.fastq \
  analyses/environmental_sequencing/ont_trial/demultiplex/reads/*${sample}_r.fastq > \
  analyses/environmental_sequencing/ont_trial/demultiplex/reads/${sample}_barcoded.fastq
done
for sample in $(cat misc_files/environmental_sequencing/ont_trial/pool2_samples.txt) ; do
 cat analyses/environmental_sequencing/ont_trial/demultiplex/reads/*${sample}_f.fastq \
  analyses/environmental_sequencing/ont_trial/demultiplex/reads/*${sample}_r.fastq > \
  analyses/environmental_sequencing/ont_trial/demultiplex/reads/${sample}_barcoded.fastq
done
```

Next, we will trim the adapter and 16S primers from the demultiplexed reads using cutadapt. We will also filter them by length so that we keep only the full length 16S sequences. There are some samples for which we had amplified a fragment of the 16S using cyanobacteria-specific primers to make sure we didn't have false negatives. However, a preliminary Kraken run showed that there is Nostocales reads in virtually all samples, so the 16S fragment is not useful anymore. This step will also put all reads in the same orientation and remove any cases where there is a mismatch between the barcodes on the 5' and 3' end. 

```sh
sbatch scripts/environmental_sequencing/ont_trial/trim_primers_and_barcodes.sh
```

We can now remove the output from demultiplex where forward and reverse reads are in separate files as well as the barcoded reads.

```sh
rm analyses/environmental_sequencing/ont_trial/demultiplex/reads/*_f.fastq
rm analyses/environmental_sequencing/ont_trial/demultiplex/reads/*_r.fastq
rm analyses/environmental_sequencing/ont_trial/demultiplex/reads/*UNKNOWN*
rm analyses/environmental_sequencing/ont_trial/demultiplex/reads/*barcoded.fastq
```

## Extracting Nostocales reads

We are now going to do taxonomic assignation on the trimmed demultiplexed reads using Kraken and the greengenes database. The script will run both Kraken and Bracken and convert the report mpa format in case we need it later on.

```sh
# File with list of samples
cat misc_files/environmental_sequencing/ont_trial/pool*_samples.txt > \
 misc_files/environmental_sequencing/ont_trial/all_trial_samples.txt
# Run Kraken2/bracken
sbatch scripts/environmental_sequencing/ont_trial/kraken_greengenes.sh

ls *.mpa > mpa_list.txt
bracken_mpa_path=$(cat mpa_list.txt)
combine_mpa.py -i ${bracken_mpa_path} -o combined_mpa.txt
```

After the trimming, samples s13_14b, s8_17b, neg1 and neg3 had too few reads for the Bracken classification, so we will exclude them from the upcoming analyses. These are array indices 17, 24, 80 and 90.

We are now going to use the Kraken2 classficication to extract the reads classified as Nostocales.

```sh
sbatch scripts/environmental_sequencing/ont_trial/extract_nostocales_reads.sh
```

This will print them in both fasta and fastq format. Also, the following samples didn't have any reads classified as Nostocales after trimming:

-s8_11t
-s8_12t
-s8_17t
-s8_11b
-s8_12b
-s8_16b
-neg2
-s13_8b
-s15_8t
-s15_11t
-neg4
-neg5

We will exclude them from the subsequent analyses (sample indices 1,2,7,11,12,16,48,84,93,96,104 and 105).

## Pooling Nostocales reads for error correction

We are going to pool the nostocales reads in three different ways:

-ALL nostocales reads from ALL 5 sites.
-Separately pooling reads from each site.
-Separately by site and layer (i.e., top vs bottom)

We first will get a list of the sites included in the trial. This will exclude the reads from 12 benchmark samples and the 5 negative controls because they didn't have nostocales reads after trimming (see above).

```sh
cat misc_files/environmental_sequencing/ont_trial/all_trial_samples.txt | head -n 105 | sed "s|_.*||" | sort | uniq | tail -n 5 > misc_files/environmental_sequencing/ont_trial/trial_sites.txt
```

Now, let's concatenate the nostocales reads by site:

```sh
# Varaibles with paths
kraken_out_path="analyses/environmental_sequencing/ont_trial/kraken/greengenes"
reads_out_path="analyses/environmental_sequencing/ont_trial/nanoclust/nostocales/reads"
sample_list="misc_files/environmental_sequencing/ont_trial/trial_samples_with_nostocales.txt"
site_list="misc_files/environmental_sequencing/ont_trial/trial_sites.txt"
# Directory for pooled reads
mkdir -p ${reads_out_path}
# Copy of nostocales reads
for sample in $(cat ${sample_list}) ; do
 cp ${kraken_out_path}/${sample}/${sample}_nostocales.fastq \
  ${reads_out_path}/${sample}_nostocales.fastq
done
# Pool ALL samples
cat ${reads_out_path}/s*_*_nostocales.fastq > ${reads_out_path}/all_nostocales.fastq
# Pool by site
for site in $(cat ${site_list}) ; do
 cat ${reads_out_path}/${site}_*_nostocales.fastq > \
  ${reads_out_path}/${site}_nostocales.fastq
done
# Pool by site and layer
for site in $(cat ${site_list}) ; do
 cat ${reads_out_path}/${site}_*t_nostocales.fastq > \
  ${reads_out_path}/${site}_t_nostocales.fastq
 cat ${reads_out_path}/${site}_*b_nostocales.fastq > \
  ${reads_out_path}/${site}_b_nostocales.fastq
done
# Print list of layer filenames
basename -a ${reads_out_path}/s*_*_nostocales.fastq | sed "s|_nostocales.fastq||" > \
 misc_files/environmental_sequencing/ont_trial/site_layer_list.txt
# Remove copies
rm ${reads_out_path}/*[0-9]t_nostocales.fastq ${reads_out_path}/*[0-9]b_nostocales.fastq ${reads_out_path}/*top_nostocales.fastq
```

And now we run NanoCLUST on the pooled nostocales reads:

```sh
sbatch scripts/environmental_sequencing/ont_trial/nanoclust_nostocales_all.sh
sbatch scripts/environmental_sequencing/ont_trial/nanoclust_nostocales_by_site.sh
sbatch scripts/environmental_sequencing/ont_trial/nanoclust_nostocales_by_layer.sh
```

Several of the pooled runs failed:

In preparation for phylogenetic placement, we are going to pool all the corrected consensus sequences inferred by Nanoclust:

```sh
# Path to seqkit for relabeling seqs
seqkit_path="/hpc/group/bio1/carlos/apps"
# Paths for label manipulation
consensus_paths="misc_files/environmental_sequencing/ont_trial/nostocales_consensus_paths.txt"
cluster_prefix="analyses/environmental_sequencing/ont_trial/nanoclust/nostocales/.*pool.*/.*/.*_nostocales/"
cluster_suffix="/consensus_medaka.fasta/consensus.fasta"
pool_prefix="analyses/environmental_sequencing/ont_trial/nanoclust/nostocales/.*pool.*/"
pool_suffix="/cluster.*/consensus_medaka.fasta/consensus.fasta"
out_seqs_path="analyses/environmental_sequencing/ont_trial/placement/nostocales/sequences"
# List of paths to consensus sequences
ls analyses/environmental_sequencing/ont_trial/nanoclust/nostocales/*pool*/*/*_nostocales/cluster*/consensus_medaka.fasta/consensus.fasta > \
 ${consensus_paths}
# Merge all consensus sequences from Nanoclust and relabel them with the pool and cluster info
for consensus in $(cat ${consensus_paths}) ; do
 cluster=$(echo ${consensus} | sed "s|${cluster_prefix}||" | sed "s|${cluster_suffix}||")
 pool=$(echo ${consensus} | sed "s|${pool_suffix}||" | sed "s|${pool_prefix}||")
 cat ${consensus} | ${seqkit_path}/seqkit replace -p .+ -r "${pool}_${cluster}" >> \
  ${out_seqs_path}/all_nostocales_consensus.fasta
done
```

## Phylogenetic placement of reads classfified as Nostocales

### Placement of Nanoclust-corrected reads

A preliminary view of the placements of the Nostocales reads on the tree shows taht the error correction was not very sensitive to the read context. That is, similar consensus sequences were inferred from the three types od read pools that we implemented. We will have a better idea once we have a placement on Nostoc, but this is encouraging so far.

Let's get a list of the sequences that were placed within Nostoc and subset the fasta file to that for placement within the Nostoc tree

```sh
#
Rscript scripts/environmental_sequencing/ont_trial/parse_corrected_nostocales_placements.R
#
seqkit_path="/hpc/group/bio1/carlos/apps"
placement_path="analyses/environmental_sequencing/ont_trial/placement"
mkdir -p ${placement_path}/nostoc/sequences
mkdir -p ${placement_path}/nostoc/alignments
mkdir -p ${placement_path}/nostoc/trees
${seqkit_path}/seqkit grep -f ${placement_path}/nostocales/nostoc_corrected_labels.txt \
 ${placement_path}/nostocales/sequences/all_nostocales_consensus.fasta > \
 ${placement_path}/nostoc/sequences/all_nostoc_consensus.fasta
```

A total of 338 sequences fell within Nostoc. Now we will do a placement of those sequences on the Nostoc tree.

```sh
sbatch scripts/environmental_sequencing/ont_trial/place_corrected_nostoc_on_ref.sh
```

We are also going to run a blast of the corrected Nostoc sequences against a database of curated Nostoc 16S sequences:

```sh
module load NCBI-BLAST/2.7.1
blast_db="analyses/placement_refs/blast_dbs/16s_db"
blast_out="analyses/environmental_sequencing/ont_trial/placement/nostoc/blast"
mkdir -p${blast_out}
blastn -query analyses/environmental_sequencing/ont_trial/placement/nostoc/sequences/all_nostoc_consensus.fasta \
 -db ${blast_db}  \
 -outfmt '6 qseqid sseqid length nident gaps pident' \
 -max_target_seqs 1 > \
 ${blast_out}/all_nostoc_consensus_blast.txt
```



### Placement of uncorrected raw reads classfified as Nostocales

Let's first create file with a list of samples for which there were reads classfified as nostocales by removing the ones we found to not have any after trimming:

```sh
sed -e '1d;2d;7d;11d;12d;16d;17d;24d;48d;80d;84d;90d;93d;96d;104d;105d' \
 misc_files/environmental_sequencing/ont_trial/all_trial_samples.txt > \
 misc_files/environmental_sequencing/ont_trial/trial_samples_with_nostocales.txt
```

Now, let's pool all reads classified as Nostocales into a single file:

```sh
# Variables with pathts
placement_path="analyses/environmental_sequencing/ont_trial/placement"
sample_list="misc_files/environmental_sequencing/ont_trial/trial_samples_with_nostocales.txt"
kraken_out_path="analyses/environmental_sequencing/ont_trial/kraken/greengenes/"
# Directory for placement
mkdir -p ${placement_path}/nostocales/sequences
mkdir -p ${placement_path}/nostocales/alignments
mkdir -p ${placement_path}/nostocales/trees
# Pool nostocales seqs
for sample in $(cat ${sample_list}) ; do
 cat ${kraken_out_path}/${sample}/${sample}_nostocales_unique.fasta >> \
 ${placement_path}/nostocales/sequences/all_unique_nostocales.fasta
done
```

A quick `grep ">" ${placement_path}/nostocales/sequences/all_unique_nostocales.fasta | wc -l` indicates that we have 95041 sequences in this file.

```sh
sbatch scripts/environmental_sequencing/ont_trial/place_nostocales_on_ref_16s.sh
```



NanoCLUST

```sh
sbatch scripts/environmental_sequencing/ont_trial/nanoclust_a.sh
sbatch scripts/environmental_sequencing/ont_trial/nanoclust_b.sh
 #--reads "${base_dir}/analyses/environmental_sequencing/ont_trial/demultiplex/reads/${sample}.fastq" \
```

rsync -av cjp47@dcc-login.oit.duke.edu:/hpc/group/bio1/carlos/nostoc_communities/analyses/environmental_sequencing/ont_trial/nanoclust/nostocales/site_pools analyses/environmental_sequencing/ont_trial/nanoclust/




```sh
sample="n1_top"
blastn -query ${kraken_out_path}/${sample}/${sample}_nostocales_unique.fasta \
 -db ../nostoc/analyses/species_delimitation/16s/16s_db \
 -outfmt '6 qseqid sseqid length nident gaps pident' \
 -max_target_seqs 1 \
 -out tmp_blast/${sample}.txt
```

```sh
blastn -query analyses/environmental_sequencing/ont_trial/nanoclust/nostocales/site_pools/s15/s15_nostocales/cluster3/consensus_medaka.fasta/consensus.fasta  -db ../nostoc/analyses/species_delimitation/16s/16s_db  -outfmt '6 qseqid sseqid length nident gaps pident' -max_target_seqs 40
```
