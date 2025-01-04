# Sequencing of 16S from *Peltigera* environmental substrates and selected *Peltigera* specimens

## Read sorting and QC reports

The PacBio reads are already demultiplexed in the raw data, so we just need to sort and label them by sample:
```sh
sbatch scripts/environmental_sequencing/pacbio_env/sort_and_label_reads_batch1.sh
sbatch scripts/environmental_sequencing/pacbio_env/sort_and_label_reads_batch2.sh
sbatch scripts/environmental_sequencing/pacbio_env/sort_and_label_reads_batch3.sh
```
***Note***: batch3 includes (i) the 16S libraries from 308 *Peltigera* specimens that represent the diversity of *Nostoc* we found with *rbcLX* sequences, and (ii) the 16S libraries from 13 env samples.


Now we will get QC reports from FastQC and nanoplot and summarize them with MultiQC

```sh
# Output directories
mkdir -p analyses/environmental_sequencing/pacbio_env/qc/fastqc
mkdir -p analyses/environmental_sequencing/pacbio_env/qc/nanoplot
mkdir -p analyses/environmental_sequencing/pacbio_env/qc/multiqc
# QC reports
sbatch scripts/environmental_sequencing/pacbio_env/fastqc_batch1.sh
sbatch scripts/environmental_sequencing/pacbio_env/fastqc_batch2.sh
sbatch scripts/environmental_sequencing/pacbio_env/fastqc_batch3.sh
sbatch scripts/environmental_sequencing/pacbio_env/nanoplot_batch1.sh
sbatch scripts/environmental_sequencing/pacbio_env/nanoplot_batch2.sh
sbatch scripts/environmental_sequencing/pacbio_env/nanoplot_batch3.sh
# Summarize
conda activate multiqc
multiqc analyses/environmental_sequencing/pacbio_env/qc/fastqc \
 analyses/environmental_sequencing/pacbio_env/qc/nanoplot \
 --outdir analyses/environmental_sequencing/pacbio_env/qc/multiqc
```
***Notes***: 
- FastQC failed on samples s1_3t, s10_23t, and s1_12t.
- NanoPlot failed on samples

## Reads to OTU table

Trim the primers, orient the reads forward and add the sample names to the sequence headers:
```sh
sbatch scripts/environmental_sequencing/pacbio_env/trim_orient_relabel_batch1.sh
sbatch scripts/environmental_sequencing/pacbio_env/trim_orient_relabel_batch2.sh
sbatch scripts/environmental_sequencing/pacbio_env/trim_orient_relabel_batch3.sh
```
Call ASVs using UNOISE3 implemente in Vsearch:

```sh
sbatch scripts/environmental_sequencing/pacbio_env/unoise_asv_call.sh
```

Remove contamination from OTU table:

```sh
sbatch scripts/environmental_sequencing/pacbio_env/decontaminate_batch123.R
```

Normalize otu_tables and plot rarefaction curves

```sh
sbatch scripts/environmental_sequencing/pacbio_env/normalize_otu_tables.R
```


Extracting Nostocales reads

```sh
sbatch scripts/environmental_sequencing/pacbio_env/kraken_pipeline_unoise.sh
sbatch scripts/environmental_sequencing/pacbio_env/kraken_pipeline_unoise_batch123.sh
```

All-to-all BLASTn of Nostocales ASVs

```sh
sbatch scripts/environmental_sequencing/pacbio_env/blast_all_to_all_nostocales_batch123.sh
```

BLASTn of Nostocales ASVs to reference dataset of *Nostoc* 16S sequences

```sh
sbatch scripts/environmental_sequencing/pacbio_env/blast_nostoc_ref_batch123.sh
```

## Curation of 16S from selected *Peltigera* specimens

The goal of the curation was to determine the 16S ASV that corresponded to the *Nostoc* ID determined with *rbcLX*. 

We first summarized the lichen 16S data by getting the top 10 most abundant (i.e., # of mapped reads) nostocalean ASVs in the 308 lichen samples that we sequenced. We also printed a table witht the *Nostoc* IDs that we assigned to those samples based on *rbcLX*:

```sh
mkdir -p analyses/lichen_sequencing/pacbio_lichen
Rscript scripts/lichen_sequencing/pacbio_lichen/prep_for_lichen_16s_curation.R
```

Then, we manually curated the lichen 16S data using (i) the results of the BLASTn seach against the *Nostoc* 16S reference sequences, (ii) the [results of EPA placement of all Nostocales ASVs on the T-BAS tree of *Nostoc*](https://tbas.cifr.ncsu.edu/tbas2_3/genetree.php?runnumber=PLU76GM7), and (iii) the *Nostoc* IDs assigned with *rbcLX*. All libraries had 16S reads that mapped to many Nostocales ASVs. This is expected because the DNA was extracted from specimens that were not surface sterilized. However, with one exception, the most abundant ASV in each sample corresponded to the taxon determined with *rbcLX*. In addition, the top ASV in each library was typically one or two orders of magnitude more abundant than the rest. This suggests that the other mappings correspond to expected epiphytic *Nostoc* or low-level contamination. The curated IDs are in `analyses/lichen_sequencing/pacbio_lichen/selected_specimens_id_table_curated.csv`

*Note*: We retrieved a newick version of the T-BAS tree with the placement and stored it in `analyses/environmental_sequencing/pacbio_env/unoise/taxonomy/tbas/tbas_nostocales_asvs.tree`.

## All-to-all BLASTn of Nostocales ASVs

```sh
sbatch scripts/environmental_sequencing/pacbio_env/blast_all_to_all_nostocales_batch123.sh
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
