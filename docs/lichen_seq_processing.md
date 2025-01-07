# Sequencing of nrITS and *rbcLX* with ONT

## Read QC, demultiplexing and trimming of primers and barcodes

Make a working copy and concatenate the reads for each pool:

```sh
mkdir -p analyses/reads/ont_lichen
sbatch scripts/lichen_sequencing/ont_lichen/cat_reads.sh
```

Get QC report for ONT reads:

```sh
sbatch scripts/lichen_sequencing/ont_lichen/nanoplot_read_qc.sh
```

Prepare the barcode files for demultiplexing the ONT reads:

```sh
Rscript scripts/lichen_sequencing/ont_lichen/sort_lichen_barcodes.R
```

Demultiplex the reads searching separately for ITS and *rbcLX* and looking for reads in both orientations:

```sh
mkdir -p analyses/lichen_sequencing/ont_lichen/demultiplex/its
mkdir -p analyses/lichen_sequencing/ont_lichen/demultiplex/rbclx
sbatch scripts/lichen_sequencing/ont_lichen/demultiplex_its.sh
sbatch scripts/lichen_sequencing/ont_lichen/demultiplex_rbclx.sh
```

Trim the barcodes and primers from the demultiplexed reads using cutadapt. We also filtered the reads by length to remove chimeras:

```sh
# Tables with barcode assignments for all samples
cat analyses/lichen_sequencing/ont_lichen/demultiplex/*its*.tsv > analyses/lichen_sequencing/ont_lichen/demultiplex/all_its_barcodes_plus.tsv
cat analyses/lichen_sequencing/ont_lichen/demultiplex/*rbclx*.tsv > analyses/lichen_sequencing/ont_lichen/demultiplex/all_rbclx_barcodes_plus.tsv
# Trim reads
sbatch scripts/lichen_sequencing/ont_lichen/trim_its_primers_and_barcodes.sh
sbatch scripts/lichen_sequencing/ont_lichen/trim_rbclx_primers_and_barcodes.sh
```

Remove the unbinned reads:

```sh
rm analyses/lichen_sequencing/ont_lichen/demultiplex/its/*UNKNOWN*
rm analyses/lichen_sequencing/ont_lichen/demultiplex/rbclx/*UNKNOWN*
```

## Read error correction

We used a strategy originally described by [Xuan et al. (2023)](https://www.biorxiv.org/content/10.1101/2023.06.19.544637v1) to correct 16S rDNA reads amplified with Unique Molecular Identifiers (UMIs), which are analogous to the barcodes we used to tag the ITS and *rbcLX* amplicons from the lichen symbionts. For each set of sequences, we first conducted a blast search against a database of ITS and rbcLX sequences from alberta from Pardo-De la Hoz et al and removed reads with <97% identity and <400 (ITS) or (700) alignment length. The purpose of this was to remove reads from organisms other than Peltigera or lichenized Nostoc, which can be amplified ocasionally with our priemers. Then, we find a centroid sequence with the USEARCH algorithm using a 95% identity threshold. This centroid sequence is then used as a seed for three rounds of Racon polishing, followed by two rounds of polishing with Medaka:

*Removing the last step of racon polishing decreased the number of unique haplotypes from 707 to 598 in ITS and from 617 to 458 in rbcLX

```sh
# Output directories
mkdir -p analyses/lichen_sequencing/ont_lichen/error_correction/its
mkdir -p analyses/lichen_sequencing/ont_lichen/error_correction/rbclx

# Make blast database with the ITS and rbclxsequences from Pardo-de la Hoz et al. 2024
conda activate blast
makeblastdb -in data/sequences/its_abmi_pardodelahoz2024.fasta \
    -dbtype nucl \
    -parse_seqids \
    -out analyses/lichen_sequencing/ont_lichen/error_correction/its/its_abmi_pardodelahoz2024
makeblastdb -in data/sequences/rbclx_abmi_pardodelahoz2024.fasta \
    -dbtype nucl \
    -parse_seqids \
    -out analyses/lichen_sequencing/ont_lichen/error_correction/rbclx/rbclx_abmi_pardodelahoz2024

# Run the eror correction pipeline
sbatch scripts/lichen_sequencing/ont_lichen/polish_its.sh
sbatch scripts/lichen_sequencing/ont_lichen/polish_rbclx.sh
```

After that, we concatenated all the polished sequences and the table of read counts used for polishing:

```sh
mkdir -p analyses/lichen_sequencing/ont_lichen/consensus
# Compile sequences
cat analyses/lichen_sequencing/ont_lichen/error_correction/its/*polished.fasta > analyses/lichen_sequencing/ont_lichen/consensus/lichen_its.fasta
cat analyses/lichen_sequencing/ont_lichen/error_correction/rbclx/*polished.fasta > analyses/lichen_sequencing/ont_lichen/consensus/lichen_rbclx.fasta
# Compile read counts
echo "sample,n_its_reads" > analyses/lichen_sequencing/ont_lichen/consensus/its_head.csv
echo "sample,n_rbclx_reads" > analyses/lichen_sequencing/ont_lichen/consensus/rbclx_head.csv
cat analyses/lichen_sequencing/ont_lichen/consensus/its_head.csv analyses/lichen_sequencing/ont_lichen/error_correction/its/*read_count.csv > analyses/lichen_sequencing/ont_lichen/consensus/its_read_counts.csv
cat analyses/lichen_sequencing/ont_lichen/consensus/rbclx_head.csv analyses/lichen_sequencing/ont_lichen/error_correction/rbclx/*read_count.csv > analyses/lichen_sequencing/ont_lichen/consensus/rbclx_read_counts.csv
rm analyses/lichen_sequencing/ont_lichen/consensus/its_head.csv
rm analyses/lichen_sequencing/ont_lichen/consensus/rbclx_head.csv
```

## Symbiont IDs and Haplotype clustering

We collapsed the consensus sequences into haplotypes with vsearch and then did a BLASTn against the ITS and *rbcLX* sequences from Alberta that we published in my paper on *Nostoc* phylogenomics:

```sh
sbatch scripts/lichen_sequencing/ont_lichen/collapse_and_blast_its.sh
sbatch scripts/lichen_sequencing/ont_lichen/collapse_and_blast_rbclx.sh
```

Then we added the IDs to the voucher table for manual curation.

```sh
# TO-DO: Rerun using updated data s1c
Rscript scripts/lichen_sequencing/ont_lichen/ids_to_voucher.R
```

***Manual curation notes***: Mark PA2065-PA2072 (P2073-P2080), and PA1561-PA1568 (PA1569-PA1576) as contaminants (pipetting error). I kept records If they had either >10 reads or if the same haplotype was found more than once and at least one of the records had >10 reads. ITS ids were based on blast to database of ABMI sequences from Pardo-De la Hoz et al. 2025. rbcLX IDs we based on blast against database of ABMI sequences from Pardo-De la Hoz et al. 2025 and verified with [T-BAS placement on the phylogenomic tree of *Nostoc*](https://tbas.cifr.ncsu.edu/tbas2_3/genetree.php?runnumber=OS6RU5JN).

After ID curation, I noticed that the *rbcLX* consensus sequences have 5' and 3' tails of 8 bp that don't seem to correspond to biological info. Therefore, I aligned all the consensus sequences and trimmed them manually in Mesquite. I also removed obvious sequencing mistakes, i.e., in-frame insertions in coding regions.

I also used ITSx to extract the ITS region (i.e., ITS1, 5.8S, and ITS2) prior to a more refined haplotype collapsing:

```sh
sbatch scripts/lichen_sequencing/ont_lichen/itsx.sh
```

After that, I re-did the haplotype collapsing but this time including the sequences from the ABMI sampling (Pardo-De la Hoz et al.) in addition to the ones we sequenced with ONT:

```sh
sbatch scripts/lichen_sequencing/multiscale/collapse_multiscale_its_rbclx.sh
```

***Notes:*** I reran ITSx and removed duplicates of P8337 and P10089 from the ITS ABMI seqs because I realized they are chimeras. I also expanded the ABMI dataset from Pardo-De la Hoz et al. by adding 392 ITS sequences from *Peltigera* and 10 *rbcLX* sequenced of *Nostoc*. The majority of the ITS sequences were generated as part of the previous study but were not published in it because the *Nostoc* sequencing had failed for those specimens. Only a few of the specimens (PA4100s) were sequenced in this study as part of the ONT lichen sequencing. This haplotype clustering only considered ITS sequences with > 450 bp and *rbcLX* sequences with > 700 bp

## Match *rbcLX* haplotypes with 16S ASVs

I matched the 16S ASVs to the revised *rbcLX* haplotypes and extrapolated the correspondence to specimens  with identical haplotypes in the local dataset. 

```sh
Rscript scripts/lichen_sequencing/multiscale/match_local_to_haplotypes_asvs.R
```

After this, we had correspondence between *rbcLX* IDs and 16S ASVs for 2253 out 2318 *Nostoc*. This script will also print a revised version of the voucher table for the local dataset (`analyses/lichen_sequencing/multiscale/pelt_detection.RData`), which includes the matched ASVs and the revised ITS and *rbcLX* haplotypes. I used this version of the voucher table for the detection plot. Note that there were 5 cases where a single *rbcLX* haplotype matched with two different 16S ASVs. Nevertheless, the matching ASVs belonged to the same clade as the *rbcLX* haplotype.