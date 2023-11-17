# Software setup

This document contains instructions to install the required software used for the study.

## NanoPlot and MultiQC

```sh
mamba create -n multiqc -c bioconda -c conda-forge multiqc
mamba create -n nanoplot -c bioconda nanoplot
```

## R environment and packages

### Use renv to automatically install all dependencies.

I used the [*renv*](https://rstudio.github.io/renv/reference/restore.html) R package to keep track of the dependencies I used. This generes a lock file that contains the versions of all packages used in the project and you can use it to automatically install them. I worked on R version 4.2.2. You can install renv by doing:

```R
install.packages("renv")
```

The lockfile (`renv.lock`) is in the home of this repo. You can install the packages by doing:

```R
renv::restore(lockfile = "renv.lock")
```

### List of package versions I installed

- Tidyverse 1.3.2
- spgs 1.0-4

## Demultiplex

```sh
mamba create -n demultiplex python=3.9
git clone https://github.com/jfjlaros/demultiplex
cd demultiplex
pip install .
```

## Cutadapt

[Cutadapt](https://cutadapt.readthedocs.io/en/stable/) finds and removes adapter sequences, primers, poly-A tails and other types of unwanted sequence from your high-throughput sequencing reads.

```sh
mamba create -n cutadapt -c bioconda -c defaults -c conda-forge cutadapt=4.5
```

## **Kraken2 and allies**

Kraken is a program for k-mer based taxonomic assignation of metagenomic and 16S amplicon reads

```
# Install kraken2 with conda
mamba create -n kraken2 -c conda-forge -c bioconda kraken2 krakentools bracken
# Path to where you want the Greengenes 16S database
greengenes_db_path=/hpc/group/bio1/carlos/apps/kraken2/greengenes/greengenes
# Build 16S database
kraken2-build --db /hpc/group/bio1/carlos/apps/kraken2/greengenes/greengenes --special greengenes
bracken-build -d /hpc/group/bio1/carlos/apps/kraken2/greengenes/greengenes -t 8 -k 35 -l 1500
```

## SeqKit, FASTX

SeqKit and FASTX are toolkits for FASTA/Q manipulation. I am using SeqKit v2.2.0 and FASTX v0.0.13. You can download the executable binaries from https://github.com/shenwei356/seqkit/releases/tag/v2.2.0 and http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2.

## **NanoCLUST**

[NanoCLUST](https://github.com/genomicsITER/NanoCLUST/tree/master) is a pipeline to correct nanopore amplicon reads of 16S. It will **cluster** the reads using HBDSCAN and UMAP projections of the reads. Then, it performs **error correction and consensus inference** using Canu, Racon, and Medaka. Finally, the consensus sequences are **classifified** by blasting against a reference 16S database.

### Dependencies

NanoCLUST is implemented as a [Nextflow](https://nextflow.io/docs/latest/getstarted.html) workflow, and it works with versions 22.10.6 or eralier. Download the binaries from the github repo:

```sh
mkdir nextflow_22.10.6
cd nextflow_22.10.6
wget https://github.com/nextflow-io/nextflow/releases/download/v22.10.6/nextflow
cd ../
```

NanoCLUST is supposed to automatically install dependencies via conda environments or Docker containers. I cannot use the Docker containers because I am working on Duke's computer cluster. I wanted to use the conda dependencies but the automatic installation never happened and I could not figure out why. Instead, I created a single conda environment with *almost* all dependencies. The distribution of NanoCLUST that I downloaded [commit 9364ddc] was using Medaka v.1.0.3, but Medaka v.1.8.0 was available so I wanted to use that instead to take advantage of the more recent error models. To do this, I had to use newer versions of several other dependencies. You can build the environment with these versions using the YAML file in the repo:

```sh
mamba env create -f scripts/software/nanoclust/nanoclust_env.yml
```

The only dependencies that are not part of the environment are Java > v8 and BLAST, which I used from modules in my HPC because the conda environments are buggy. 

### Main pipeline

Now you can download the NanoCLUST pipeline code and the 16S database for the blast searches:

```sh
git clone https://github.com/genomicsITER/NanoCLUST.
mkdir db db/taxdb
wget https://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz && tar -xzvf 16S_ribosomal_RNA.tar.gz -C db
wget https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz && tar -xzvf taxdb.tar.gz -C db/taxdb
```

I made some slight modifications to the `main.nf` file, which has the main code for the pipeline. First, I added a command for the pipeline to pause during the canu error correction until the corrected reads are generated:

```sh
while [ ! -f corrected_reads.correctedReads.fasta.gz ] ; do sleep 1 ; done
```

This solved an issue that was arising because it tried to find the reads for the next step when Canu was still running. I also modified the canu call to make it use a specific partition from Duke's cluster (`gridOptions="--partition=scavenger"`). Finally, I modified the medaka call so it used the correct error model for my ONT reads (`r1041_e82_400bps_hac_v4.2.0`).

With these modifications (except the medaka error model), I ran the NanoCLUST dataset succesfully:

```sh
cd NanoCLUST
cp scripts/software/nanoclust/main.nf . # replace the main.nf file
module load Java/11.0.8 # Requires Java >8
module load NCBI-BLAST/2.12.0-rhel8
conda activate nanoclust
../nextflow_22.10.6/nextflow run main.nf -profile test,conda
cd ../
```

- **Demultiplex**

```sh
mamba create -n demultiplex python=3.9
git clone https://github.com/jfjlaros/demultiplex
cd demultiplex
pip install .
```

