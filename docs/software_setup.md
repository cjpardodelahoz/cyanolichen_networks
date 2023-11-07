# Software setup

This document contains instructions to install the required software used for the study.

## **Kraken2 and allies**

```
# Install kraken2 with conda
mamba create -n kraken2 -c conda-forge -c bioconda kraken2 krakentools bracken
# Path to where you want the Greengenes 16S database
greengenes_db_path=/hpc/group/bio1/carlos/apps/kraken2/greengenes/greengenes
# Build 16S database
kraken2-build --db /hpc/group/bio1/carlos/apps/kraken2/greengenes/greengenes --special greengenes
```

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

