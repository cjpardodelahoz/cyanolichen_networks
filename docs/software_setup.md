# Software setup

This document contains instructions to install the required software used for the study.

### Conda environments

- **Kraken2 and allies**

```
# Install kraken2 with conda
mamba create -n kraken2 -c conda-forge -c bioconda kraken2 krakentools bracken
# Path to where you want the Greengenes 16S database
greengenes_db_path=/hpc/group/bio1/carlos/apps/kraken2/greengenes/greengenes
# Build 16S database
kraken2-build --db /hpc/group/bio1/carlos/apps/kraken2/greengenes/greengenes --special greengenes
```

- **NanoCLUST**

```sh
# Install nextflow version 22.10.6
mkdir nextflow_22.10.6
wget https://github.com/nextflow-io/nextflow/releases/download/v22.10.6/nextflow
git clone https://github.com/genomicsITER/NanoCLUST.
mamba env create -f scripts/software/nanoclust/nanoclust_env.yml
mkdir db db/taxdb
wget https://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz && tar -xzvf 16S_ribosomal_RNA.tar.gz -C db
wget https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz && tar -xzvf taxdb.tar.gz -C db/taxdb

canu -correct -p corrected_reads -nanopore-raw subset.fastq genomeSize=1.5k stopOnLowCoverage=1 minInputCoverage=2 minReadLength=500 minOverlapLength=200
canu -correct -p corrected_reads -nanopore-raw subset.fastq genomeSize=1.5k minInputCoverage=2 minReadLength=500 minOverlapLength=200
gridOptions="--partition=scavenger"
 
# r10.4.1_e8.2_400bps_hac@v4.2.0
# r1041_e82_400bps_hac_v4.2.0

module load Java/11.0.8 # Requires Java >8
module load NCBI-BLAST/2.12.0-rhel8
conda activate nanoclust
../nextflow_22.10.6/nextflow run main.nf -profile test,conda
```

- **Demultiplex**

```sh
mamba create -n demultiplex python=3.9
git clone https://github.com/jfjlaros/demultiplex
cd demultiplex
pip install .
```

