# Environmental sequencing trial with 16S and ONT reads


Make a working copy and concatenate the reads:

```sh
sbatch scripts/environmental_sequencing/ont_trial/cat_reads.sh
```

Get QC report for ONT reads:

```sh
mkdir -p analyses/environmental_sequencing/ont_trial/qc
sbatch scripts/environmental_sequencing/ont_trial/fastqc.sh
```

