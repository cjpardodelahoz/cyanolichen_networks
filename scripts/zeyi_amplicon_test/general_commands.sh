#!/bin/bash

# Install demultiplex
mamba create -n demultiplex python=3.9
git clone https://github.com/jfjlaros/demultiplex
cd demultiplex
pip install .
# Use demultiplex
conda activate demultiplex
demultiplex match -m 1 BARCODES READS