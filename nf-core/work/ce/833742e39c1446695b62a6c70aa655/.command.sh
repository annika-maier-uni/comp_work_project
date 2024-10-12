#!/bin/bash -ue
# List the reference genome file to ensure it exists
ls -l genome.fa

# Create HISAT2 index files from the reference genome
hisat2-build genome.fa hisat_build
