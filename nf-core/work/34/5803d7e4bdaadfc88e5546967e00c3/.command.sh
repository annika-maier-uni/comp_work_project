#!/bin/bash -ue
# Run FastQC
fastqc SRR6357076_1.fastq.gz SRR6357076_2.fastq.gz --threads 4 --memory 16.GB

# Capture the versions of the tools used in the process for reproducibility
cat <<-END_VERSIONS > versions.yml
"FASTQC":
    fastqc: $( fastqc --version | sed '/FastQC v/!d; s/.*v//' )
END_VERSIONS
