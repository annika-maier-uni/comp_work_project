#!/bin/bash -ue
trim_galore --paired SRR6357073_1.fastq.gz SRR6357073_2.fastq.gz --fastqc  

# Capture version information for Trim Galore and Cutadapt, and write them to versions.yml
cat <<-END_VERSIONS > versions.yml
    "TRIMMING":
        trimgalore: $(echo $(trim_galore --version 2>&1) | sed 's/^.*version //; s/Last.*$//')
        cutadapt: $(cutadapt --version)
END_VERSIONS
