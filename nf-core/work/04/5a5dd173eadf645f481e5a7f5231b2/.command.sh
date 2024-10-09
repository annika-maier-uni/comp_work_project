#!/bin/bash -ue
STAR --genomeDir ./data/genome.fa       --readFilesCommand "gzip -d -c -f"     --readFilesIn SRR6357070_1.fastq.gz SRR6357070_2.fastq.gz     --outSAMmode Full     --outSAMattributes Standard     --outSAMunmapped None     --outReadsUnmapped Fastx     --outFilterMismatchNoverLmax 0.02     --runThreadN params.cpus
