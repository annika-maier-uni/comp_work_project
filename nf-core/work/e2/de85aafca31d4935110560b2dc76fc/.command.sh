#!/bin/bash -ue
STAR --genomeDir ./data/genome.fa     --readFilesIn SRR6357070_1.fastq.gz SRR6357070_2.fastq.gz     --runThreadN 4
