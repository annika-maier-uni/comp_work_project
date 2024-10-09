#!/bin/bash -ue
echo null
trim_galore [] --cores 8 --paired --gzip null_1.fastq.gz null_2.fastq.gz
