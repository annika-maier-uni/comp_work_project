#!/bin/bash -ue
# this creates name.1.ht - name.8.ht
less ./modules/hisat2/align/genome.fa
hisat2-build ./modules/hisat2/align/genome.fa hisat_build
