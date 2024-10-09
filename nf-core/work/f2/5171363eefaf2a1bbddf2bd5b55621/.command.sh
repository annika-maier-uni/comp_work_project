#!/bin/bash -ue
# this creates name.1.ht - name.8.ht
less ./data/genome.fa
hisat2-build ./data/genome.fa hisat_build
