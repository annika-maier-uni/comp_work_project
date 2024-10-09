#!/bin/bash -ue
# this creates name.1.ht - name.8.ht
less ./genome.fa
hisat2-build ./genome.fa hisat_build
