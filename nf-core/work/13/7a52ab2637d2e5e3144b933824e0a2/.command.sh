#!/bin/bash -ue
python validation_samplesheet.py -file samplesheet_8_samples.csv

if grep -q "failed" "validation.txt"; then
    echo "Samplesheet validation failed!"
    exit 1
else
    echo "Samplesheet validation passed!"
fi
