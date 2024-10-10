#!/usr/bin/env nextflow

// Enable DSL2 syntax
nextflow.enable.dsl=2




process SAMPLESHEET_VALIDATION {


    input:
    path python
    path samplesheet

    output:
    stdout

    script:
    """
    python $python $samplesheet
    """
}