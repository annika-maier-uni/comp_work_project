#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.output = './results/fastqc'

process FASTQC {


    container "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
      
      // TODO
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    // 'https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0' :
    // 'biocontainers/fastqc:0.12.1--hdfd78af_0' }"
    // publishDir "${params.output}", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    path("*_fastqc.html")
    path("*_fastqc.zip")
    //versions 

    script:
    """
    fastqc ${reads}
    """
}

