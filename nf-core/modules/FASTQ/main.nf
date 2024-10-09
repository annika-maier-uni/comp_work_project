#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.output = './results/fastqc'

process FASTQC {
    container "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
    publishDir "${params.output}", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    path("*_fastqc.html")
    path("*_fastqc.zip")
    path  "versions.yml"           , emit: versions

    script:
    """
    fastqc ${reads}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed '/FastQC v/!d; s/.*v//' )
    END_VERSIONS
    """
}

