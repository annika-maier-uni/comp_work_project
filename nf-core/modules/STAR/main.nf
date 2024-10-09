#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.output = './results/star'
params.cpus = 4
params.memory = '7.GB'

process STAR {
    container "quay.io/biocontainers/star:2.7.11b--h43eeafb_2"
    publishDir "${params.output}", mode: 'copy'

    input:
    tuple val(meta), path(reads)
    path(reference)

    output:
    tuple val("${meta}"), file("*.bam"), emit: bams
    file("software_versions.STAR.txt")


    script:

    """
    STAR --genomeDir ${reference} \
    --readFilesIn ${reads} \
    --runThreadN ${params.cpus}
    """
}
