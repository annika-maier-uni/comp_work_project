#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.output = './results/trimgalore'


process TRIMMING {
    container "quay.io/biocontainers/trim-galore:0.6.7--hdfd78af_0"
    publishDir "${params.output}", mode: 'copy'


    input:
    tuple val(meta), path(reads)           // channel: [ val(meta), [ reads ] ]

    output:
    tuple val(meta), path("*{3prime,5prime,trimmed,val}*.fq.gz"), emit: reads
    tuple val(meta), path("*report.txt")                        , emit: log     , optional: true
    tuple val(meta), path("*.html")                             , emit: html    , optional: true
    tuple val(meta), path("*.zip")                              , emit: zip     , optional: true


    //trim_galore [options] <filename(s)>
    script:
    def prefix = "${meta.sample}"
    """
    trim_galore --paired ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz --fastqc
    """


}
