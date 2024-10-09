#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reads = '/Users/weronikajaskowiak/Documents/GitHub/comp_work_project/nf-core/data/samplesheet.csv'
params.output = '/Users/weronikajaskowiak/Desktop/results'

process FASTQC {
    container "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
    publishDir "${params.output}", mode: 'copy'

    input:
    tuple val(sampleID), path(read1), path(read2)

    output:
    path("*_fastqc.html")
    path("*_fastqc.zip")

    script:
    """
    fastqc ${read1} ${read2}
    """
}

workflow {

    // Define the fastqc input channel
    reads_in = Channel.fromPath(params.reads)
        .splitCsv(header: true)
        .map { row -> [row.sample, file(row.fastq_1), file(row.fastq_2)] }

    // Run the fastqc step with the reads_in channel
    FASTQC(reads_in)

}
