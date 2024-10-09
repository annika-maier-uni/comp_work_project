#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reads = '../../assets/samplesheet.csv'
params.output = '../../results'


process TRIMMING {
    container "biocontainers/trim-galore:0.6.7--hdfd78af_0"
    publishDir "${params.output}", mode: 'copy'


    input:
    tuple val(meta), path(reads)           // channel: [ val(meta), [ reads ] ]

    output:
    tuple val(meta), path("*{3prime,5prime,trimmed,val}*.fq.gz"), emit: reads


    //trim_galore [options] <filename(s)>
    script:
    def prefix = "${meta.id}"
    """
    echo ${prefix}
    trim_galore --paired ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz
    """


}


workflow {
    reads_in = Channel
        .fromPath('../../assets/samplesheet.csv')
        .splitCsv(header: true)
        .map { row -> [["sample": row.sample, "strandedness": row.strandedness],[file(row.fastq_1), file(row.fastq_2)]]}
        .view()

    TRIMMING(reads_in)

}