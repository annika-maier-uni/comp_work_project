#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reads = '../../assets/samplesheet.csv'
params.output = '../../results'


process TRIMMING {
    container "community.wave.seqera.io/library/trim-galore:0.6.10--e1d78c153f940cdf"
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
    trim_galore $args --cores 8 --paired --gzip ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz
    """


}


workflow {
    reads_in = Channel
        .fromPath(params.reads)
        .splitCsv(header: true)
        .map {
            meta, fastq_1, fastq_2 ->
                if (!fastq_2) {
                    throw new IllegalArgumentException("Sample ${meta.id} is missing a paired-end file! Only paired-end files are allowed.")
                } else {
                    return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                }
        }
        .groupTuple()
        .view()
    //TRIMMING(reads_in)

}