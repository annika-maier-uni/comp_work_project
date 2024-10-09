process TRIMGALORE {

    input:
    tuple val(meta), path(reads)           // channel: [ val(meta), [ reads ] ]

    output:
    tuple val(meta), path("*{3prime,5prime,trimmed,val}*.fq.gz"), emit: reads
    tuple val(meta), path("*report.txt")                        , emit: log     , optional: true
    tuple val(meta), path("*unpaired*.fq.gz")                   , emit: unpaired, optional: true
    tuple val(meta), path("*.html")                             , emit: html    , optional: true
    tuple val(meta), path("*.zip")                              , emit: zip     , optional: true
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    //trim_galore [options] <filename(s)>
    script:
    def prefix = "${meta.id}"
    """
    trim_galore $args --cores 8 --paired --gzip ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz
    """
}
