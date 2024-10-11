

process SAMTOOLS_SORT_AND_INDEX {

    publishDir "${params.outdir}/ssmtools", mode: 'copy'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_0' :
        'biocontainers/samtools:1.20--h50ea8bc_0' }"

    input:
    path input_sam

    output:
    path "sorted.sam"             , emit: sam
    path  "versions.yml"          , emit: versions


    script:

    """
    samtools sort ${input_sam} -o sorted.sam \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}