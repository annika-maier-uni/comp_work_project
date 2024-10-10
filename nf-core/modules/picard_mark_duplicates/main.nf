nextflow.enable.dsl = 2

params.output = './results/picard'


process PICARD_MARKDUPLICATES {

    publishDir "${params.output}", mode: 'copy'

    //conda "${moduleDir}/environment.yml"

    container 'docker.io/broadinstitute/picard'
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //      'https://depot.galaxyproject.org/singularity/picard:3.2.0--hdfd78af_0' :
    //      'biocontainers/picard:3.1.1--hdfd78af_0' }"

    input:
    path input_sam

    output: 
    path "marked_duplicates.bam"
    path "metrics.txt"
    path "versions.yml"          , emit: versions

    script:
    """
    java -Xmx12g -jar /usr/picard/picard.jar MarkDuplicates \
        -I ${input_sam} \
        -O marked_duplicates.bam \
        -M metrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(java -Xmx12g -jar /usr/picard/picard.jar MarkDuplicates --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """    
}
