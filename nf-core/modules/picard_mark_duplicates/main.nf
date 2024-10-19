/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simple RNA Sequencing Pipeline
# Author: Weronika Jaśkowiak, Maike Nägele, Tabea Attig, Annika Maier
# Date: 11.10.2024
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Nextflow Process: PICARD_MARKDUPLICATES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Description:
# The PICARD_MARKDUPLICATES process identifies and marks duplicate reads
# in a BAM file generated from an alignment step. 
#
#
# Input:
# - path(input_sam): Path to the input SAM/BAM file from which duplicates need to be marked.
#
# Output:
# - path "marked_duplicates.bam": The output BAM file containing reads with duplicates marked,
#                                 which can be used in subsequent analysis steps.
#
# - path "metrics.txt": A metrics file generated by Picard detailing the number and percentage
#                       of duplicate reads identified in the input file.
#
# - path "versions.yml": YAML file containing version information of
#                        Picard used in the duplication marking process.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process PICARD_MARKDUPLICATES {

    publishDir "${params.outdir}/picard", mode: 'copy'

    // using offical docker container from broadinstitute - RNAseq container not accessible for the Python version we use
    container 'docker.io/broadinstitute/picard'

    input:
    path input_sam

    output: 
    path "marked_duplicates.bam"
    path "metrics.txt"
    path "versions.yml"          , emit: versions

    script:
    """
    # Run Picard MarkDuplicates to mark duplicate reads in the input SAM file
    java -Xmx12g -jar /usr/picard/picard.jar MarkDuplicates \
        -I ${input_sam} \
        -O marked_duplicates.bam \
        -M metrics.txt

    # Capture the versions of the tools used in the process for reproducibility
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(java -Xmx12g -jar /usr/picard/picard.jar MarkDuplicates --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """    
}
