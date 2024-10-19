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
# Nextflow Process: SAMTOOLS_SORT_AND_INDEX
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Description:
# The SAMTOOLS_SORT_AND_INDEX process sorts a SAM file generated from
# an alignment step and then creates an index for the sorted file using
# Samtools. 
#
#
# Input:
# - path(input_sam): Path to the input SAM file that needs to be sorted and indexed.
#
# Output:
# - path "sorted.sam": The sorted SAM file produced by Samtools after 
#                      sorting the input SAM file.
#
# - path "versions.yml": YAML file containing the version information of Samtools
#                        used in the sorting and indexing process.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process SAMTOOLS_SORT_AND_INDEX {

    publishDir "${params.outdir}/samtools", mode: 'copy'

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
    # Sort the SAM file using Samtools and save the output as 'sorted.sam'
    samtools sort ${input_sam} -o sorted.sam \\

    # Capture the versions of the tools used in the process for reproducibility
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
