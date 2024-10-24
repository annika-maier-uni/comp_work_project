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
# Nextflow Process: TRIMMING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Description:
# This process utilizes Trim Galore to perform quality trimming on paired-end FASTQ files. 
# It generates quality control reports using FastQC and captures tool 
# versions for reproducibility.
#
#
# Input:
# - Tuple: val(meta), path(reads): Metadata and paired-end FASTQ files to be trimmed.
#
# Output:
# - Tuple: val(meta), path("*{3prime,5prime,trimmed,val}*.fq.gz"): Trimmed FASTQ files with various
#                                                                  suffixes indicating trimming types.
#
# - Tuple: val(meta), path("*report.txt"): Report file generated by Trim Galore, detailing 
#                                          the trimming process.
#
# - Tuple: val(meta), path("*.html"): HTML report file generated by FastQC, 
#                                     summarizing the quality of the trimmed reads.
#
# - Tuple: val(meta), path("*.zip"): Compressed file of the output, if applicable.
#
# - Path: versions.yml: YAML file containing the versions of Trim Galore and Cutadapt used in the process.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process TRIMMING {

    // Conda environment
    conda "${moduleDir}/environment.yml"

    // Select container based on the engine (Singularity or Docker)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trim-galore:0.6.7--hdfd78af_0' :
        'biocontainers/trim-galore:0.6.7--hdfd78af_0' }"

    // Copy results to output directory
    publishDir "${params.outdir}/trimgalore", mode: 'copy'

    input:
    tuple val(meta), path(reads)  // Input: meta info and reads (paired-end FASTQ files)

    output:
    tuple val(meta), path("*{3prime,5prime,trimmed,val}*.fq.gz"), emit: reads  // Output: trimmed reads
    tuple val(meta), path("*report.txt")                        , emit: log     , optional: true  // Output: report (optional)
    tuple val(meta), path("*.html")                             , emit: html    , optional: true  // Output: HTML report (optional)
    tuple val(meta), path("*.zip")                              , emit: zip     , optional: true  // Output: compressed file (optional)
    path  "versions.yml"                                        , emit: versions  // Output: version file

    script:
    // Use the sample name from meta info as a prefix for input files
    def prefix = "${meta.sample}"

    // Initialize empty variables to store adapter and length trimming parameters
    def adapter_trimming = ''
    def length_trimmedreads = ''

    // If adapter trimming is specified in parameters, append the adapter option
    if (params.adapter == adapter_trimming) {
        adapter_trimming = adapter_trimming ?: '--adapter ${params.adapter}'
    }
    // If minimum read length trimming is specified in parameters, append the length option
    if (params.length == length_trimmedreads) {
        length_trimmedreads = length_trimmedreads ?: '--length ${params.length}'
    }

    // Run Trim Galore on paired reads with FastQC, and apply adapter trimming and length trimming if specified
    """
    trim_galore --paired ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz --fastqc ${adapter_trimming} ${length_trimmedreads}

    # Capture version information for Trim Galore and Cutadapt, and write them to versions.yml
    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            trimgalore: \$(echo \$(trim_galore --version 2>&1) | sed 's/^.*version //; s/Last.*\$//')
            cutadapt: \$(cutadapt --version)
    END_VERSIONS
    """
}
