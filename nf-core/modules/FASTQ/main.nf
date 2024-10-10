#!/usr/bin/env nextflow

// Enable DSL2 syntax
nextflow.enable.dsl=2



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Nextflow Process: FASTQC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Description:
# The FASTQC process runs FastQC on paired-end FASTQ files
# to assess the quality of sequencing data.
#
# Input:
# - tuple val(meta), path(reads): Metadata and paired-end
#   FASTQ files for quality assessment.
#
# Output:
# - path("*_fastqc.html"): FastQC HTML report files
#   generated for each FASTQ file.

# - path("*_fastqc.zip"): FastQC zip files containing
#   detailed results.

# - path "versions.yml"
#   recording the version of FastQC used during the analysis.
#
# Author: Annika, Maike, Tabea, Weronika
# Date: 09.10.2024
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process FASTQC {

    // Use Conda environment
    conda "${moduleDir}/environment.yml"

    // Select container based on the engine (Singularity or Docker)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0' :
        'biocontainers/fastqc:0.12.1--hdfd78af_0' }"

    // Copy results to output directory
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(meta), path(reads)  // Input: meta info and reads (paired-end FASTQ files)

    output:
    path("*_fastqc.html")  // FastQC report
    path("*_fastqc.zip")   // FastQC zip file
    path "versions.yml", emit: versions  // FastQC version info


    script:
    """
    # Run FastQC
    fastqc ${reads} --threads $params.max_cpus --memory $params.max_memory

    # Save HISAT2 and Samtools versions
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed '/FastQC v/!d; s/.*v//' )
    END_VERSIONS
    """
}

