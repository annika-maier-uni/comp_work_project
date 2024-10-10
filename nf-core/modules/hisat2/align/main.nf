#!/usr/bin/env nextflow

// Enable DSL2 syntax
nextflow.enable.dsl=2


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Nextflow Process: HISAT2_BUILD and HISAT2_ALIGN
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Description:
# The HISAT2_BUILD process generates a HISAT2 index from
# a reference genome., and the HISAT2_ALIGN process aligns
# paired-end FASTA/FASTQ reads to the generated index.
#
# HISAT2_BUILD:
# Input:
# - path(reference): Reference genome file to build the index.
#
# Output:
# - path "${name}.*.ht2": HISAT2 index files generated from
#   the reference genome (name.1.ht2 to name.8.ht2).
#
# HISAT2_ALIGN:
# Input:
# - path files: HISAT2 index files created from the HISAT2_BUILD process.

# - tuple val(meta), path(fasta1), path(fasta2):
#   Paired-end FASTA/FASTQ files for alignment.
#
# Output:
# - path "output_hisat2_aligned_sam_file.sam": Aligned SAM
#   file resulting from the alignment of paired-end reads.

# - path "versions.yml"
#   recording the versions of HISAT2 and Samtools used in
#   the alignment process.
#
# Author: Annika, Maike, Tabea, Weronika
# Date: 09.10.2024
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


process HISAT2_BUILD {

    // Use Conda environment
    conda "${moduleDir}/environment.yml"

    // Select container based on the engine (Singularity or Docker)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2cdf6bf1e92acbeb9b2834b1c58754167173a410-0' :
        'biocontainers/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2cdf6bf1e92acbeb9b2834b1c58754167173a410-0' }"

    // Copy results to output directory
    publishDir "${params.outdir}/build", mode: 'copy'

    input:
    path(reference)  // Input reference genome

    output:
    path "${name}.*.ht2"  // HISAT2 index output

    script:
    name = "hisat_build"  // Index prefix
    println "Running HISAT2_BUILD"
    """
    # Create HISAT2 index files (name.1.ht2 to name.8.ht2)
    ls -l ${reference}
    hisat2-build ${reference} ${name}
    """
}

process HISAT2_ALIGN {

    // Use Conda environment
    conda "${moduleDir}/environment.yml"

    // Select container based on the engine (Singularity or Docker)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2cdf6bf1e92acbeb9b2834b1c58754167173a410-0' :
        'biocontainers/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2cdf6bf1e92acbeb9b2834b1c58754167173a410-0' }"


    // Copy results to output directory
    publishDir "${params.outputdir}/align", mode: 'copy'

    input:
    path files  // HISAT2 index files
    tuple val(meta), path(fastas) // Input FASTA/FASTQ files for alignment 

    output:
    path "output_hisat2_aligned_sam_file.sam", emit: sam_file  // Aligned SAM file
    // TODO - probably have to change it to 
    // path "${meta.sample}_aligned.sam"

    path "versions.yml", emit: versions  // Output versions

    script:
    full_name = files[0].getName()  // Get the index file name
    name = full_name.split("\\.")[0]  // Extract base name

    // Split the `fastas` list into _1 and _2 files
    def fastq1_list = fastas.findAll { it.name.endsWith('_1.fastq.gz') || it.name.endsWith('_1.fq.gz') }.sort()
    def fastq2_list = fastas.findAll { it.name.endsWith('_2.fastq.gz') || it.name.endsWith('_2.fq.gz') }.sort()

    // Combine paired lists into a single string for HISAT2
    def fasta1 = fastq1_list.join(',')
    def fasta2 = fastq2_list.join(',')


     // Extract read group information from the `meta` object
    def RGID = meta.sample  
    def RGSM = meta.sample 
    def RGLB = "${meta.sample}_lib"
    def RGPL = "unknown" // our data does not have the platform information
    def RGPU = "${meta.sample}_unit" 

    """
    # Align paired-end reads with HISAT2
    hisat2 --fast -x ${name} -1 ${fasta1} -2 ${fasta2} -S output_hisat2_aligned_sam_file.sam \
        --rg-id ${RGID} \
        --rg "SM:${RGSM}" \
        --rg "LB:${RGLB}" \
        --rg "PL:${RGPL}" \
        --threads ${params.max_cpus} \
        --rg "PU:${RGPU}"

    # Save HISAT2 and Samtools versions
    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            hisat2: \$(hisat2 --version | grep -o 'version [^ ]*' | cut -d ' ' -f 2)
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}