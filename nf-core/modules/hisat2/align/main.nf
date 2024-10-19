/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 Simple RNA Sequencing Pipeline
# Author: Weronika Jaśkowiak, Maike Nägele, Tabea Attig, Annika Maier
# Date: 11.10.2024
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
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
#
# HISAT2_BUILD:
# Input:
# - path(reference): Reference genome file to build the index.
#
# Output:
# - path "${name}.*.ht2": HISAT2 index files generated from
#   the reference genome (name.1.ht2 to name.8.ht2).
#
#
# HISAT2_ALIGN:
# Input:
# - path files: HISAT2 index files created from the HISAT2_BUILD process.
#
# - tuple val(meta), path(fasta1), path(fasta2):
#   Paired-end FASTA/FASTQ files for alignment.
#
# Output:
# - path "output_hisat2_aligned_sam_file.sam": Aligned SAM file resulting 
#                                              from the alignment of paired-end reads.
#
# - path "versions.yml": recording the versions of HISAT2 and Samtools used in
#                        the alignment process.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process HISAT2_BUILD {

    // Use Conda environment specified in the environment.yml file
    conda "${moduleDir}/environment.yml"

    // Select the appropriate container based on the workflow's container engine
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2cdf6bf1e92acbeb9b2834b1c58754167173a410-0' :
        'biocontainers/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2cdf6bf1e92acbeb9b2834b1c58754167173a410-0' }"

    // Specify the directory to publish the resulting index files
    publishDir "${params.outdir}/build", mode: 'copy'

    input:
    path(reference)  // Input reference genome file for index generation

    output:
    path "${name}.*.ht2"  // Output HISAT2 index files (name.1.ht2 to name.8.ht2)

    script:
    name = "hisat_build"  // Prefix for the generated index files
    println "Running HISAT2_BUILD"  // Log the start of the process
    """
    # List the reference genome file to ensure it exists
    ls -l ${reference}

    # Create HISAT2 index files from the reference genome
    hisat2-build ${reference} ${name}
    """
}

process HISAT2_ALIGN {

    // Use Conda environment specified in the environment.yml file
    conda "${moduleDir}/environment.yml"

    // Select the appropriate container based on the workflow's container engine
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2cdf6bf1e92acbeb9b2834b1c58754167173a410-0' :
        'biocontainers/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2cdf6bf1e92acbeb9b2834b1c58754167173a410-0' }"

    // Specify the directory to publish the resulting alignment files
    publishDir "${params.outdir}/align", mode: 'copy'

    input:
    path files  // HISAT2 index files generated from the HISAT2_BUILD process
    tuple val(meta), path(fastas) // Input paired-end FASTA/FASTQ files for alignment

    output:
    path "output_hisat2_aligned_sam_file.sam", emit: sam  // Aligned SAM file output
    // TODO - probably have to change it to
    // path "${meta.sample}_aligned.sam"

    path "versions.yml", emit: versions  // Output versions file

    script:
    full_name = files[0].getName()  // Get the name of the first index file
    name = full_name.split("\\.")[0]  // Extract the base name from the index file name

    // Identify paired-end FASTA/FASTQ files for alignment
    def fastq1_list = fastas.findAll { it.name.endsWith('_1.fastq.gz') || it.name.endsWith('_1.fq.gz') }.sort() // List of R1 files
    def fastq2_list = fastas.findAll { it.name.endsWith('_2.fastq.gz') || it.name.endsWith('_2.fq.gz') }.sort() // List of R2 files

    // Combine the lists of R1 and R2 files into comma-separated strings for HISAT2
    def fasta1 = fastq1_list.join(',')
    def fasta2 = fastq2_list.join(',')

    // Extract read group information from the meta object for inclusion in the alignment
    def RGID = meta.sample
    def RGSM = meta.sample
    def RGLB = "${meta.sample}_lib" // Library name
    def RGPL = "unknown" // Platform information is not available in the data
    def RGPU = "${meta.sample}_unit" // Unit name

    """
    # Align paired-end reads with HISAT2 using the specified index and input files
    hisat2 --fast -x ${name} -1 ${fasta1} -2 ${fasta2} -S output_hisat2_aligned_sam_file.sam \
        --rg-id ${RGID} \
        --rg "SM:${RGSM}" \
        --rg "LB:${RGLB}" \
        --rg "PL:${RGPL}" \
        --threads ${params.max_cpus} \
        --rg "PU:${RGPU}"

    # Capture the versions of the tools used in the process for reproducibility
    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            hisat2: \$(hisat2 --version | grep -o 'version [^ ]*' | cut -d ' ' -f 2)
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
