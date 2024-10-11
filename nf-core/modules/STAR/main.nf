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
# Nextflow Process: INDEX_FILE and STAR_ALIGN
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Description:
# The INDEX_FILE process generates a STAR index from a reference genome and GTF file,
# which is required for alignment in the STAR_ALIGN process.
#
# INDEX_FILE:
# Input:
# - path(genome): Path to the reference genome FASTA file.
# - path(gtf): Path to the GTF file containing gene annotations.
#
# Output:
# - path "index_dir/": Directory containing the STAR index files.
#
# STAR_ALIGN:
# Input:
# - tuple path(fasta1), path(fasta2): Paired-end FASTQ files for alignment.
# - path(indexDir): Path to the directory containing the STAR index files (output from INDEX_FILE).
# - path(gtf): Path to the GTF file containing gene annotations.
#
# Output:
# - path "*.bam": BAM file with aligned reads, sorted by coordinate.
# - path "versions.yml": YAML file recording the versions of STAR, Samtools, and GAWK used in the process.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


process INDEX_FILE {

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:ded3841da0194af2701c780e9b3d653a85d27489-0' :
        'biocontainers/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:ded3841da0194af2701c780e9b3d653a85d27489-0' }"

    publishDir "${params.outdir}/STAR/INDEX", mode: 'copy'

    input:
    path genome
    path gtf

    output:
    path "index_dir/", emit: index

    script:
    """
    mkdir "index_dir"
    STAR --runMode genomeGenerate \\
      --genomeDir index_dir/ \\
      --genomeFastaFiles ${genome}


    """
}


process STAR_ALIGN {

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:ded3841da0194af2701c780e9b3d653a85d27489-0' :
        'biocontainers/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:ded3841da0194af2701c780e9b3d653a85d27489-0' }"

    publishDir "${params.outdir}/STAR/ALIGN", mode: 'copy'

    input:
    tuple path(fasta1), path(fasta2)
    path indexDir
    path gtf

    output:
    path("*.bam"), emit: bam
    path  "versions.yml"                                        , emit: versions  // Output: version file

    script:
    """
    STAR \\
        --runThreadN ${params.max_cpus} \\
        --readFilesIn ${fasta1} ${fasta2} \\
        --readFilesCommand zcat \
        --sjdbGTFfile ${gtf} \\
        --outSAMtype BAM SortedByCoordinate \\
        --genomeDir ${indexDir}

    # Capture the versions of the tools used in the process for reproducibility
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        gawk: \$(echo \$(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')
    END_VERSIONS
    """
}