

process STAR_ALIGN {

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:ded3841da0194af2701c780e9b3d653a85d27489-0' :
        'biocontainers/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:ded3841da0194af2701c780e9b3d653a85d27489-0' }"

    publishDir "${params.outdir}/STAR/ALIGN", mode: 'copy'

    input:
    // (trim_fq, IDX.out, gtf)
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


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        gawk: \$(echo \$(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')
    END_VERSIONS
    """
}