params.output = './results/STAR/alignment'

process STAR_ALIGN {

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:ded3841da0194af2701c780e9b3d653a85d27489-0' :
        'biocontainers/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:ded3841da0194af2701c780e9b3d653a85d27489-0' }"

    publishDir "${params.output}", mode: 'copy'

    input:
    // (trim_fq, IDX.out, gtf)
    tuple path(fasta1), path(fasta2)
    path indexDir
    path gtf

    output:
    path("*.bam"), emit: align_bam
    path  "versions.yml"                                        , emit: versions  // Output: version file

    script:

    """

    STAR  \\
        --readFilesIn ${fasta1} ${fasta2} \\
        --outSAMtype BAM SortedByCoordinate \\
        --sjdbGTFfile ${gtf} \\
        --genomeDir ${indexDir} \\
        --runThreadN 16 \\
        --outFileNamePrefix 557.

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        gawk: \$(echo \$(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')
    END_VERSIONS
    """
}