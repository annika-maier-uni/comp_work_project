#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.output_build = './results/hisat2_build'
params.output_align = './results/hisat2_align'

process HISAT2_BUILD{
    container "quay.io/biocontainers/hisat2:2.2.1--h87f3376_5"

    input: 
    path(reference)

    output:
    path "${name}.*.ht2"

    script:

    name = "hisat_build"

    println "Running HISAT2_BUILD"
    """
    # this creates name.1.ht - name.8.ht
    ls -l ${reference}
    hisat2-build ${reference} ${name}
    """
}


process HISAT2_ALIGN {

    debug true
    input:

    path all_fasta1
    path all_fasta2


    output:
    path "output_hisat2_aligned_sam_file.sam"

    script:

    name = "hisat_build"
    """
    echo "running hisat2"

    # We don't accept unpaired reads
    # --fast for testing 
    hisat2 --fast -x name -1 ${all_fasta1} -2 ${all_fasta2} -S output_hisat2_aligned_sam_file.sam
    """


}
