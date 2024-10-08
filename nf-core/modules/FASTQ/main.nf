#!/usr/bin/env nextflow

process FASTQC{

    debug true

    input:
    val name
    val full_example_reference_genome_path

    output:
    path "${name}.*.ht2"

    script:
    """
    # this creates name.1.ht - name.8.ht
    hisat2-build ${full_example_reference_genome_path} ${name}
    """

}


workflow{

     // hardcoded paths for testing
    annikas_path = "/home/satan/Documents/computational_workflows/comp_work_project/nf-core/"

    example_reference_genome_path = "data/genome.fa"

    full_example_reference_genome_path = annikas_path + example_reference_genome_path

    FASTQC("testName2", full_example_reference_genome_path)
}