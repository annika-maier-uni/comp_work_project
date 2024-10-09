process HIASAT2_BUILD{

    input: 
    val full_example_reference_genome_path

    output:
    path "${name}.*.ht2"

    script:

    name = name = "hisat_build"

    println "Running HISAT2_BUILD"
    """
    # this creates name.1.ht - name.8.ht
    hisat2-build ${full_example_reference_genome_path} ${name}
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

workflow{


     // hardcoded paths for testing
    annikas_path = "/home/satan/Documents/computational_workflows/comp_work_project/nf-core/"    
    example_reference_genome_path = "data/genome.fa" 
    full_example_reference_genome_path = annikas_path + example_reference_genome_path

    // Output is not necessary I think.
    // hisat2 fetches them based on name
    eight_reference_files = HIASAT2_BUILD(full_example_reference_genome_path)

    // currently these are one file each for testing
    // later I need multiple files
    // Specifically, they need ot be comma separated

    // Comma-separated list of files containing mate 1s 
    // (or read 1s, filename usually includes _1), e.g. 
    // -1 flyA_1.fq,flyB_1.fq. Sequences specified with this option 
    // must correspond file-for-file and read-for-read with those 
    // specified in -2. 
    // Reads may be a mix of different lengths. 

    fasta1 = annikas_path + "data/SRR6357070_1.fastq.gz" 
    fasta2 = annikas_path + "data/SRR6357070_2.fastq.gz" 

    // HISAT2_ALIGN(name, fasta1, fasta2)
}