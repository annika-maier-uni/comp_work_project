/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Load the required modules for FastQC, trimming, and HISAT2 alignment
include { FASTQC } from './modules/FASTQ/'
include { TRIMMING } from './modules/trimgalore/'
include { HISAT2_BUILD } from './modules/hisat2/align/'
include { HISAT2_ALIGN } from './modules/hisat2/align/'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Define default paths for samplesheet and reference genome if not provided by the user
params.samplesheet = params.samplesheet ?: "./data/samplesheet.csv"
params.reference  = params.reference  ?: "./data/genome.fa"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    // 1. Reading the samplesheet and preparing data channels

    // Create a channel from the samplesheet file
    reads_channel = Channel
        .fromPath(params.samplesheet)   // Load samplesheet
        .splitCsv(header: true)         // Split the CSV into rows, using the first row as header
        .map { row ->                   // Map each row to structured format for further processing
            [
                ["sample": row.sample, "strandedness": row.strandedness],
                [file(row.fastq_1), file(row.fastq_2)]
            ]
        }
        //.view()  // Debug view for channel content

    // Channel for the reference genome file
    reference_channel = Channel
        .fromPath(params.reference)   // Load reference genome

    // Create a channel for FASTQ files
    tuple_channel = Channel
        .fromPath(params.samplesheet)   // Load samplesheet
        .splitCsv(header: true)         // Split the CSV into rows
        .map { row ->                   // Map each row to structured format with FASTQ file paths
            [
                ["sample": row.sample, "strandedness": row.strandedness],
                [file(row.fastq_1)],     // Path to FASTQ 1
                [file(row.fastq_2)]      // Path to FASTQ 2
            ]
        }
        //.view()  // Debug view for channel content

    // 2. Perform quality control on the FASTQ files using FastQC
    FASTQC(reads_channel)

    // 3. Trim the FASTQ files for quality using Trimgalore
    TRIMMING(reads_channel)
    reads = TRIMMING.out.reads

    reads
        .map { meta, fastq -> fastq }
        .map {fastq1,fastq2 -> tuple([fastq1],[fastq2])}
        .view()
        .set{trimmed_reads}


    // 4. Align reads using HISAT2
    // First, build the HISAT2 index for the reference genome
    output = HISAT2_BUILD(reference_channel)

    // Then, align the reads with the built index
    HISAT2_ALIGN(output,trimmed_reads)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
