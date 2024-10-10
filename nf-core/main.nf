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
include { PICARD_MARKDUPLICATES } from './modules/picard_mark_duplicates/'
include { SAMTOOLS_SORT_AND_INDEX } from './modules/samtools/'
include { INDEX_FILE } from './modules/STAR/align/create_index_file/'
include { STAR_ALIGN} from './modules/STAR/align/'
include {SAMPLESHEET_VALIDATION} from './modules/samplesheet_validation'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
params.samplesheet = params.samplesheet ?: './data/samplesheet.csv'
params.fasta  = params.fasta  ?: "./data/genome.fa"
params.gtf  = params.gtf  ?: "./data/genes.gtf"
params.python_file = './bin/validation_samplesheet.py'
params.outdir      = params.outdir ?: "./results"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    // 1. Reading the samplesheet and preparing data channels
    samplesheet_channel = Channel
        .fromPath(params.samplesheet)

    python_channel = Channel
        .fromPath(params.python_file)

    SAMPLESHEET_VALIDATION(python_channel,samplesheet_channel)

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
    fasta_channel = Channel
        .fromPath(params.fasta)   // Load reference genome
        .view()


    gtf_channel = Channel
        .fromPath(params.gtf)
        .view()

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
    reads.view()

    reads
        .map { meta, fastq -> fastq }
        .map {fastq1,fastq2 -> tuple([fastq1],[fastq2])}
        .view()
        .set{trimmed_reads}

    // 4. Align reads using HISAT2
    // First, build the HISAT2 index for the reference genome

    hiasat2_build_index = HISAT2_BUILD(fasta_channel)

    // Then, align the reads with the built index
    HISAT2_ALIGN(hiasat2_build_index, reads)
    aligned_sam = HISAT2_ALIGN.out.sam_file

    // Use samtools to sort and index
    SAMTOOLS_SORT_AND_INDEX(aligned_sam)
    sorted_sam = SAMTOOLS_SORT_AND_INDEX.out.sorted

    // Mark duplicates
    PICARD_MARKDUPLICATES(sorted_sam)


    INDEX_FILE(fasta_channel,gtf_channel)
    STAR_ALIGN(trimmed_reads,INDEX_FILE.out.index, gtf_channel)


}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/







