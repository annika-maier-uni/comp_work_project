/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simple RNA Sequencing Pipeline
# Author: Weronika Jaśkowiak, Maike Nägele, Tabea Attig, Annika Maier
# Date: 11.10.2024
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Load the required modules for FastQC, trimming, and alignment with HISAT2 and STAR
include { FASTQC } from './modules/FASTQ/'
include { TRIMMING } from './modules/trimgalore/'
include { HISAT2_BUILD } from './modules/hisat2/align/'
include { HISAT2_ALIGN } from './modules/hisat2/align/'
include { PICARD_MARKDUPLICATES } from './modules/picard_mark_duplicates/'
include { SAMTOOLS_SORT_AND_INDEX } from './modules/samtools/'
include { INDEX_FILE } from './modules/STAR/'
include { STAR_ALIGN } from './modules/STAR/'
include { SAMPLESHEET_VALIDATION } from './modules/samplesheet_validation'
// include { VALIDATION_SUCCESS } from './modules/samplesheet_validation'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Set default parameters for input files and tools
params.samplesheet = params.samplesheet ?: './data/samplesheets/samplesheet.csv'   // Path to samplesheet file (CSV format) containing sample metadata
params.fasta  = params.fasta  ?: "./data/genome.fa"                   // Path to reference genome file in FASTA format
params.gtf  = params.gtf  ?: "./data/genes.gtf"                       // Path to gene annotation file in GTF format
params.python_file = './bin/validation_samplesheet.py'                // Path to the Python script for validating the samplesheet
params.align = params.align ?: 'HISAT2'                               // Default alignment tool: HISAT2 (can be overridden)
params.outdir = params.outdir  ?: 'results'                           // Default output directory where results will be stored
params.adapter = params.adapter ?: ''                                 // Default adapter sequence for trimming (empty if not provided)
params.length = params.length ?: ''                                   // Default minimum read length for trimming (empty if not provided)
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    // CHANNELS: Define channels for different input files and processed data

    // 1. Channel for the samplesheet file
    samplesheet_channel = Channel
        .fromPath(params.samplesheet)

    // 2. Channel for the Python validation script file
    python_channel = Channel
        .fromPath(params.python_file)



    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        EXECUTION: Validate samplesheet, perform quality control, trimming, alignment, sorting and marking duplicates
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // 1. Samplesheet validation using the provided Python script
    samplesheet_channel = SAMPLESHEET_VALIDATION(python_channel, samplesheet_channel)

    // 3. Channel for reading data from the samplesheet file
    samplesheet_channel
        .splitCsv(header: true)      // Split CSV file using header for column names
        .map { row ->                // Map each row to structured format
            [
                ["sample": row.sample, "strandedness": row.strandedness],
                [file(row.fastq_1), file(row.fastq_2)] // Paired FASTQ file paths
            ]
        }
        .set {reads_channel}

    // 4. Channel for the reference genome (FASTA) file
    fasta_channel = Channel
        .fromPath(params.fasta)

    // 5. Channel for the GTF (Gene Transfer Format) file
    gtf_channel = Channel
        .fromPath(params.gtf)

    // 4. Perform quality control using FastQC
    FASTQC(reads_channel)

    // 5. Trim the FASTQ files for quality using Trim Galore
    TRIMMING(reads_channel)
    reads = TRIMMING.out.reads  // Capture trimmed reads for alignment

    // 6. Align reads using HISAT2 if chosen
    if (params.align == 'HISAT2') {
        // 6a. Build HISAT2 index from the reference genome
        hisat2_build_index = HISAT2_BUILD(fasta_channel)

        // 6b. Align the reads using HISAT2
        HISAT2_ALIGN(hisat2_build_index, reads)
        aligned_sam = HISAT2_ALIGN.out.sam

        // 6c. Sort and index the aligned reads using Samtools
        SAMTOOLS_SORT_AND_INDEX(aligned_sam)
        sorted_sam = SAMTOOLS_SORT_AND_INDEX.out.sam

        // 7. Mark duplicates using Picard
        PICARD_MARKDUPLICATES(sorted_sam)
    }

    // 8. Alternatively, align reads using STAR
    else {
        // Prepare reads channel for STAR
        reads
            .map { meta, fastq -> fastq }
            .map { fastq1, fastq2 -> tuple([fastq1], [fastq2]) }
            .set { trimmed_reads }

        // 8a. Build STAR index
        INDEX_FILE(fasta_channel, gtf_channel)

        // 8b. Align reads using STAR with the generated index
        STAR_ALIGN(trimmed_reads, INDEX_FILE.out.index, gtf_channel)
        sorted_bam = STAR_ALIGN.out.bam
    }
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/