/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULES
//
include { FASTQC } from './modules/FASTQ/'
include {TRIMMING} from './modules/trimgalore/'
include {STAR} from './modules/STAR/'
include {HISAT2_BUILD} from './modules/hisat2/align/'
include {HISAT2_ALIGN} from './modules/hisat2/align/'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
params.samplesheet = params.samplesheet ?: "./modules/trimgalore/samplesheet_new.csv"
params.reference =  params.reference ?: "./modules/hisat2/align/genome.fa"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow {
    // 1. Reading in a samplesheet
    // Create a channel to read the samplesheet
    reads_channel = Channel
        // Load the sample sheet specified by the user as a parameter (params.samplesheet)
        .fromPath(params.samplesheet)
        // Split the CSV file into rows.
        // The 'header: true' option tells Nextflow that the first row of the CSV contains column names.
        .splitCsv(header: true)
        // Map each row into a structured format:
        // Each row contains the sample name, strandedness, and two FASTQ file paths
        .map { row ->
            [
                ["sample": row.sample, "strandedness": row.strandedness],
                [file(row.fastq_1), file(row.fastq_2)]
            ]
        }
        .view()

        //.view()  // View the contents of the channel for debugging

    // Channel for reference genome files
    reference_channel = Channel
        .fromPath(params.reference)

    // Tuple channel for FASTQ files
    tuple_channel = Channel
        // Load the sample sheet specified by the user as a parameter (params.samplesheet)
        .fromPath(params.samplesheet)
        // Split the CSV file into rows.
        // The 'header: true' option tells Nextflow that the first row of the CSV contains column names.
        .splitCsv(header: true)
        // Map each row into a structured format:
        // Each row contains the sample name, strandedness, and two FASTQ file paths
        .map { row ->
            [
                ["sample": row.sample, "strandedness": row.strandedness],
                [file(row.fastq_1)],
                [file(row.fastq_2)]
            ]
        }
        .view()




    // 2. FastQC: Quality control for the FASTQ files
    FASTQC(reads_channel)

    // 3. Trimming: Trimming the FASTQ files for quality
    TRIMMING(reads_channel)

    // 4. Aligning
    // Build the HISAT2 index for the reference genome
    HISAT2_BUILD(reference_channel)

    HISAT2_ALIGN(tuple_channel)

  // Access the first FASTQ files

}





/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
