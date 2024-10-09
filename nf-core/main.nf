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
include {HIASAT2_BUILD} from './modules/hisat2/align/'

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
    //1. Reading in a samplesheet
    //create a channel to read the samplesheet
    reads_channel = Channel
        // Load the sample sheet specified by the user as a parameter (params.samplesheet)
        .fromPath(params.samplesheet)
        // Split the CSV file into rows.
        //The 'header: true' option tells Nextflow that the first row of the CSV contains column names.
        .splitCsv(header: true)
        // Map each row into a structured format. For each row in the CSV:
        .map { row -> [["sample": row.sample, "strandedness": row.strandedness],[file(row.fastq_1), file(row.fastq_2)]]}
        .view()

    //2. Fastqc
    FASTQC(reads_channel)

    //3. Trimming
    TRIMMING(reads_channel)

    //4.Aligning
    //STAR(reads_channel, params.reference)
    HIASAT2_BUILD(params.reference)
    //HISAT2_ALIGN()





}






/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
