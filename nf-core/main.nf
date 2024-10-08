/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { ALIGN_STAR    } from './subworkflows/local/align_star'

//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//
include { FASTQ_ALIGN_HISAT2               } from '../../subworkflows/nf-core/fastq_align_hisat2'
include { FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS              } from '../../subworkflows/nf-core/fastq_qc_trim_filter_setstrandedness'
