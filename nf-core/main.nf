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
include { FASTQ_ALIGN_HISAT2               } from './subworkflows/fastq_align_hisat2'
include { FASTQ              } from '../../subworkflows/FASTQ/'
include { TRIMMING              } from '../../modules/trimgalore/'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow little_RNASEQ {

    take:
    ch_samplesheet       // channel: path(sample_sheet.csv)
    ch_gtf               // channel: path(genome.gtf)



    main:

    ch_multiqc_files = Channel.empty()
    ch_trim_status = Channel.empty()
    ch_map_status = Channel.empty()
    ch_strand_status = Channel.empty()
    fastqc_zip  = Channel.empty()

    //
    // Create channel from input file provided through params.input
    //
    Channel
        .fromPath("C:\Users\tabat\MasterBioinformatik\Semester3\ComputationalWorkflows\Project\comp_work_project\nf-core\assets\samplesheet.csv")
        .splitCsv(header: true)
        .map {
            meta, fastq_1, fastq_2 ->
                if (!fastq_2) {
                    throw new IllegalArgumentException("Sample ${meta.id} is missing a paired-end file! Only paired-end files are allowed.")
                } else {
                    return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                }
        }
        .groupTuple()
        //.map {
        //    checkSamplesAfterGrouping(it)
        //}
        .set{ ch_fastq }

    //
    // Run RNA-seq FASTQ preprocessing subworkflow
    //

    // The subworkflow only has to do Salmon indexing if it discovers 'auto'
    // samples, and if we haven't already made one elsewhere
    salmon_index_available = params.salmon_index || (!params.skip_pseudo_alignment && params.pseudo_aligner == 'salmon')

    FASTQ (
        ch_fastq,
        ch_fasta,
        ch_transcript_fasta,
        ch_gtf,
        ch_salmon_index,
        ch_sortmerna_index,
        ch_bbsplit_index,
        ch_ribo_db,
        params.skip_bbsplit,
        params.skip_fastqc || params.skip_qc,
        params.skip_trimming,
        params.skip_umi_extract,
        !salmon_index_available,
        !params.sortmerna_index && params.remove_ribo_rna,
        params.trimmer,
        params.min_trimmed_reads,
        params.save_trimmed,
        params.remove_ribo_rna,
        params.with_umi,
        params.umi_discard_read,
        params.stranded_threshold,
        params.unstranded_threshold
    )

    ch_multiqc_files                  = ch_multiqc_files.mix(FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.multiqc_files)
    ch_versions                       = ch_versions.mix(FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.versions)
    ch_strand_inferred_filtered_fastq = FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.reads
    ch_trim_read_count                = FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.trim_read_count

    ch_trim_status = ch_trim_read_count
        .map {
            meta, num_reads ->
                return [ meta.id, num_reads > params.min_trimmed_reads.toFloat() ]
        }


    TRIMMING (ch_fastq)

    trim_unpaired = TRIMGALORE.out.unpaired
    trim_html     = TRIMGALORE.out.html
    trim_zip      = TRIMGALORE.out.zip
    trim_log      = TRIMGALORE.out.log
    ch_versions   = ch_versions.mix(TRIMGALORE.out.versions.first())


     TRIMMING
            .out
            .reads
            .join(trim_log, remainder: true)
            .map {
                meta, reads, trim_log ->
                    if (trim_log) {
                        num_reads = getTrimGaloreReadsAfterFiltering(meta.single_end ? trim_log : trim_log[-1])
                        [ meta, reads, num_reads ]
                    } else {
                        [ meta, reads, min_trimmed_reads.toFloat() + 1 ]
                    }
            }
            .set { ch_num_trimmed_reads }

    ch_num_trimmed_reads
            .filter { meta, reads, num_reads -> num_reads >= min_trimmed_reads.toFloat() }
            .map { meta, reads, num_reads -> [ meta, reads ] }
            .set { trim_reads }

    ch_num_trimmed_reads
            .map { meta, reads, num_reads -> [ meta, num_reads ] }
            .set { trim_read_count }
    }

    ch_filtered_reads = trim_reads
    ch_trim_read_count  =  trim_read_count
    ch_versions = ch_versions
        ch_multiqc_files = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.fastqc_zip
            .mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_zip)
            .mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_log)
            .mix(ch_multiqc_files)

    emit:
    reads = trim_reads // channel: [ val(meta), [ reads ] ]

    fastqc_html        // channel: [ val(meta), [ html ] ]
    fastqc_zip         // channel: [ val(meta), [ zip ] ]

    umi_log            // channel: [ val(meta), [ log ] ]

    trim_unpaired      // channel: [ val(meta), [ reads ] ]
    trim_html          // channel: [ val(meta), [ html ] ]
    trim_zip           // channel: [ val(meta), [ zip ] ]
    trim_log           // channel: [ val(meta), [ txt ] ]
    trim_read_count    // channel: [ val(meta), val(count) ]

    versions = ch_versions.ifEmpty(null) // channel: [ versions.yml ]



    ch_multiqc_files                  = ch_multiqc_files.mix(FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.multiqc_files)
    ch_versions                       = ch_versions.mix(FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.versions)
    ch_strand_inferred_filtered_fastq = FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.reads
    ch_trim_read_count                = FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS.out.trim_read_count

    ch_trim_status = ch_trim_read_count
        .map {
            meta, num_reads ->
                return [ meta.id, num_reads > params.min_trimmed_reads.toFloat() ]
        }









    //
    // Alignment with STAR and gene/transcript quantification with Salmon
    //
    ch_genome_bam          = Channel.empty()
    ch_genome_bam_index    = Channel.empty()
    ch_star_log            = Channel.empty()
    ch_unaligned_sequences = Channel.empty()
    if (!params.skip_alignment && params.aligner == 'star_salmon') {
        // Check if an AWS iGenome has been provided to use the appropriate version of STAR
        def is_aws_igenome = false
        if (params.fasta && params.gtf) {
            if ((file(params.fasta).getName() - '.gz' == 'genome.fa') && (file(params.gtf).getName() - '.gz' == 'genes.gtf')) {
                is_aws_igenome = true
            }
        }

        ALIGN_STAR (
            ch_strand_inferred_filtered_fastq,
            ch_star_index.map { [ [:], it ] },
            ch_gtf.map { [ [:], it ] },
            params.star_ignore_sjdbgtf,
            '',
            params.seq_center ?: '',
            is_aws_igenome,
            ch_fasta.map { [ [:], it ] }
        )
        ch_genome_bam          = ALIGN_STAR.out.bam
        ch_genome_bam_index    = ALIGN_STAR.out.bai
        ch_transcriptome_bam   = ALIGN_STAR.out.bam_transcript
        ch_star_log            = ALIGN_STAR.out.log_final
        ch_unaligned_sequences = ALIGN_STAR.out.fastq
        ch_multiqc_files = ch_multiqc_files.mix(ALIGN_STAR.out.stats.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(ALIGN_STAR.out.flagstat.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(ALIGN_STAR.out.idxstats.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(ch_star_log.collect{it[1]})

        if (params.bam_csi_index) {
            ch_genome_bam_index = ALIGN_STAR.out.csi
        }
        ch_versions = ch_versions.mix(ALIGN_STAR.out.versions)

        
    //
    // Alignment with HISAT2
    //
    if (!params.skip_alignment && params.aligner == 'hisat2') {
        FASTQ_ALIGN_HISAT2 (
            ch_strand_inferred_filtered_fastq,
            ch_hisat2_index.map { [ [:], it ] },
            ch_splicesites.map { [ [:], it ] },
            ch_fasta.map { [ [:], it ] }
        )
        ch_genome_bam          = FASTQ_ALIGN_HISAT2.out.bam
        ch_genome_bam_index    = FASTQ_ALIGN_HISAT2.out.bai
        ch_unaligned_sequences = FASTQ_ALIGN_HISAT2.out.fastq
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_HISAT2.out.stats.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_HISAT2.out.flagstat.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_HISAT2.out.idxstats.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_HISAT2.out.summary.collect{it[1]})

        if (params.bam_csi_index) {
            ch_genome_bam_index = FASTQ_ALIGN_HISAT2.out.csi
        }
        ch_versions = ch_versions.mix(FASTQ_ALIGN_HISAT2.out.versions)

        

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
