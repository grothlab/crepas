//
// Read QC, UMI extraction and trimming
//

include { FASTQC                            } from '../../../modules/nf-core/fastqc/main'
include { FASTQ_EXTRACT_SPIKEIN_BARCODES } from '../../../modules/local/fastq_extract_spikein_barcodes/main'
include { UMITOOLS_EXTRACT                  } from '../../../modules/nf-core/umitools/extract/main'
include { UMITRANSFER                       } from '../../../modules/local/umitransfer/main'
include { TRIMGALORE as TRIMGALORE          } from '../../../modules/nf-core/trimgalore/main'
include { TRIMGALORE as TRIMGALORE_HARDTRIM } from '../../../modules/nf-core/trimgalore/main'

//
// Function that parses TrimGalore log output file to get total number of reads after trimming
//
def getTrimGaloreReadsAfterFiltering(log_file) {
    def total_reads = 0
    def filtered_reads = 0
    log_file.eachLine { line ->
        def total_reads_matcher = line =~ /([\d\.]+)\ssequences processed in total/
        def filtered_reads_matcher = line =~ /shorter than the length cutoff[^:]+:\s([\d\.]+)/
        if (total_reads_matcher) total_reads = total_reads_matcher[0][1].toFloat()
        if (filtered_reads_matcher) filtered_reads = filtered_reads_matcher[0][1].toFloat()
    }
    return total_reads - filtered_reads
}

workflow FASTQ_FASTQC_UMITOOLS_UMITRANSFER_TRIMGALORE {
    take:
    ch_reads             // channel: [ val(meta), [ reads ] ]
    skip_fastqc       // boolean: true/false
    skip_spikein_barcode_extract // boolean: true/false
    ch_spikein_barcode_table  // channel: [ val(meta), path(spikein_barcode_table.tsv) ]
    skip_trimming     // boolean: true/false
    umi_discard_read  // integer: 0, 1 or 2
    min_trimmed_reads // integer: > 0
    hardtrim3_length  // integer: > 0
    hardtrim5_length  // integer: > 0


    main:
    ch_versions = channel.empty()
    fastqc_html = channel.empty()
    fastqc_zip  = channel.empty()
    if (!skip_fastqc) {
        FASTQC (ch_reads)
        fastqc_html = FASTQC.out.html
        fastqc_zip  = FASTQC.out.zip
    }

    ch_sep_umi_fq       = channel.empty()
    ch_extract_umi    = channel.empty()
    sep_umi_fq_log             = channel.empty()
    no_sep_umi_fq_log          = channel.empty()

    ch_barcode_counts = channel.empty()
    if (!skip_spikein_barcode_extract) {
        FASTQ_EXTRACT_SPIKEIN_BARCODES (
            ch_reads,
            ch_spikein_barcode_table
        )
        ch_barcode_counts = FASTQ_EXTRACT_SPIKEIN_BARCODES.out.counts
        ch_versions = ch_versions.mix(FASTQ_EXTRACT_SPIKEIN_BARCODES.out.versions.first())
    }


    // split ch_umi_reads channel into the ones that have meta.sep_umi_fq and the ones that don't
    ch_reads
        .branch { meta, read -> 
            sep_umi_fq: meta.sep_umi_fq
            extract_umi: !meta.sep_umi_fq && meta.extract_umi
            other: true
        }
        .set { ch_reads }

    UMITRANSFER (ch_reads.sep_umi_fq)
    ch_sep_umi_fq   = UMITRANSFER.out.reads
    sep_umi_fq_log  = UMITRANSFER.out.log
    ch_versions     = ch_versions.mix(UMITRANSFER.out.versions.first())

    UMITOOLS_EXTRACT (ch_reads.extract_umi)
    ch_extract_umi    = UMITOOLS_EXTRACT.out.reads
    no_sep_umi_fq_log   = UMITOOLS_EXTRACT.out.log

    // Discard R1 / R2 if required
    if (umi_discard_read in [1,2]) {
        UMITOOLS_EXTRACT
            .out
            .reads
            .map {
                meta, read ->
                    meta.single_end ? [ meta, read ] : [ meta + ['single_end': true], read[umi_discard_read % 2] ]
            }
            .set { ch_extract_umi }
    }

    ch_reads = ch_reads.other.mix(ch_extract_umi).mix(ch_sep_umi_fq)

    ch_trim_reads   = ch_reads
    ch_htrim_reads  = ch_reads
    trim_unpaired   = channel.empty()
    trim_html       = channel.empty()
    trim_zip        = channel.empty()
    trim_log        = channel.empty()
    trim_read_count = channel.empty()
    if (!skip_trimming) {
        TRIMGALORE (ch_reads)
        trim_unpaired = TRIMGALORE.out.unpaired
        trim_html     = TRIMGALORE.out.html
        trim_zip      = TRIMGALORE.out.zip
        trim_log      = TRIMGALORE.out.log

        //
        // Filter FastQ files based on minimum trimmed read count after adapter trimming
        //
        TRIMGALORE
            .out
            .reads
            .join(trim_log, remainder: true)
            .map {
                meta, read, log ->
                    if (log) {
                        def num_reads = getTrimGaloreReadsAfterFiltering(meta.single_end ? log : log[-1])
                        [ meta, read, num_reads ]
                    } else {
                        [ meta, read, min_trimmed_reads.toFloat() + 1 ]
                    }
            }
            .set { ch_num_trimmed_reads }

        ch_num_trimmed_reads
            .filter { meta, read, num_reads -> num_reads >= min_trimmed_reads.toFloat() }
            .map { meta, read, num_reads -> [ meta, read ] }
            .set { ch_trim_reads }

        ch_num_trimmed_reads
            .map { meta, read, num_reads -> [ meta, num_reads ] }
            .set { trim_read_count }

        htrim_unpaired   = channel.empty()
        htrim_html       = channel.empty()
        htrim_zip        = channel.empty()
        htrim_log        = channel.empty()
        ch_htrim_reads   = ch_trim_reads
        if (hardtrim3_length || hardtrim5_length) {
            TRIMGALORE_HARDTRIM (ch_trim_reads)
            ch_htrim_reads  = TRIMGALORE_HARDTRIM.out.reads
            htrim_unpaired  = TRIMGALORE_HARDTRIM.out.unpaired
            htrim_html      = TRIMGALORE_HARDTRIM.out.html
            htrim_zip       = TRIMGALORE_HARDTRIM.out.zip
            htrim_log       = TRIMGALORE_HARDTRIM.out.log
        }

    }

    emit:
    reads          = ch_htrim_reads // channel: [ val(meta), [ reads ] ]

    barcode_counts = ch_barcode_counts // channel: [ val(meta), path(*.tsv) ]

    fastqc_html        // channel: [ val(meta), [ html ] ]
    fastqc_zip         // channel: [ val(meta), [ zip ] ]

    sep_umi_fq_log            // channel: [ val(meta), [ log ] ]
    no_sep_umi_fq_log         // channel: [ val(meta), [ log ] ]

    trim_unpaired      // channel: [ val(meta), [ reads ] ]
    trim_html          // channel: [ val(meta), [ html ] ]
    trim_zip           // channel: [ val(meta), [ zip ] ]
    trim_log           // channel: [ val(meta), [ txt ] ]
    trim_read_count    // channel: [ val(meta), val(count) ]

    htrim_unpaired    // channel: [ val(meta), [ reads ] ]
    htrim_html        // channel: [ val(meta), [ html ] ]
    htrim_zip         // channel: [ val(meta), [ zip ] ]
    htrim_log         // channel: [ val(meta), [ txt ] ]

    versions = ch_versions // channel: [ versions.yml ]
}
