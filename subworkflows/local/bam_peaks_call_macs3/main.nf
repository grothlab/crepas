//
// Call peaks with MACS3, and optionally build signal-over-control bigWig tracks with MACS3 bdgcmp
//
include { MACS3_CALLPEAK            } from '../../../modules/nf-core/macs3/callpeak/main'
include { MACS3_BDGCMP              } from '../../../modules/local/macs3/bdgcmp/main'
include { BEDTOOLS_SLOP             } from '../../../modules/nf-core/bedtools/slop/main'
include { AWK_FIX_MACS3_BDGCMP      } from '../../../modules/local/awk_fix_macs3_bdgcmp/main'
include { UCSC_BEDCLIP              } from '../../../modules/nf-core/ucsc/bedclip/main'
include { FILE_SORT                 } from '../../../modules/local/file_sort/main'
include { BIGTOOLS_BEDGRAPHTOBIGWIG } from '../../../modules/local/bigtools/bedgraphtobigwig/main'

workflow BAM_PEAKS_CALL_MACS3 {
    take:
    ch_bam                            // channel: [ val(meta), [ ip_bam ], [ control_bam ] ]
    ch_chrom_sizes                    // channel: [ bed ]
    macs_gsize                        // integer: value for --macs_gsize parameter
    skip_bdgcmp                       // boolean: true/false


    main:

    //
    // MODULE: Call peaks with MACS3
    //
    MACS3_CALLPEAK (
        ch_bam,
        macs_gsize
    )

    //
    // Filter out samples with 0 MACS3 peaks called
    //
    MACS3_CALLPEAK
        .out
        .peak
        .filter {
            meta, peaks ->
                peaks.size() > 0
        }
        .set { ch_macs3_peaks }

    MACS3_CALLPEAK
        .out
        .pileup_bdg
        .join(MACS3_CALLPEAK.out.lambda_bdg, by: [0])
        .set { ch_bdgs }


    if (!skip_bdgcmp) {
        //
        // MODULE: Fold enrichment of the IP pileup over the control lambda
        //
        MACS3_BDGCMP (
            ch_bdgs
        )

        //
        // MODULE: Clip the fold-enrichment intervals at the chromosome ends
        //
        BEDTOOLS_SLOP (
            MACS3_BDGCMP.out.bdg,
            ch_chrom_sizes.map{ it -> it[1] }
        )

        //
        // MODULE: Remove intervals with an invalid (-1) end coordinate
        //
        AWK_FIX_MACS3_BDGCMP (
            BEDTOOLS_SLOP.out.bed,
        )

        //
        // MODULE: Clip intervals to the chromosome sizes
        //
        UCSC_BEDCLIP (
            AWK_FIX_MACS3_BDGCMP.out.bed,
            ch_chrom_sizes.map{ it -> it[1] }
        )

        //
        // MODULE: Sort the bedGraph
        //
        // TODO: maybe a whole module for this is overkill
        FILE_SORT (
            UCSC_BEDCLIP.out.bedgraph,
            'sorted'
        )

        //
        // MODULE: Convert the fold-enrichment bedGraph to bigWig
        //
        BIGTOOLS_BEDGRAPHTOBIGWIG (
            FILE_SORT.out.sorted,
            ch_chrom_sizes
        )
    }

    // Create channels: [ meta, ip_bam, peaks ]
    ch_bam
        .join(ch_macs3_peaks, by: [0])
        .map {
            meta, ip_bam, _control_bam, peaks ->
                [ meta, ip_bam, peaks ]
        }
        .set { ch_bam_peaks }

    emit:
    peaks        = ch_macs3_peaks                // channel: [ val(meta), path(peaks) ]
    bam_peaks    = ch_bam_peaks                  // channel: [ val(meta), path(ip_bam), path(peaks) ]
    qc_peaks     = ch_macs3_peaks                // channel: [ val(meta), path(peaks) ]
    xls          = MACS3_CALLPEAK.out.xls        // channel: [ val(meta), path(xls) ]
    gapped_peaks = MACS3_CALLPEAK.out.gapped     // channel: [ val(meta), path(gapped_peak) ]
    bed          = MACS3_CALLPEAK.out.bed        // channel: [ val(meta), path(bed) ]
    pileup_bdg   = MACS3_CALLPEAK.out.pileup_bdg // channel: [ val(meta), path(bedgraph) ]
    lambda_bdg   = MACS3_CALLPEAK.out.lambda_bdg // channel: [ val(meta), path(bedgraph) ]
}
