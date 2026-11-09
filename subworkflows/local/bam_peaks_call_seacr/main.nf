//
// Call peaks with SEACR
//

include { SEACR_CALLPEAK           } from '../../../modules/nf-core/seacr/callpeak/main'

workflow BAM_PEAKS_CALL_SEACR {

    take:
    ch_bedgraph_control               // channel: [ val(meta), [ ip_bedgraph ], [ ipcontrol_bedgraph ] ]
    seacr_peak_threshold
    ch_bam                            // channel: [ val(meta), [ ip_bam ], [ control_bam ] ]

    main:

    //
    // Call peaks with SEACR
    //
    SEACR_CALLPEAK (
        ch_bedgraph_control,
        seacr_peak_threshold
    )

    //
    // Filter out samples with 0 SEACR peaks called
    //
    SEACR_CALLPEAK
        .out
        .bed
        .filter {
            meta, peaks ->
                peaks.size() > 0
        }
        .set { ch_seacr_peaks }

    // Create channels: [ meta, ip_bam, peaks ]
    ch_bam
        .map { meta, ip_bam, _control_bam -> [ meta.id, ip_bam ] }
        .combine(
            ch_seacr_peaks.map { meta, peaks -> [ meta.id, meta, peaks ] },
            by: 0
        )
        .map { _id, ip_bam, meta, peaks ->
            [ meta, ip_bam, peaks ]
        }
        .set { ch_bam_peaks }

    emit:
    peaks     = ch_seacr_peaks  // channel: [ val(meta), path(peaks) ]
    bam_peaks = ch_bam_peaks    // channel: [ val(meta), path(ip_bam), path(peaks) ]
}
