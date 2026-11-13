//
// Call peaks with LanceOtron, on 1 bp coverage of the BAM files
//

include { DEEPTOOLS_BAMCOVERAGE as LANCEOTRON_BAMCOVERAGE } from '../../../modules/nf-core/deeptools/bamcoverage/main'
include { LANCEOTRON_CALLPEAKS                            } from '../../../modules/local/lanceotron/callpeaks/main'
include { LANCEOTRON_CALLPEAKSINPUT                       } from '../../../modules/local/lanceotron/callpeaksinput/main'
include { LANCEOTRON_TO_PEAK                              } from '../../../modules/local/lanceotron_to_peak/main'

workflow BAM_PEAKS_CALL_LANCEOTRON {

    take:
    ch_bam_bai             // channel: [ val(meta), path(bam), path(bai) ], IP and input control BAM files

    main:

    // LanceOtron scores the shape of the coverage, so it gets its own 1 bp, read-extended coverage track
    // (the coverage tracks are binned at --coverage_bin_size). One track per BAM file: input controls shared
    // by several antibodies are one file.
    ch_bam_bai
        .unique { _meta, bam, _bai -> bam }
        .set { ch_bam_bai_unique }

    //
    // MODULE: 1 bp coverage of each BAM file
    //
    LANCEOTRON_BAMCOVERAGE (
        ch_bam_bai_unique,
        [],
        [],
        channel.value([[:], []])
    )

    // Create channel: [ meta, bigwig ], for every sample entry of a BAM file
    ch_bam_bai_unique
        .map { meta, bam, _bai -> [ meta, bam ] }
        .join(LANCEOTRON_BAMCOVERAGE.out.bigwig, by: 0)
        .map { _meta, bam, bigwig -> [ bam, bigwig ] }
        .set { ch_bam_bigwig }

    ch_bam_bai
        .map { meta, bam, _bai -> [ bam, meta ] }
        .combine(ch_bam_bigwig, by: 0)
        .map { _bam, meta, bigwig -> [ meta, bigwig ] }
        .set { ch_bigwig }

    // Split the coverage tracks into IPs with an input control, IPs without one, and input controls
    // (which are not called themselves)
    ch_bigwig
        .branch { meta, bigwig ->
            ips_with_ipcontrol: meta.input_control
                return [ meta.input_control, meta.antibody, meta, bigwig ]
            ips_wo_ipcontrol: !meta.input_control && !meta.is_input_control
                return [ meta, bigwig ]
            ipcontrols: !meta.input_control && meta.is_input_control
                return [ meta.id, meta.input_control_of_antibody, meta, bigwig ]
        }
        .set { ch_bigwig_by_type }

    // Create channel: [ meta, ip_bigwig, ipcontrol_bigwig ], each IP with its input control
    ch_bigwig_by_type
        .ips_with_ipcontrol
        .combine(ch_bigwig_by_type.ipcontrols, by: [ 0, 1 ])
        .map { _ipcontrol_id, _antibody, ip_meta, ip_bigwig, _ipcontrol_meta, ipcontrol_bigwig ->
            [ ip_meta, ip_bigwig, ipcontrol_bigwig ]
        }
        .set { ch_bigwig_ip_ipcontrol }

    //
    // MODULE: Call peaks with LanceOtron, for IPs with an input control (adds a p-value against the control)
    //
    LANCEOTRON_CALLPEAKSINPUT (
        ch_bigwig_ip_ipcontrol
    )

    //
    // MODULE: Call peaks with LanceOtron, for IPs without an input control
    //
    LANCEOTRON_CALLPEAKS (
        ch_bigwig_by_type.ips_wo_ipcontrol
    )

    //
    // MODULE: Keep the candidates above the score cutoff, as narrowPeak, LanceOtron table and QC scores
    //
    LANCEOTRON_TO_PEAK (
        LANCEOTRON_CALLPEAKSINPUT.out.candidates.mix(LANCEOTRON_CALLPEAKS.out.candidates)
    )

    // Drop empty peak files
    LANCEOTRON_TO_PEAK
        .out
        .narrowpeak
        .filter { _meta, peaks -> peaks.size() > 0 }
        .set { ch_lanceotron_peaks }

    // Create channel: [ meta, scores ], for the peak QC plots, for the samples with peaks
    LANCEOTRON_TO_PEAK
        .out
        .scores
        .join(ch_lanceotron_peaks, by: 0)
        .map { meta, scores, _peaks -> [ meta, scores ] }
        .set { ch_lanceotron_scores }

    // Create channel: [ meta, ip_bam, peaks ], for FRiP scores
    ch_bam_bai
        .filter { meta, _bam, _bai -> !meta.is_input_control }
        .map { meta, ip_bam, _bai -> [ meta.id, ip_bam ] }
        .combine(
            ch_lanceotron_peaks.map { meta, peaks -> [ meta.id, meta, peaks ] },
            by: 0
        )
        .map { _id, ip_bam, meta, peaks ->
            [ meta, ip_bam, peaks ]
        }
        .set { ch_bam_peaks }

    emit:
    peaks     = ch_lanceotron_peaks                 // channel: [ val(meta), path(narrowPeak) ]
    bam_peaks = ch_bam_peaks                        // channel: [ val(meta), path(ip_bam), path(narrowPeak) ]
    qc_peaks  = ch_lanceotron_scores                // channel: [ val(meta), path(scores) ]
    table     = LANCEOTRON_TO_PEAK.out.tsv          // channel: [ val(meta), path(tsv) ]
    bigwig    = LANCEOTRON_BAMCOVERAGE.out.bigwig   // channel: [ val(meta), path(bigwig) ]
}
