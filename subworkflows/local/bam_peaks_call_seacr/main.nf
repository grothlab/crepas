//
// Call peaks with SEACR, on the coverage tracks and on 1 bp coverage of the BAM files
//

include { BEDTOOLS_GENOMECOV } from '../../../modules/nf-core/bedtools/genomecov/main'
include { BEDGRAPH_NORMALIZE } from '../../../modules/local/bedgraph_normalize/main'
include { SEACR_CALLPEAK     } from '../../../modules/nf-core/seacr/callpeak/main'

workflow BAM_PEAKS_CALL_SEACR {

    take:
    ch_bam                 // channel: [ val(meta), path(bam) ], IP and input control BAM files
    ch_bedgraph            // channel: [ val(meta), path(bedgraph) ], coverage tracks, one per BAM file and normalisation, with meta.norm_factor_type and meta.norm_factor_val
    coverage_bin_size      // integer: bin size of the coverage tracks
    seacr_peak_threshold   // float: fraction of top peaks to keep when there is no input control

    main:

    // SEACR runs at two bin sizes (meta.seacr_bin_size): the coverage tracks' bin size and 1 bp.
    // A BAM file and its coverage tracks are matched by [ id, genome, input_control_of_antibody ].

    // 1 bp coverage: once per BAM file that has coverage tracks, unless the tracks are 1 bp already
    ch_bam
        .filter { _meta, _bam -> coverage_bin_size != 1 }
        .map { meta, bam -> [ [ meta.id, meta.genome, meta.input_control_of_antibody ?: '' ], meta, bam ] }
        .combine(
            ch_bedgraph.map { meta, _bdg -> [ meta.id, meta.genome, meta.input_control_of_antibody ?: '' ] }.unique().map { key -> [ key ] },
            by: 0
        )
        .map { _key, meta, bam -> [ meta, bam ] }
        .unique { _meta, bam -> bam }          // input controls shared by several antibodies are one file
        .map { meta, bam -> [ meta, bam, 1 ] } // unscaled; scaled per normalisation below
        .set { ch_bam_unique }

    //
    // MODULE: 1 bp coverage of each BAM file
    //
    BEDTOOLS_GENOMECOV (
        ch_bam_unique,
        [],
        'bedGraph',
        false
    )

    // Create channel: [ bam, bedgraph_1bp ]
    ch_bam_unique
        .map { meta, bam, _scale -> [ meta, bam ] }
        .join(BEDTOOLS_GENOMECOV.out.genomecov, by: 0)
        .map { _meta, bam, bdg -> [ bam, bdg ] }
        .set { ch_bedgraph_1bp }

    // Create channel: [ meta, bedgraph ], each coverage track at its own bin size and at 1 bp
    ch_bedgraph
        .map { meta, _bdg -> [ [ meta.id, meta.genome, meta.input_control_of_antibody ?: '' ], meta ] }
        .combine(
            ch_bam.map { meta, bam -> [ [ meta.id, meta.genome, meta.input_control_of_antibody ?: '' ], bam ] },
            by: 0
        )
        .map { _key, meta, bam -> [ bam, meta ] }
        .combine(ch_bedgraph_1bp, by: 0)
        .map { _bam, meta, bdg -> [ meta + [ seacr_bin_size: 1, peak_qc_variant: "1bp_${meta.norm_factor_type}" ], bdg ] }
        .mix(
            ch_bedgraph.map { meta, bdg -> [ meta + [ seacr_bin_size: coverage_bin_size, peak_qc_variant: meta.norm_factor_type ], bdg ] }
        )
        .set { ch_seacr_bedgraph }

    //
    // MODULE: Scale the 1 bp coverage per normalisation, and add a header line to all bedGraphs (SEACR skips the first line)
    //
    BEDGRAPH_NORMALIZE (
        ch_seacr_bedgraph,
        'bedGraph'
    )

    // Create channel: [ meta, ip_bedgraph, ipcontrol_bedgraph ], each IP with the input control of the same bin size
    // and normalisation; IPs without a control are called with a numeric threshold
    BEDGRAPH_NORMALIZE
        .out
        .normalized
        .branch { meta, bdg ->
            ips_with_ipcontrol: meta.input_control
                return [ meta.input_control, meta.antibody, meta.norm_factor_type, meta.seacr_bin_size, meta, bdg ]
            ips_wo_ipcontrol: !meta.input_control && !meta.is_input_control
                return [ meta, bdg, [] ]
            ipcontrols: !meta.input_control && meta.is_input_control
                return [ meta.id, meta.input_control_of_antibody, meta.norm_factor_type, meta.seacr_bin_size, meta, bdg ]
        }
        .set { ch_bedgraph_by_type }

    ch_bedgraph_by_type
        .ips_with_ipcontrol
        .combine(ch_bedgraph_by_type.ipcontrols, by: [ 0, 1, 2, 3 ])
        .map { _ipcontrol_id, _antibody, _norm_factor_type, _seacr_bin_size, ip_meta, ip_bdg, _ipcontrol_meta, ipcontrol_bdg ->
            [ ip_meta, ip_bdg, ipcontrol_bdg ]
        }
        .mix(ch_bedgraph_by_type.ips_wo_ipcontrol)
        .set { ch_bedgraph_control }

    //
    // MODULE: Call peaks with SEACR
    //
    SEACR_CALLPEAK (
        ch_bedgraph_control,
        seacr_peak_threshold
    )

    // Drop empty peak files
    SEACR_CALLPEAK
        .out
        .bed
        .filter { _meta, peaks -> peaks.size() > 0 }
        .set { ch_seacr_peaks }

    // Create channel: [ meta, ip_bam, peaks ], for FRiP scores
    ch_bam
        .filter { meta, _bam -> !meta.is_input_control }
        .map { meta, ip_bam -> [ meta.id, ip_bam ] }
        .combine(
            ch_seacr_peaks.map { meta, peaks -> [ meta.id, meta, peaks ] },
            by: 0
        )
        .map { _id, ip_bam, meta, peaks ->
            [ meta, ip_bam, peaks ]
        }
        .set { ch_bam_peaks }

    emit:
    peaks     = ch_seacr_peaks                       // channel: [ val(meta), path(peaks) ]
    bam_peaks = ch_bam_peaks                         // channel: [ val(meta), path(ip_bam), path(peaks) ]
    bedgraph  = BEDGRAPH_NORMALIZE.out.normalized    // channel: [ val(meta), path(bedgraph) ], the bedGraphs given to SEACR
}
