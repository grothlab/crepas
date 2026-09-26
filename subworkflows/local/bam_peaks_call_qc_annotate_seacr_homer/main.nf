//
// Call peaks with SEACR, annotate with HOMER and perform downstream QC
//

include { SEACR_CALLPEAK           } from '../../../modules/nf-core/seacr/callpeak/main'
include { HOMER_ANNOTATEPEAKS      } from '../../../modules/nf-core/homer/annotatepeaks/main'
include { FRIP_SCORE               } from '../../../modules/local/frip_score/main'
include { MULTIQC_CUSTOM_PEAKS     } from '../../../modules/local/multiqc_custom_peaks/main'
include { PLOT_HOMER_ANNOTATEPEAKS } from '../../../modules/local/plot_homer_annotatepeaks/main'

workflow BAM_PEAKS_CALL_QC_ANNOTATE_SEACR_HOMER {

    take:
    ch_bedgraph_control               // channel: [ val(meta), [ ip_bedgraph ], [ ipcontrol_bedgraph ] ]
    seacr_peak_threshold
    ch_bam                             // channel: [ val(meta), [ ip_bam ], [ control_bam ] ]
    ch_fasta                          // channel: [ fasta ]
    ch_gtf                            // channel: [ gtf ]
    annotate_peaks_suffix             //  string: suffix for input HOMER annotate peaks files to be trimmed off
    ch_peak_count_header_multiqc      // channel: [ header_file ]
    ch_frip_score_multiqc             // channel: [ header_file ]
    ch_peak_annotation_header_multiqc // channel: [ header_file ]
    skip_peak_annotation              // boolean: true/false
    skip_peak_qc                      // boolean: true/false

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

    //
    // Calculate FRiP score
    //
    FRIP_SCORE (
        ch_bam_peaks
    )

    // Create channels: [ meta, peaks, frip ]
    ch_bam_peaks
        .join(FRIP_SCORE.out.txt, by: [0])
        .map {
            meta, ip_bam, peaks, frip ->
                [ meta, peaks, frip ]
        }
        .set { ch_bam_peak_frip }

    //
    // FRiP score custom content for MultiQC
    //
    MULTIQC_CUSTOM_PEAKS (
        ch_bam_peak_frip,
        ch_peak_count_header_multiqc,
        ch_frip_score_multiqc
    )

    ch_homer_annotatepeaks          = channel.empty()
    ch_plot_homer_annotatepeaks_txt = channel.empty()
    ch_plot_homer_annotatepeaks_pdf = channel.empty()
    ch_plot_homer_annotatepeaks_tsv = channel.empty()
    if (!skip_peak_annotation) {
        //
        // Annotate peaks with HOMER
        //
        HOMER_ANNOTATEPEAKS (
            ch_seacr_peaks,
            ch_fasta.map { it -> it[1] },
            ch_gtf.map { it -> it[1] }
        )
        ch_homer_annotatepeaks = HOMER_ANNOTATEPEAKS.out.txt

        if (!skip_peak_qc) {

            // Create channels: [ meta, [ anns ] ]
            // Where meta = [ id:exp_type, exp_type:exp_type ]
            HOMER_ANNOTATEPEAKS.out.txt
                .map {
                    meta, anns ->
                        [ meta.exp_type, meta, anns ]
                }
                .groupTuple(by: 0)
                .map {
                    exp_type, metas, anns ->
                        def meta_new = metas[0].clone()
                        meta_new.id = exp_type
                        def sorted_anns = anns.sort { it -> it.name }
                        [ meta_new, sorted_anns ]
                }
                .set { ch_homer_annotatepeaks_grouped }
            //
            // Peak annotation QC plots with R
            //
            PLOT_HOMER_ANNOTATEPEAKS (
                ch_homer_annotatepeaks_grouped,
                ch_peak_annotation_header_multiqc,
                annotate_peaks_suffix
            )
            ch_plot_homer_annotatepeaks_txt = PLOT_HOMER_ANNOTATEPEAKS.out.txt
            ch_plot_homer_annotatepeaks_pdf = PLOT_HOMER_ANNOTATEPEAKS.out.pdf
            ch_plot_homer_annotatepeaks_tsv = PLOT_HOMER_ANNOTATEPEAKS.out.tsv
        }
    }

    emit:
    peaks                        = ch_seacr_peaks                   // channel: [ val(meta), [ peaks ] ]

    frip_txt                     = FRIP_SCORE.out.txt               // channel: [ val(meta), [ txt ] ]

    frip_multiqc                 = MULTIQC_CUSTOM_PEAKS.out.frip    // channel: [ val(meta), [ frip ] ]
    peak_count_multiqc           = MULTIQC_CUSTOM_PEAKS.out.count   // channel: [ val(meta), [ counts ] ]

    homer_annotatepeaks          = ch_homer_annotatepeaks           // channel: [ val(meta), [ txt ] ]

    plot_homer_annotatepeaks_txt = ch_plot_homer_annotatepeaks_txt  // channel: [ txt ]
    plot_homer_annotatepeaks_pdf = ch_plot_homer_annotatepeaks_pdf  // channel: [ pdf ]
    plot_homer_annotatepeaks_tsv = ch_plot_homer_annotatepeaks_tsv  // channel: [ tsv ]
}
