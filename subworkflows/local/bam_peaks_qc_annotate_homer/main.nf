//
// Peak QC and annotation shared by every peak caller: FRiP scores, peak counts, HOMER annotation and QC plots
//

include { FRIP_SCORE               } from '../../../modules/local/frip_score/main'
include { MULTIQC_CUSTOM_PEAKS     } from '../../../modules/local/multiqc_custom_peaks/main'
include { HOMER_ANNOTATEPEAKS      } from '../../../modules/nf-core/homer/annotatepeaks/main'
include { PLOT_MACS3_QC            } from '../../../modules/local/plot_macs3_qc/main'
include { PLOT_HOMER_ANNOTATEPEAKS } from '../../../modules/local/plot_homer_annotatepeaks/main'

workflow BAM_PEAKS_QC_ANNOTATE_HOMER {

    take:
    ch_bam_peaks          // channel: [ val(meta), path(ip_bam), path(peaks) ], one per IP BAM, for FRiP scores and peak counts
    ch_peaks              // channel: [ val(meta), path(peaks) ], one per peak file, for HOMER annotation
    ch_qc_peaks           // channel: [ val(meta), path(peaks) ], MACS3-like peak files (fold enrichment, p- and q-value columns) for the peak QC plots
    ch_peak_count_frip_headers  // channel: [ val(peak_caller), val(norm_factor_type), path(peak_count_header), path(frip_score_header) ], norm_factor_type is '' for callers not run per normalisation
    ch_peak_annotation_headers  // channel: [ val(peak_caller), path(peak_annotation_header), val(annotate_peaks_suffix) ]
    ch_fasta              // channel: [ val(meta), path(fasta) ]
    ch_gtf                // channel: [ val(meta), path(gtf) ]
    skip_peak_annotation  // boolean: true/false
    skip_peak_qc          // boolean: true/false

    main:

    // Every meta carries meta.peak_caller, which selects the MultiQC headers and the output names/paths.
    // Peaks of callers run on each normalisation of a sample (SEACR) also carry meta.norm_factor_type,
    // and get the peak count and FRiP score headers of that normalisation.

    //
    // Calculate FRiP score
    //
    FRIP_SCORE (
        ch_bam_peaks
    )

    // Create channel: [ meta, peaks, frip, peak_count_header, frip_score_header ]
    ch_bam_peaks
        .join(FRIP_SCORE.out.txt, by: 0)
        .map { meta, _ip_bam, peaks, frip -> [ meta.peak_caller, meta.norm_factor_type ?: '', meta, peaks, frip ] }
        .combine(ch_peak_count_frip_headers, by: [ 0, 1 ])
        .map { _peak_caller, _norm_factor_type, meta, peaks, frip, peak_count_header, frip_score_header ->
            [ meta, peaks, frip, peak_count_header, frip_score_header ]
        }
        .set { ch_peaks_frip_headers }

    //
    // FRiP score and peak count custom content for MultiQC
    //
    MULTIQC_CUSTOM_PEAKS (
        ch_peaks_frip_headers
    )

    ch_homer_annotatepeaks          = channel.empty()
    ch_plot_peak_qc_txt             = channel.empty()
    ch_plot_peak_qc_pdf             = channel.empty()
    ch_plot_homer_annotatepeaks_txt = channel.empty()
    ch_plot_homer_annotatepeaks_pdf = channel.empty()
    ch_plot_homer_annotatepeaks_tsv = channel.empty()
    if (!skip_peak_annotation) {
        //
        // Annotate peaks with HOMER
        //
        HOMER_ANNOTATEPEAKS (
            ch_peaks,
            ch_fasta.map { it -> it[1] },
            ch_gtf.map { it -> it[1] }
        )
        ch_homer_annotatepeaks = HOMER_ANNOTATEPEAKS.out.txt

        if (!skip_peak_qc) {

            // Create channel: [ meta, [ peaks ] ], one per peak caller and experiment type
            // meta is that of the first sample, with id = exp_type
            ch_qc_peaks
                .map { meta, peaks -> [ [ meta.peak_caller, meta.exp_type ], meta, peaks ] }
                .groupTuple(by: 0)
                .map { key, metas, peaks ->
                    [ metas[0] + [ id: key[1] ], peaks.sort { it -> it.name } ]
                }
                .set { ch_qc_peaks_grouped }

            //
            // Peak QC plots (fold enrichment, p-value, q-value, length) with R
            //
            PLOT_MACS3_QC (
                ch_qc_peaks_grouped
            )
            ch_plot_peak_qc_txt = PLOT_MACS3_QC.out.txt
            ch_plot_peak_qc_pdf = PLOT_MACS3_QC.out.pdf

            // Create channel: [ meta, [ anns ], peak_annotation_header, annotate_peaks_suffix ], one per peak caller and experiment type
            HOMER_ANNOTATEPEAKS.out.txt
                .map { meta, anns -> [ [ meta.peak_caller, meta.exp_type ], meta, anns ] }
                .groupTuple(by: 0)
                .map { key, metas, anns ->
                    [ key[0], metas[0] + [ id: key[1] ], anns.sort { it -> it.name } ]
                }
                .combine(ch_peak_annotation_headers, by: 0)
                .map { _peak_caller, meta, anns, peak_annotation_header, suffix ->
                    [ meta, anns, peak_annotation_header, suffix ]
                }
                .set { ch_homer_annotatepeaks_grouped }

            //
            // Peak annotation QC plots with R
            //
            PLOT_HOMER_ANNOTATEPEAKS (
                ch_homer_annotatepeaks_grouped
            )
            ch_plot_homer_annotatepeaks_txt = PLOT_HOMER_ANNOTATEPEAKS.out.txt
            ch_plot_homer_annotatepeaks_pdf = PLOT_HOMER_ANNOTATEPEAKS.out.pdf
            ch_plot_homer_annotatepeaks_tsv = PLOT_HOMER_ANNOTATEPEAKS.out.tsv
        }
    }

    emit:
    frip_txt                     = FRIP_SCORE.out.txt               // channel: [ val(meta), path(txt) ]
    frip_multiqc                 = MULTIQC_CUSTOM_PEAKS.out.frip    // channel: [ val(meta), path(tsv) ]
    peak_count_multiqc           = MULTIQC_CUSTOM_PEAKS.out.count   // channel: [ val(meta), path(tsv) ]
    homer_annotatepeaks          = ch_homer_annotatepeaks           // channel: [ val(meta), path(txt) ]
    plot_peak_qc_txt             = ch_plot_peak_qc_txt              // channel: [ val(meta), path(txt) ]
    plot_peak_qc_pdf             = ch_plot_peak_qc_pdf              // channel: [ val(meta), path(pdf) ]
    plot_homer_annotatepeaks_txt = ch_plot_homer_annotatepeaks_txt  // channel: [ val(meta), path(txt) ]
    plot_homer_annotatepeaks_pdf = ch_plot_homer_annotatepeaks_pdf  // channel: [ val(meta), path(pdf) ]
    plot_homer_annotatepeaks_tsv = ch_plot_homer_annotatepeaks_tsv  // channel: [ val(meta), path(tsv) ]
}
