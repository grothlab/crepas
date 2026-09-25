
//
// Call peaks with MACS3, annotate with HOMER and perform downstream QC
//
include { MACS3_CALLPEAK           } from '../../../modules/nf-core/macs3/callpeak/main'
include { MACS3_BDGCMP             } from '../../../modules/local/macs3/bdgcmp/main'
include { BEDTOOLS_SLOP                  } from '../../../modules/nf-core/bedtools/slop/main'
include { AWK_FIX_MACS3_BDGCMP       } from '../../../modules/local/awk_fix_macs3_bdgcmp/main'
include { UCSC_BEDCLIP               } from '../../../modules/nf-core/ucsc/bedclip/main'
include { FILE_SORT                  } from '../../../modules/local/file_sort/main'
include { BIGTOOLS_BEDGRAPHTOBIGWIG } from '../../../modules/local/bigtools/bedgraphtobigwig/main'
include { HOMER_ANNOTATEPEAKS      } from '../../../modules/nf-core/homer/annotatepeaks/main'
include { FRIP_SCORE               } from '../../../modules/local/frip_score/main'
include { MULTIQC_CUSTOM_PEAKS     } from '../../../modules/local/multiqc_custom_peaks/main'
include { PLOT_MACS3_QC            } from '../../../modules/local/plot_macs3_qc/main'
include { PLOT_HOMER_ANNOTATEPEAKS } from '../../../modules/local/plot_homer_annotatepeaks/main'

workflow BAM_PEAKS_CALL_QC_ANNOTATE_MACS3_HOMER {
    take:
    ch_bam                            // channel: [ val(meta), [ ip_bam ], [ control_bam ] ]
    ch_fasta                          // channel: [ fasta ]
    ch_gtf                            // channel: [ gtf ]
    ch_chrom_sizes                    // channel: [ bed ]
    macs_gsize                        // integer: value for --macs_gsize parameter
    annotate_peaks_suffix             //  string: suffix for input HOMER annotate peaks files to be trimmed off
    ch_peak_count_header_multiqc      // channel: [ header_file ]
    ch_frip_score_multiqc             // channel: [ header_file ]
    ch_peak_annotation_header_multiqc // channel: [ header_file ]
    is_narrow_peak                    // boolean: true/false
    skip_peak_annotation              // boolean: true/false
    skip_peak_qc                      // boolean: true/false
    skip_bdgcmp                       // boolean: true/false


    main:

    //
    // Call peaks with MACS3
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
        MACS3_BDGCMP (
            ch_bdgs
        )

        BEDTOOLS_SLOP (
            MACS3_BDGCMP.out.bdg,
            ch_chrom_sizes.map{ it -> it[1] }
        )

        AWK_FIX_MACS3_BDGCMP (
            BEDTOOLS_SLOP.out.bed,
        )

        UCSC_BEDCLIP (
            AWK_FIX_MACS3_BDGCMP.out.bed,
            ch_chrom_sizes.map{ it -> it[1] }
        )

        // TODO: maybe a whole module for this is overkill
        FILE_SORT (
            UCSC_BEDCLIP.out.bedgraph,
            'sorted'
        )

        BIGTOOLS_BEDGRAPHTOBIGWIG (
            FILE_SORT.out.sorted,
            ch_chrom_sizes
        )
    }



    // Create channels: [ meta, ip_bam, peaks ]
    ch_bam
        .join(ch_macs3_peaks, by: [0])
        .map {
            meta, ip_bam, control_bam, peaks ->
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
    ch_plot_macs3_qc_txt            = channel.empty()
    ch_plot_macs3_qc_pdf            = channel.empty()
    ch_plot_homer_annotatepeaks_txt = channel.empty()
    ch_plot_homer_annotatepeaks_pdf = channel.empty()
    ch_plot_homer_annotatepeaks_tsv = channel.empty()
    if (!skip_peak_annotation) {
        //
        // Annotate peaks with HOMER
        //
        HOMER_ANNOTATEPEAKS (
            ch_macs3_peaks,
            ch_fasta.map { it -> it[1] },
            ch_gtf.map { it -> it[1] }
        )
        ch_homer_annotatepeaks = HOMER_ANNOTATEPEAKS.out.txt

        if (!skip_peak_qc) {

            // Create channels: [ meta, [ peaks ] ]
            // Where meta = [ id:exp_type, exp_type:exp_type ]
            ch_macs3_peaks
                .map {
                    meta, peaks ->
                        [ meta.exp_type, meta, peaks ]
                }
                .groupTuple()
                .map {
                    exp_type, metas, peaks ->
                        def meta_new = metas[0].clone()
                        meta_new.id = exp_type
                        def sorted_peaks = peaks.sort { it -> it.name }
                        [ meta_new, sorted_peaks ]
                }
                .set { ch_macs3_peaks_grouped }

            //
            // MACS3 QC plots with R
            //
            PLOT_MACS3_QC (
                ch_macs3_peaks_grouped,
                is_narrow_peak
            )
            ch_plot_macs3_qc_txt = PLOT_MACS3_QC.out.txt
            ch_plot_macs3_qc_pdf = PLOT_MACS3_QC.out.pdf

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
    peaks                        = ch_macs3_peaks                   // channel: [ val(meta), [ peaks ] ]
    xls                          = MACS3_CALLPEAK.out.xls           // channel: [ val(meta), [ xls ] ]
    gapped_peaks                 = MACS3_CALLPEAK.out.gapped        // channel: [ val(meta), [ gapped_peak ] ]
    bed                          = MACS3_CALLPEAK.out.bed           // channel: [ val(meta), [ bed ] ]
    pileup_bdg                   = MACS3_CALLPEAK.out.pileup_bdg           // channel: [ val(meta), [ bedgraph ] ]
    lambda_bdg                   = MACS3_CALLPEAK.out.lambda_bdg           // channel: [ val(meta), [ bedgraph ] ]

    frip_txt                     = FRIP_SCORE.out.txt               // channel: [ val(meta), [ txt ] ]

    frip_multiqc                 = MULTIQC_CUSTOM_PEAKS.out.frip    // channel: [ val(meta), [ frip ] ]
    peak_count_multiqc           = MULTIQC_CUSTOM_PEAKS.out.count   // channel: [ val(meta), [ counts ] ]

    homer_annotatepeaks          = ch_homer_annotatepeaks           // channel: [ val(meta), [ txt ] ]

    plot_macs3_qc_txt            = ch_plot_macs3_qc_txt             // channel: [ txt ]
    plot_macs3_qc_pdf            = ch_plot_macs3_qc_pdf             // channel: [ pdf ]

    plot_homer_annotatepeaks_txt = ch_plot_homer_annotatepeaks_txt  // channel: [ txt ]
    plot_homer_annotatepeaks_pdf = ch_plot_homer_annotatepeaks_pdf  // channel: [ pdf ]
    plot_homer_annotatepeaks_tsv = ch_plot_homer_annotatepeaks_tsv  // channel: [ tsv ]
}
