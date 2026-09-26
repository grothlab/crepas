//
// Call peaks with MACE, annotate with HOMER and perform downstream QC
//
include { MACE_PREPROCESSOR        } from '../../../modules/local/mace/preprocessor/main'
include { UCSC_WIGTOBIGWIG         } from '../../../modules/nf-core/ucsc/wigtobigwig/main'
include { MACE_MACE                } from '../../../modules/local/mace/mace/main'
include { HOMER_ANNOTATEPEAKS      } from '../../../modules/nf-core/homer/annotatepeaks/main'
include { FRIP_SCORE               } from '../../../modules/local/frip_score/main'
include { MULTIQC_CUSTOM_PEAKS     } from '../../../modules/local/multiqc_custom_peaks/main'
// include { PLOT_MACS3_QC as PLOT_MACE_QC } from '../../../modules/local/plot_macs3_qc/main'
include { PLOT_HOMER_ANNOTATEPEAKS } from '../../../modules/local/plot_homer_annotatepeaks/main'

workflow BAM_PEAKS_CALL_QC_ANNOTATE_MACE_HOMER {
    take:
    ch_bam_bai                        // channel: [ val(meta), bam, bai ]
    ch_fasta                          // channel: [ fasta  ]
    ch_gtf                            // channel: [ gtf ]
    ch_chrom_sizes                    // channel: [ val(meta), path(chrom_sizes) ]
    annotate_peaks_suffix             //  string: suffix for input HOMER annotate peaks files to be trimmed off
    ch_peak_count_header_multiqc      // channel: [ header_file ]
    ch_frip_score_multiqc             // channel: [ header_file ]
    ch_peak_annotation_header_multiqc // channel: [ header_file ]
    skip_peak_annotation              // boolean: true/false
    skip_peak_qc                      // boolean: true/false

    main:

    // Create channel: [ meta, [bams_merged_reps] ]
    ch_bam_bai
        .map { meta, bam, bai ->
            def meta_clone = meta.clone()
            def antibody = meta_clone.input_control_of_antibody ?: meta_clone.antibody
            meta_clone.id = meta_clone.id - ~/_bRep_.*$/
            meta_clone.input_control = meta_clone.input_control - ~/_bRep_.*$/
            [meta_clone.id, antibody, meta_clone, bam, bai]
        }
        .groupTuple(by: [0, 1])
        .map { id, antibody, metas, bams, bais ->
            [metas[0], bams.flatten(), bais.flatten()]
        }
        .set { ch_bam_merged_reps }

    // TODO: Print to file for debuggin
    ch_bam_merged_reps
        .map { meta, bams, bais ->
            "${meta}\t${bams}\t${bais}"
        }
        .collectFile(name: 'ch_bam_merged_reps.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_PEAKS_CALL_QC_ANNOTATE_MACE_HOMER")

    //
    // MODULE: Preprocess BAM files for MACE peak caller
    //
    MACE_PREPROCESSOR (
        ch_bam_merged_reps,
        ch_chrom_sizes
    )

    // Add strand to the meta information
    MACE_PREPROCESSOR.out.forward_wig.map { meta, forward_wig ->
        def meta_clone = meta.clone()
        meta_clone.strand = 'forward'
        [meta_clone, forward_wig]
    }.set { ch_fwd_wig }

    MACE_PREPROCESSOR.out.reverse_wig.map { meta, reverse_wig ->
        def meta_clone = meta.clone()
        meta_clone.strand = 'reverse'
        [meta_clone, reverse_wig]
    }.set { ch_rwd_wig }

    // Merge forward and reverse strands into one channel
    ch_wig = ch_fwd_wig.mix(ch_rwd_wig)

    UCSC_WIGTOBIGWIG (
        ch_wig,
        ch_chrom_sizes.map { it -> it[1] }
    )

    // Split channel by strand and create channels: [ meta, bw ]
    UCSC_WIGTOBIGWIG.out.bw.map { meta, bw ->
        def meta_clone = meta.clone()
        meta_clone.remove('strand')
        [meta_clone, meta, bw]
    }.branch { meta_clone, meta, bw ->
        forward: meta.strand == 'forward'
        return [meta_clone, bw]
        reverse: meta.strand == 'reverse'
        return [meta_clone, bw]
    }.set { ch_bw_by_strand }

    // Create channel: [ meta, forward_bw, reverse_bw ]
    ch_bw_by_strand.forward.combine(ch_bw_by_strand.reverse, by: 0).set { ch_bw }

    //
    // MODULE: Border detection and border pairing with MACE
    //
    MACE_MACE (
        ch_bw,
        ch_chrom_sizes
    )

    //
    // Filter out samples with 0 MACE peaks called
    //
    MACE_MACE.out.border_pair_peak.filter { meta, peaks ->
        peaks.size() > 0
    }.set { ch_mace_peaks }

    // Create channels: [ meta, bam, peaks ]
    ch_bam_bai
        .join(ch_mace_peaks, by: 0)
        .map { meta, bam, bai, peaks ->
            [meta, bam, bai, peaks]
        }
        .transpose()
        .map { meta, bam, bai, peaks ->
            def meta_clone = meta.clone()
            meta_clone.id = bam.getSimpleName()
            [meta_clone, bam, peaks]
        }
        .set { ch_bam_peak }

    //
    // Calculate FRiP score
    //
    FRIP_SCORE (
        ch_bam_peak
    )

    // Create channels: [ meta, peaks, frip ]
    ch_bam_peak
        .join(FRIP_SCORE.out.txt, by: 0)
        .map { meta, bam, peaks, frip ->
            [meta, peaks, frip]
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

    ch_homer_annotatepeaks = channel.empty()
    // ch_plot_mace_qc_txt            = channel.empty()
    // ch_plot_mace_qc_pdf            = channel.empty()
    ch_plot_homer_annotatepeaks_txt = channel.empty()
    ch_plot_homer_annotatepeaks_pdf = channel.empty()
    ch_plot_homer_annotatepeaks_tsv = channel.empty()
    if (!skip_peak_annotation) {
        //
        // Annotate peaks with HOMER
        //
        HOMER_ANNOTATEPEAKS (
            ch_mace_peaks,
            ch_fasta.map { it -> it[1] },
            ch_gtf.map { it -> it[1] }
        )
        ch_homer_annotatepeaks = HOMER_ANNOTATEPEAKS.out.txt

        if (!skip_peak_qc) {

            // Create channels: [ meta, [ peaks ] ]
            // Where meta = [ id:exp_type, exp_type:exp_type ]
            // ch_mace_peaks
            //     .map {
            //         meta, peaks ->
            //             [ meta.exp_type, meta.genome, peaks ]
            //     }
            //     .groupTuple(by: [0, 1])
            //     .map {
            //         exp_type, genome, peaks ->
            //             def meta_new = [:]
            //             meta_new.id = exp_type
            //             meta_new.exp_type = exp_type
            //             meta_new.genome = genome
            //             [ meta_new, peaks ]
            //     }
            //     .set { ch_mace_peaks_grouped }

            // Create channels: [ meta, [ anns ] ]
            // Where meta = [ id:exp_type, exp_type:exp_type ]
            HOMER_ANNOTATEPEAKS.out.txt
                .map { meta, anns ->
                    [meta.exp_type, meta.genome, anns]
                }
                .groupTuple(by: [0, 1])
                .map { exp_type, genome, anns ->
                    def meta_new = [:]
                    meta_new.id = exp_type
                    meta_new.exp_type = exp_type
                    meta_new.genome = genome
                    [meta_new, anns]
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
    peaks                        = ch_mace_peaks // channel: [ val(meta), [ peaks ] ]
    border                       = MACE_MACE.out.border // channel: [ val(meta), [ bed ] ]
    border_cluster               = MACE_MACE.out.border_cluster // channel: [ val(meta), [ bed ] ]
    border_pair_elite            = MACE_MACE.out.border_pair_elite // channel: [ val(meta), [ bed ] ]
    border_pair                  = MACE_MACE.out.border_pair // channel: [ val(meta), [ bed ] ]
    border_pair_peak             = MACE_MACE.out.border_pair_peak // channel: [ val(meta), [ bed ] ]
    frip_txt                     = FRIP_SCORE.out.txt // channel: [ val(meta), [ txt ] ]
    frip_multiqc                 = MULTIQC_CUSTOM_PEAKS.out.frip // channel: [ val(meta), [ frip ] ]
    peak_count_multiqc           = MULTIQC_CUSTOM_PEAKS.out.count // channel: [ val(meta), [ counts ] ]
    homer_annotatepeaks          = ch_homer_annotatepeaks // channel: [ val(meta), [ txt ] ]
    plot_homer_annotatepeaks_txt = ch_plot_homer_annotatepeaks_txt // channel: [ txt ]
    plot_homer_annotatepeaks_pdf = ch_plot_homer_annotatepeaks_pdf // channel: [ pdf ]
    plot_homer_annotatepeaks_tsv = ch_plot_homer_annotatepeaks_tsv // channel: [ tsv ]
}
