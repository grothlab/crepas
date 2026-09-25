
//
// Call peaks with epic2, annotate with HOMER and perform downstream QC
//

include { EPIC2           } from '../../../modules/local/epic2/epic2/main'
include { HOMER_ANNOTATEPEAKS      } from '../../../modules/nf-core/homer/annotatepeaks/main'
include { FRIP_SCORE               } from '../../../modules/local/frip_score/main'
include { MULTIQC_CUSTOM_PEAKS     } from '../../../modules/local/multiqc_custom_peaks/main'
//include { PLOT_MACS3_QC as PLOT_EPIC2_QC } from '../../../modules/local/plot_macs3_qc/main'
include { PLOT_HOMER_ANNOTATEPEAKS } from '../../../modules/local/plot_homer_annotatepeaks/main'

workflow BAM_PEAKS_CALL_QC_ANNOTATE_EPIC2_HOMER {
    take:
    ch_bam_bai                            // channel: [ val(meta), bam, bai ]
    ch_fasta                          // channel: [ fasta  ]
    ch_gtf                            // channel: [ gtf ]
    ch_chrom_sizes                     // channel: [ bed ]
    efective_gfraction                // float: effective genome fraction
    annotate_peaks_suffix             //  string: suffix for input HOMER annotate peaks files to be trimmed off
    ch_peak_count_header_multiqc      // channel: [ header_file ]
    ch_frip_score_multiqc             // channel: [ header_file ]
    ch_peak_annotation_header_multiqc // channel: [ header_file ]
    skip_peak_annotation              // boolean: true/false
    skip_peak_qc                      // boolean: true/false

    main:

    // Branch channels based on if input control is present
    ch_bam_bai
        .branch { meta, bam, bai ->
            ips_with_ipcontrol: meta.input_control
                return [ meta.input_control, meta.antibody, meta, bam, bai ]
            ips_wo_ipcontrol: !meta.input_control && !meta.is_input_control
                return [ meta.id, meta.antibody, meta, bam, bai ]
            ipcontrols: !meta.input_control && meta.is_input_control
                return [ meta.id, meta.input_control_of_antibody, meta, bam, bai ]
        }
        .set { ch_bam_by_type }

    // Create channel: [ meta, [ip_bams_merged_reps], [ip_bais_merged_reps], [ipcontrol_bams_merged_reps], [ipcontrol_bais_merged_reps] ]
    ch_bam_by_type
        .ips_with_ipcontrol
        .combine(ch_bam_by_type.ipcontrols, by: [0, 1])
        .map { ipcontrol_id, antibody, ip_meta, ip_bam, ip_bai, ipcontrol_meta, ipcontrol_bam, ipcontrol_bai ->
            [ ipcontrol_id, antibody, ip_meta, ip_bam, ip_bai, ipcontrol_bam, ipcontrol_bai ]
        }
        .mix(ch_bam_by_type.ips_wo_ipcontrol)
        // ips_wo_ipcontrol do not have ipcontrol_bam (it[5]) and ipcontrol_bai (it[6])
        .map { it ->
            def meta_clone = it[2].clone()
            meta_clone.id = meta_clone.id - ~/_bRep_.*$/
            meta_clone.input_control = meta_clone.input_control - ~/_bRep_.*$/
            [ meta_clone.id, it[1], meta_clone, it[3], it[4], it[5] ?: [], it[6] ?: [] ]
        }
        .groupTuple(by: [0, 1])
        .map {
            id, antibody, metas, ip_bams, ip_bais, ipcontrol_bams, ipcontrol_bais ->
                [ metas[0], ip_bams.flatten(), ip_bais.flatten(), ipcontrol_bams.flatten(), ipcontrol_bais.flatten() ]
        }
        .set { ch_ip_control_bam_bai_merged_reps }

    // TODO: Print to file for debuggin
    ch_ip_control_bam_bai_merged_reps
        .map {
            meta, ip_bams, ip_bais, ipcontrol_bams, ipcontrol_bais ->
                "${meta}\t${ip_bams}\t${ip_bais}\t${ipcontrol_bams}\t${ipcontrol_bais}"
        }
        .collectFile( name: 'ch_ip_control_bam_merged_reps.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_PEAKS_CALL_QC_ANNOTATE_EPIC2_HOMER" )


    //
    // Call peaks with epic2
    //
    EPIC2 (
        ch_ip_control_bam_bai_merged_reps,
        ch_chrom_sizes,
        efective_gfraction

    )

    //
    // Filter out samples with 0 epic2 peaks called
    //
    EPIC2
        .out
        .peak
        .filter {
            meta, peaks ->
                peaks.size() > 0
        }
        .set { ch_epic2_peaks }

    // Create channels: [ meta, ip_bam, peaks ]
    ch_bam_bai
        .join(ch_epic2_peaks, by: 0)
        .map {
            meta, ip_bam, control_bam, peaks ->
                [ meta, ip_bam, peaks ]
        }
        // Split channel by ip_bam
        .transpose()
        .map {
            meta, ip_bam, peaks ->
                def meta_clone = meta.clone()
                meta_clone.id = ip_bam.getSimpleName()
                [ meta_clone, ip_bam, peaks ]
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
    //ch_plot_epic2_qc_txt            = channel.empty()
    //ch_plot_epic2_qc_pdf            = channel.empty()
    ch_plot_homer_annotatepeaks_txt = channel.empty()
    ch_plot_homer_annotatepeaks_pdf = channel.empty()
    ch_plot_homer_annotatepeaks_tsv = channel.empty()
    if (!skip_peak_annotation) {
        //
        // Annotate peaks with HOMER
        //
        HOMER_ANNOTATEPEAKS (
            ch_epic2_peaks,
            ch_fasta.map{ it -> it[1] },
            ch_gtf.map{ it -> it[1] }
        )
        ch_homer_annotatepeaks = HOMER_ANNOTATEPEAKS.out.txt

        if (!skip_peak_qc) {

            // Create channels: [ meta, [ peaks ] ]
            // Where meta = [ id:exp_type, exp_type:exp_type ]
            // ch_epic2_peaks
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
            //     .set { ch_epic2_peaks_grouped }

            //
            // epic2 QC plots with R
            //
            // PLOT_EPIC2_QC (
            //     ch_epic2_peaks_grouped,
            //     is_narrow_peak
            // )
            // ch_plot_epic2_qc_txt = PLOT_EPIC2_QC.out.txt
            // ch_plot_epic2_qc_pdf = PLOT_EPIC2_QC.out.pdf
            // ch_versions = ch_versions.mix(PLOT_EPIC2_QC.out.versions)

            // Create channels: [ meta, [ anns ] ]
            // Where meta = [ id:exp_type, exp_type:exp_type ]
            HOMER_ANNOTATEPEAKS.out.txt
                .map {
                    meta, anns ->
                        [ meta.exp_type, meta.genome, anns ]
                }
                .groupTuple(by: [0, 1])
                .map {
                    exp_type, genome, anns ->
                        def meta_new = [:]
                        meta_new.id = exp_type
                        meta_new.exp_type = exp_type
                        meta_new.genome = genome
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

    peaks                        = ch_epic2_peaks                   // channel: [ val(meta), [ peaks ] ]
    frip_txt                     = FRIP_SCORE.out.txt               // channel: [ val(meta), [ txt ] ]
    frip_multiqc                 = MULTIQC_CUSTOM_PEAKS.out.frip    // channel: [ val(meta), [ frip ] ]

    peak_count_multiqc           = MULTIQC_CUSTOM_PEAKS.out.count   // channel: [ val(meta), [ counts ] ]
    homer_annotatepeaks          = ch_homer_annotatepeaks           // channel: [ val(meta), [ txt ] ]

    //plot_epic2_qc_txt               = ch_plot_epic2_qc_txt             // channel: [ txt ]
    //plot_epic2_qc_pdf               = ch_plot_epic2_qc_pdf             // channel: [ pdf ]

    plot_homer_annotatepeaks_txt = ch_plot_homer_annotatepeaks_txt  // channel: [ txt ]
    plot_homer_annotatepeaks_pdf = ch_plot_homer_annotatepeaks_pdf  // channel: [ pdf ]
    plot_homer_annotatepeaks_tsv = ch_plot_homer_annotatepeaks_tsv  // channel: [ tsv ]
}
