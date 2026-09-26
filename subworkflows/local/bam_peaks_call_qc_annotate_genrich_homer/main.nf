//
// Call peaks with Genrich, annotate with HOMER and perform downstream QC
//

include { SAMTOOLS_SORT                    } from '../../../modules/nf-core/samtools/sort/main'
include { GENRICH                          } from '../../../modules/nf-core/genrich/main'
include { HOMER_ANNOTATEPEAKS              } from '../../../modules/nf-core/homer/annotatepeaks/main'
include { FRIP_SCORE                       } from '../../../modules/local/frip_score/main'
include { MULTIQC_CUSTOM_PEAKS             } from '../../../modules/local/multiqc_custom_peaks/main'
include { PLOT_MACS3_QC as PLOT_GENRICH_QC } from '../../../modules/local/plot_macs3_qc/main'
include { PLOT_HOMER_ANNOTATEPEAKS         } from '../../../modules/local/plot_homer_annotatepeaks/main'

workflow BAM_PEAKS_CALL_QC_ANNOTATE_GENRICH_HOMER {
    take:
    ch_bam                            // channel: [ val(meta), [ ip_bam ], [ ipcontrol_bam ] ]
    ch_fasta_fai                      // channel: [ val(meta), path(fasta), path(fai) ]
    ch_gtf                            // channel: [ gtf ]
    ch_blacklist                      // channel: [ bed ]
    annotate_peaks_suffix             //  string: suffix for input HOMER annotate peaks files to be trimmed off
    ch_peak_count_header_multiqc      // channel: [ header_file ]
    ch_frip_score_multiqc             // channel: [ header_file ]
    ch_peak_annotation_header_multiqc // channel: [ header_file ]
    is_narrow_peak                    // boolean: true/false
    skip_peak_annotation              // boolean: true/false
    skip_peak_qc                      // boolean: true/false

    main:

    SAMTOOLS_SORT ( ch_bam, ch_fasta_fai, '' )

    SAMTOOLS_SORT
        .out
        .bam
        .branch { meta, bam ->
            ips_with_ipcontrol: meta.input_control
                return [meta.input_control, meta.antibody, meta, bam]
            ips_wo_ipcontrol: !meta.input_control && !meta.is_input_control
                return [meta.id, meta.antibody, meta, bam]
            ipcontrols: !meta.input_control && meta.is_input_control
                return [meta.id, meta.input_control_of_antibody, meta, bam]
        }
        .set { ch_bam_by_type }

    // Create channel: [ meta, [ip_bams_merged_reps], [ipcontrol_bams_merged_reps] ]
    ch_bam_by_type
        .ips_with_ipcontrol
        .combine(ch_bam_by_type.ipcontrols, by: [0, 1])
        .map { ipcontrol_id, antibody, ip_meta, ip_bam, ipcontrol_meta, ipcontrol_bam ->
            [ ipcontrol_id, antibody, ip_meta, ip_bam, ipcontrol_bam ]
        }
        .mix(ch_bam_by_type.ips_wo_ipcontrol)
        .map { it ->
            def meta_clone = it[2].clone()
            meta_clone.id = meta_clone.id - ~/_bRep_.*$/
            meta_clone.input_control = meta_clone.input_control - ~/_bRep_.*$/
            [meta_clone.id, it[1], meta_clone, it[3], it[4] ?: []]
        }
        .groupTuple(by: [0, 1])
        .map { id, antibody, metas, ip_bams, ipcontrol_bams ->
            [metas[0], ip_bams.flatten(), ipcontrol_bams.flatten()]
        }
        .set { ch_ip_ipcontrol_bam_merged_reps }

    // TODO: Print to file for debuggin
    ch_ip_ipcontrol_bam_merged_reps
        .map { meta, ip_bams, ipcontrol_bams ->
            "${meta}\t${ip_bams}\t${ipcontrol_bams}"
        }
        .collectFile(name: 'ch_ip_ipcontrol_bam_merged_reps.txt', newLine: true, sort: false, storeDir: "${params.outdir}")


    //
    // Call peaks with Genrich
    //
    GENRICH (
        ch_ip_ipcontrol_bam_merged_reps,
        ch_blacklist.map { it -> it[1] }
    )

    //
    // Filter out samples with 0 Genrich peaks called
    //
    GENRICH.out.peak.filter { meta, peaks ->
        peaks.size() > 0
    }.set { ch_gr_peaks }

    // Create channels: [ meta, ip_bam, peaks ]
    ch_bam
        .join(ch_gr_peaks, by: 0)
        .map { meta, ip_bam, ipcontrol_bam, peaks ->
            [meta, ip_bam, peaks]
        }
        .transpose()
        .map { meta, ip_bam, peaks ->
            def meta_clone = meta.clone()
            meta_clone.id = ip_bam.getSimpleName()
            [meta_clone, ip_bam, peaks]
        }
        .set { ch_bam_peak }

    //
    // Calculate FRiP score
    //
    FRIP_SCORE(
        ch_bam_peak
    )

    // Create channels: [ meta, peaks, frip ]
    ch_bam_peak
        .join(FRIP_SCORE.out.txt, by: 0)
        .map { meta, ip_bam, peaks, frip ->
            [meta, peaks, frip]
        }
        .set { ch_bam_peak_frip }

    //
    // FRiP score custom content for MultiQC
    //
    MULTIQC_CUSTOM_PEAKS(
        ch_bam_peak_frip,
        ch_peak_count_header_multiqc,
        ch_frip_score_multiqc
    )

    ch_homer_annotatepeaks = channel.empty()
    ch_plot_gr_qc_txt = channel.empty()
    ch_plot_gr_qc_pdf = channel.empty()
    ch_plot_homer_annotatepeaks_txt = channel.empty()
    ch_plot_homer_annotatepeaks_pdf = channel.empty()
    ch_plot_homer_annotatepeaks_tsv = channel.empty()
    if (!skip_peak_annotation) {
        //
        // Annotate peaks with HOMER
        //
        HOMER_ANNOTATEPEAKS(
            ch_gr_peaks,
            ch_fasta_fai.map { it -> it[1] },
            ch_gtf.map { it -> it[1] }
        )
        ch_homer_annotatepeaks = HOMER_ANNOTATEPEAKS.out.txt

        if (!skip_peak_qc) {

            // Create channels: [ meta, [ peaks ] ]
            // Where meta = [ id:exp_type, exp_type:exp_type ]
            ch_gr_peaks
                .map { meta, peaks ->
                    [meta.exp_type, meta.genome, peaks]
                }
                .groupTuple(by: [0, 1])
                .map { exp_type, genome, peaks ->
                    def meta_new = [:]
                    meta_new.id = exp_type
                    meta_new.exp_type = exp_type
                    meta_new.genome = genome
                    [meta_new, peaks]
                }
                .set { ch_gr_peaks_grouped }

            //
            // Genrich QC plots with R
            //
            PLOT_GENRICH_QC(
                ch_gr_peaks_grouped,
                is_narrow_peak
            )
            ch_plot_gr_qc_txt = PLOT_GENRICH_QC.out.txt
            ch_plot_gr_qc_pdf = PLOT_GENRICH_QC.out.pdf

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
            PLOT_HOMER_ANNOTATEPEAKS(
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
    peaks                        = ch_gr_peaks // channel: [ val(meta), [ peaks ] ]
    bed_intervals                = GENRICH.out.bed_intervals // channel: [ val(meta), [ bed ] ]
    bedgraph_pileup              = GENRICH.out.bedgraph_pileup // channel: [ val(meta), [ bedgraph ] ]
    bedgraph_pvalues             = GENRICH.out.bedgraph_pvalues // channel: [ val(meta), [ bedgraph ] ]
    duplicates                   = GENRICH.out.duplicates // channel: [ val(meta), [ bed ] ]
    frip_txt                     = FRIP_SCORE.out.txt // channel: [ val(meta), [ txt ] ]
    frip_multiqc                 = MULTIQC_CUSTOM_PEAKS.out.frip // channel: [ val(meta), [ frip ] ]
    peak_count_multiqc           = MULTIQC_CUSTOM_PEAKS.out.count // channel: [ val(meta), [ counts ] ]
    homer_annotatepeaks          = ch_homer_annotatepeaks // channel: [ val(meta), [ txt ] ]
    plot_gr_qc_txt               = ch_plot_gr_qc_txt // channel: [ txt ]
    plot_gr_qc_pdf               = ch_plot_gr_qc_pdf // channel: [ pdf ]
    plot_homer_annotatepeaks_txt = ch_plot_homer_annotatepeaks_txt // channel: [ txt ]
    plot_homer_annotatepeaks_pdf = ch_plot_homer_annotatepeaks_pdf // channel: [ pdf ]
    plot_homer_annotatepeaks_tsv = ch_plot_homer_annotatepeaks_tsv // channel: [ tsv ]
}
