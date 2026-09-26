include { BAM_PEAKS_CALL_QC_ANNOTATE_EPIC2_HOMER                        } from '../../../subworkflows/local/bam_peaks_call_qc_annotate_epic2_homer/main'
include { BAM_PEAKS_CALL_QC_ANNOTATE_MACS3_HOMER                        } from '../../../subworkflows/local/bam_peaks_call_qc_annotate_macs3_homer/main'
include { BAM_PEAKS_CALL_QC_ANNOTATE_GENRICH_HOMER                      } from '../../../subworkflows/local/bam_peaks_call_qc_annotate_genrich_homer/main'
include { BAM_PEAKS_CALL_QC_ANNOTATE_MACE_HOMER                         } from '../../../subworkflows/local/bam_peaks_call_qc_annotate_mace_homer/main'
include { BAM_PEAKS_CALL_QC_ANNOTATE_DANPOS2_HOMER                      } from '../../../subworkflows/local/bam_peaks_call_qc_annotate_danpos2_homer/main'
include { BAM_PEAKS_CALL_QC_ANNOTATE_SEACR_HOMER                        } from '../../../subworkflows/local/bam_peaks_call_qc_annotate_seacr_homer/main'
include { BAM_PEAKS_CALL_QC_ANNOTATE_CONSENRICH_ROCCO_HOMER             } from '../../../subworkflows/local/bam_peaks_call_qc_annotate_consenrich_rocco_homer/main'
include { EDD                                                           } from '../../../modules/local/edd/main'
include { DENOPA                                                        } from '../../../modules/local/denopa/main'
include { BED_CONSENSUS_QUANTIFY_QC_BEDTOOLS_FEATURECOUNTS_DESEQ2       } from '../../../subworkflows/local/bed_consensus_quantify_qc_bedtools_featurecounts_deseq2/main'

workflow CALL_PEAKS {

    take:
    ch_bam_index              // channel: [ val(meta), [ bam ], [ index ] ]
    ch_bigwig_norm
    ch_bedgraph_for_seacr
    ch_peak_callers         // channel: [ 'macs3', 'genrich' ]
    ch_fasta_fai            // channel: [ val(meta), path(fasta), path(fai) ]
    ch_gtf
    ch_effective_gfraction
    ch_chrom_sizes_endo
    ch_blacklist
    ch_sparsebed
    ch_active_regions
    ch_rocco_params
    ch_effective_gsize
    ch_epic2_peak_count_header
    ch_epic2_frip_score_header
    ch_epic2_peak_annotation_header
    ch_gr_peak_count_header
    ch_gr_frip_score_header
    ch_gr_peak_annotation_header
    ch_mace_peak_count_header
    ch_mace_frip_score_header
    ch_mace_peak_annotation_header
    ch_macs3_peak_count_header
    ch_macs3_frip_score_header
    ch_macs3_peak_annotation_header
    ch_seacr_peak_count_header
    ch_seacr_frip_score_header
    ch_seacr_peak_annotation_header
    ch_deseq2_pca_header
    ch_deseq2_clustering_header
    narrow_peak
    skip_peak_annotation
    skip_peak_qc
    skip_bdgcmp
    skip_consensus_peaks
    skip_deseq2_qc
    skip_consensus_plotprofile
    input_cisrpm_in_plotprofile
    seacr_peak_threshold

    main:

    ch_multiqc_files = channel.empty()

    // Create channels per peak_caller
    ch_bam_index
        .combine(ch_peak_callers)
        .map { meta, bam, index, peak_caller ->
            [ meta + [ 'peak_caller': peak_caller ], bam, index ]
        }
        .set { ch_bam_index }

    ch_bedgraph_for_seacr
        .combine(ch_peak_callers)
        .map { meta, bdg, peak_caller ->
            [ meta + [ 'peak_caller': peak_caller ], bdg ]
        }
        .filter { it -> it[0].peak_caller == 'seacr' }
        .set { ch_bedgraph_for_seacr }

    // Create channels without index
    ch_bam = ch_bam_index.map { meta, bam, index -> [ meta, bam ] }
    ch_fasta = ch_fasta_fai.map { meta, fasta, _fai -> [ meta, fasta ] }

    //
    // SUBWORKFLOW: Call consensus regions with Consenrich and ROCCO
    //
    ch_consenrich_tracks = channel.empty()
    ch_rocco_peaks = channel.empty()
    BAM_PEAKS_CALL_QC_ANNOTATE_CONSENRICH_ROCCO_HOMER (
        ch_bam_index.filter { it -> it[0].peak_caller == 'consenrich' },
        ch_chrom_sizes_endo,
        ch_blacklist,
        ch_sparsebed.ifEmpty([[:], []]),
        ch_active_regions.ifEmpty([[:], []]),
        ch_rocco_params,
        ch_effective_gsize
    )
    ch_consenrich_tracks = BAM_PEAKS_CALL_QC_ANNOTATE_CONSENRICH_ROCCO_HOMER.out.consenrich_signal
    ch_consenrich_tracks = ch_consenrich_tracks.mix(BAM_PEAKS_CALL_QC_ANNOTATE_CONSENRICH_ROCCO_HOMER.out.consenrich_residuals)
    ch_consenrich_tracks = ch_consenrich_tracks.mix(BAM_PEAKS_CALL_QC_ANNOTATE_CONSENRICH_ROCCO_HOMER.out.consenrich_eratio)
    ch_rocco_peaks = BAM_PEAKS_CALL_QC_ANNOTATE_CONSENRICH_ROCCO_HOMER.out.rocco_peaks

    //
    // SUBWORKFLOW: Call peaks with epic2, annotate with HOMER and perform downstream QC
    //
    ch_epic2_peaks = channel.empty()
    ch_epic2_frip_multiqc = channel.empty()
    ch_epic2_peak_count_multiqc = channel.empty()
    ch_epic2_plot_homer_annotatepeaks_tsv = channel.empty()
    BAM_PEAKS_CALL_QC_ANNOTATE_EPIC2_HOMER (
        ch_bam_index.filter { it -> it[0].peak_caller == 'epic2' },
        ch_fasta,
        ch_gtf,
        ch_chrom_sizes_endo,
        ch_effective_gfraction,
        ".annotatePeaks.txt",
        ch_epic2_peak_count_header,
        ch_epic2_frip_score_header,
        ch_epic2_peak_annotation_header,
        skip_peak_annotation,
        skip_peak_qc
    )
    ch_epic2_peaks = BAM_PEAKS_CALL_QC_ANNOTATE_EPIC2_HOMER.out.peaks
    ch_epic2_frip_multiqc = BAM_PEAKS_CALL_QC_ANNOTATE_EPIC2_HOMER.out.frip_multiqc
    ch_epic2_peak_count_multiqc = BAM_PEAKS_CALL_QC_ANNOTATE_EPIC2_HOMER.out.peak_count_multiqc
    ch_epic2_plot_homer_annotatepeaks_tsv = BAM_PEAKS_CALL_QC_ANNOTATE_EPIC2_HOMER.out.plot_homer_annotatepeaks_tsv

    //
    // SUBWORKFLOW: Call peaks with Genrich, annotate with HOMER and perform downstream QC
    //
    ch_genrich_peaks = channel.empty()
    BAM_PEAKS_CALL_QC_ANNOTATE_GENRICH_HOMER (
        ch_bam.filter { it -> it[0].peak_caller == 'genrich' },
        ch_fasta_fai,
        ch_gtf,
        ch_blacklist,
        ".annotatePeaks.txt",
        ch_gr_peak_count_header,
        ch_gr_frip_score_header,
        ch_gr_peak_annotation_header,
        narrow_peak,
        skip_peak_annotation,
        skip_peak_qc
    )
    ch_genrich_peaks = BAM_PEAKS_CALL_QC_ANNOTATE_GENRICH_HOMER.out.peaks
    ch_multiqc_files = ch_multiqc_files.mix(BAM_PEAKS_CALL_QC_ANNOTATE_GENRICH_HOMER.out.frip_multiqc.collect { it -> it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(BAM_PEAKS_CALL_QC_ANNOTATE_GENRICH_HOMER.out.peak_count_multiqc.collect { it -> it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(BAM_PEAKS_CALL_QC_ANNOTATE_GENRICH_HOMER.out.plot_homer_annotatepeaks_tsv.collect { it -> it[1] })

    //
    // SUBWORKFLOW: Call peaks with MACE (for ChIP-exo samples)
    //
    ch_mace_peaks = channel.empty()
    BAM_PEAKS_CALL_QC_ANNOTATE_MACE_HOMER(
        ch_bam_index.filter { it -> it[0].peak_caller == 'mace' },
        ch_fasta,
        ch_gtf,
        ch_chrom_sizes_endo,
        ".annotatePeaks.txt",
        ch_mace_peak_count_header,
        ch_mace_frip_score_header,
        ch_mace_peak_annotation_header,
        skip_peak_annotation,
        skip_peak_qc
    )
    ch_mace_peaks = BAM_PEAKS_CALL_QC_ANNOTATE_MACE_HOMER.out.peaks
    ch_multiqc_files = ch_multiqc_files.mix(BAM_PEAKS_CALL_QC_ANNOTATE_MACE_HOMER.out.frip_multiqc.collect { it -> it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(BAM_PEAKS_CALL_QC_ANNOTATE_MACE_HOMER.out.peak_count_multiqc.collect { it -> it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(BAM_PEAKS_CALL_QC_ANNOTATE_MACE_HOMER.out.plot_homer_annotatepeaks_tsv.collect { it -> it[1] })

    //
    // SUBWORKFLOW: Call peaks with DANPOS2
    //
    BAM_PEAKS_CALL_QC_ANNOTATE_DANPOS2_HOMER (
        ch_bam.filter { it -> ['dpeak', 'dpos', 'dregion', 'dtriple'].contains(it[0].peak_caller)}
    )

    //
    // Create channel for downstream processes: [ meta, [ ip_bam, ipcontrol_bam ] [ ip_index, ipcontrol_index ] ]
    // (Excluding ips_wo_ipcontrol as they don't need to be compared to anything)
    //
    ch_bam_index
        .branch { meta, bam, index ->
            ips_with_ipcontrol: meta.input_control
                return [meta.input_control, meta.antibody, meta.peak_caller, meta, bam, index]
            ips_wo_ipcontrol: !meta.input_control && !meta.is_input_control
                return [meta, bam, index]
            ipcontrols: !meta.input_control && meta.is_input_control
                return [meta.id, meta.input_control_of_antibody, meta.peak_caller, meta, bam, index]
        }
        .set { ch_bam_index_by_type }

    ch_bam_index_by_type
        .ips_with_ipcontrol
        .combine(ch_bam_index_by_type.ipcontrols, by: [0,1,2])
        .map { _ipcontrol_id, _antibody, _peak_caller, ip_meta, ip_bam, ip_index, _ipcontrol_meta, ipcontrol_bam, ipcontrol_index ->
            [ ip_meta, [ip_bam] + [ipcontrol_bam], [ip_index] + [ipcontrol_index] ]
        }
        .set { ch_ip_and_ipcontrols_bam_index }

    // Create channels: [ meta, ip_bam, ipcontrol_bam ]
    // Including ips_wo_ipcontrol as they will be used for peak calling without control
    ch_bam_index_by_type
        .ips_wo_ipcontrol
        .map { meta, bam, index -> [meta, [bam], [index]] }
        .mix(ch_ip_and_ipcontrols_bam_index)
        // ips_wo_ipcontrol do not have ipcontrol_bam
        .map { meta, bams, indices ->
            [meta, bams[0], (bams[1] ?: [])]
        }
        .set { ch_all_ip_and_controls }


    ch_edd_peaks = channel.empty()
    //
    // MODULE: Call peaks with EDD
    //
    EDD (
        ch_all_ip_and_controls.filter { it -> it[0].peak_caller == 'edd' },
        ch_chrom_sizes_endo,
        ch_blacklist
    )
    ch_edd_peaks = EDD.out.peaks


    ch_denopa_peaks = channel.empty()
    //
    // MODULE: Call peaks with denopa
    //
    DENOPA (
        ch_bam_index.filter { it -> it[0].peak_caller == 'denopa' }
    )
    ch_denopa_peaks = DENOPA.out.arers


    ch_macs3_peaks = channel.empty()
    //
    // SUBWORKFLOW: Call peaks with MACS3, annotate with HOMER and perform downstream QC
    //
    BAM_PEAKS_CALL_QC_ANNOTATE_MACS3_HOMER (
        ch_all_ip_and_controls.filter { it -> it[0].peak_caller == 'macs3' },
        ch_fasta,
        ch_gtf,
        ch_chrom_sizes_endo,
        ch_effective_gsize,
        "_peaks.annotatePeaks.txt",
        ch_macs3_peak_count_header,
        ch_macs3_frip_score_header,
        ch_macs3_peak_annotation_header,
        narrow_peak,
        skip_peak_annotation,
        skip_peak_qc,
        skip_bdgcmp
    )
    ch_macs3_peaks = BAM_PEAKS_CALL_QC_ANNOTATE_MACS3_HOMER.out.peaks
    ch_multiqc_files = ch_multiqc_files.mix(BAM_PEAKS_CALL_QC_ANNOTATE_MACS3_HOMER.out.frip_multiqc.collect { it -> it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(BAM_PEAKS_CALL_QC_ANNOTATE_MACS3_HOMER.out.peak_count_multiqc.collect { it -> it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(BAM_PEAKS_CALL_QC_ANNOTATE_MACS3_HOMER.out.plot_homer_annotatepeaks_tsv.collect { it -> it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(BAM_PEAKS_CALL_QC_ANNOTATE_MACS3_HOMER.out.xls.collect { it -> it[1] })


    //
    // SUBWORKFLOW: Call peaks with SEACR
    //
    ch_seacr_peaks = channel.empty()

    // Create channels: [ meta, ip_bam, ipcontrol_bam ]
    // Including ips_wo_ipcontrol as they will be used for peak calling without control
    ch_bedgraph_for_seacr
        .branch { meta, bdg ->
            ips_with_ipcontrol: meta.input_control
                return [meta.input_control, meta.antibody, meta.norm_factor_type, meta, bdg]
            ips_wo_ipcontrol: !meta.input_control && !meta.is_input_control
                return [meta, bdg, []]
            ipcontrols: !meta.input_control && meta.is_input_control
                return [meta.id, meta.input_control_of_antibody, meta.norm_factor_type, meta, bdg]
        }
        .set { ch_bedgraph_by_type }

    ch_bedgraph_by_type
        .ips_with_ipcontrol
        .combine(ch_bedgraph_by_type.ipcontrols, by: [0,1,2])
        .map { ipcontrol_id, antibody, norm_factor_type, ip_meta, ip_bdg, ipcontrol_meta, ipcontrol_bdg ->
            [ ip_meta, ip_bdg, ipcontrol_bdg ]
        }
        .set { ch_ip_and_ipcontrols_bdg }

    ch_bedgraph_by_type
        .ips_wo_ipcontrol
        .mix(ch_ip_and_ipcontrols_bdg)
        .set { ch_all_bdg_ip_and_controls }


    BAM_PEAKS_CALL_QC_ANNOTATE_SEACR_HOMER (
        ch_all_bdg_ip_and_controls.filter { it -> it[0].peak_caller == 'seacr' },
        seacr_peak_threshold,
        ch_all_ip_and_controls.filter { it -> it[0].peak_caller == 'seacr' },
        ch_fasta,
        ch_gtf,
        ".annotatePeaks.txt",
        ch_seacr_peak_count_header,
        ch_seacr_frip_score_header,
        ch_seacr_peak_annotation_header,
        skip_peak_annotation,
        skip_peak_qc
    )
    ch_seacr_peaks = BAM_PEAKS_CALL_QC_ANNOTATE_SEACR_HOMER.out.peaks
    ch_multiqc_files = ch_multiqc_files.mix(BAM_PEAKS_CALL_QC_ANNOTATE_SEACR_HOMER.out.frip_multiqc.collect { it -> it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(BAM_PEAKS_CALL_QC_ANNOTATE_SEACR_HOMER.out.peak_count_multiqc.collect { it -> it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(BAM_PEAKS_CALL_QC_ANNOTATE_SEACR_HOMER.out.plot_homer_annotatepeaks_tsv.collect { it -> it[1] })


    //
    //  Consensus peaks analysis
    //
    ch_consensus_bed = channel.empty()
    ch_consensus_txt = channel.empty()
    if (!skip_consensus_peaks) {
        // Create channels: [ meta, ip_bam ]
        ch_all_ip_and_controls
            .filter { it -> it[0].peak_caller == 'macs3' }
            .map { meta, ip_bam, _ipcontrol_bam ->
                [meta, ip_bam]
            }
            .set { ch_bam_for_consensus }

        BED_CONSENSUS_QUANTIFY_QC_BEDTOOLS_FEATURECOUNTS_DESEQ2(
            ch_macs3_peaks,
            ch_bam_for_consensus,
            ch_bigwig_norm,
            ch_fasta.map { it -> it[1] },
            ch_gtf.map { it -> it[1] },
            ch_deseq2_pca_header,
            ch_deseq2_clustering_header,
            narrow_peak,
            skip_peak_annotation,
            skip_deseq2_qc,
            skip_consensus_plotprofile,
            input_cisrpm_in_plotprofile
        )
        ch_consensus_bed = BED_CONSENSUS_QUANTIFY_QC_BEDTOOLS_FEATURECOUNTS_DESEQ2.out.consensus_bed
        ch_consensus_txt = BED_CONSENSUS_QUANTIFY_QC_BEDTOOLS_FEATURECOUNTS_DESEQ2.out.consensus_txt
        ch_multiqc_files = ch_multiqc_files.mix(BED_CONSENSUS_QUANTIFY_QC_BEDTOOLS_FEATURECOUNTS_DESEQ2.out.multiqc_files)
    }

    emit:

    edd_peaks           = ch_edd_peaks              // channel: [ meta, peaks ]
    macs3_peaks         = ch_macs3_peaks
    genrich_peaks       = ch_genrich_peaks           // channel: [ meta, peaks ]
    epic2_peaks         = ch_epic2_peaks            // channel: [ meta
    seacr_peaks         = ch_seacr_peaks            // channel: [ meta, peaks ]
    denopa_peaks        = ch_denopa_peaks           // channel: [ meta,
    consensus_bed       = ch_consensus_bed          // channel: [ antibody, consensus_bed ]
    consensus_txt       = ch_consensus_txt          // channel: [ antibody, consensus_txt ]
    mace_peaks          = ch_mace_peaks             // channel: [ meta, peaks ]
    consenrich_tracks   = ch_consenrich_tracks      // channel: [ meta, consenrich_signal ], [ meta, consenrich_residuals ], [ meta, consenrich_eratio ]
    rocco_peaks         = ch_rocco_peaks            // channel: [ meta,
    multiqc_files       = ch_multiqc_files          // channel: [ path(mqc_files) ]
}
