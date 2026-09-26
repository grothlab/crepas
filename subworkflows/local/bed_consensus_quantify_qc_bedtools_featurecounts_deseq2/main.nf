
//
// Call consensus peaks with BEDTools and custom scripts, annotate with HOMER, quantify with featureCounts and QC with DESeq2
//

include { HOMER_ANNOTATEPEAKS    } from '../../../modules/nf-core/homer/annotatepeaks/main'
include { SUBREAD_FEATURECOUNTS  } from '../../../modules/nf-core/subread/featurecounts/main'

include { MACS3_CONSENSUS        } from '../../../modules/local/macs3_consensus/main'
include { ANNOTATE_BOOLEAN_PEAKS } from '../../../modules/local/annotate_boolean_peaks/main'
include { DESEQ2_QC              } from '../../../modules/local/deseq2_qc/main'

include { DEEPTOOLS_COMPUTEMATRIX as DEEPTOOLS_COMPUTEMATRIX_PEAKS } from '../../../modules/nf-core/deeptools/computematrix/main'
include { DEEPTOOLS_PLOTPROFILE as DEEPTOOLS_PLOTPROFILE_PEAKS } from '../../../modules/nf-core/deeptools/plotprofile/main'
include { DEEPTOOLS_PLOTHEATMAP as DEEPTOOLS_PLOTHEATMAP_PEAKS } from '../../../modules/nf-core/deeptools/plotheatmap/main'


workflow BED_CONSENSUS_QUANTIFY_QC_BEDTOOLS_FEATURECOUNTS_DESEQ2 {
    take:
    ch_peaks                            // channel: [ val(meta), [ peaks ] ]
    ch_bam                              // channel: [ val(meta), ip_bam ]
    ch_bigwigs                          // channel: [ val(meta), [ bigwigs ] ]
    ch_fasta                            // channel: [ fasta ]
    ch_gtf                              // channel: [ gtf ]
    ch_deseq2_pca_header_multiqc        // channel: [ header_file ]
    ch_deseq2_clustering_header_multiqc // channel: [ header_file ]
    is_narrow_peak                      // boolean: true/false
    skip_peak_annotation                // boolean: true/false
    skip_deseq2_qc                      // boolean: true/false
    skip_consensus_plotprofile          // boolean: true/false
    input_cisrpm_in_plotprofile       // boolean: true/false

    main:

    ch_multiqc_files = channel.empty()

    //TODO: print to fiule for debugging
    ch_peaks
        .map {
            meta, peaks ->
                "${meta}\t${peaks}"
        }
        .collectFile( name: 'ch_peaks.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BED_CONSENSUS_QUANTIFY_QC_BEDTOOLS_FEATURECOUNTS_DESEQ2" )


    // Create channels: [ meta , [ peaks ] ]
    // Where meta = [ id:antibody, exp_type, multiple_groups:true/false, replicates_exist:true/false ]
    ch_peaks
        .map {
            meta, peak ->
                [ meta.antibody, meta.exp_type, meta.id - ~/_bRep_.*$/, meta, peak ]
        }
        .groupTuple(by: [0, 1])
            .map {
                antibody, exp_type, groups, metas, peaks ->
                [
                    antibody,
                    exp_type,
                    metas,
                    groups.groupBy().collectEntries { [(it.key) : it.value.size()] },
                    peaks
                ]
            }
            .map {
                antibody, exp_type, metas, groups, peaks ->
                def meta_clone = metas[0].clone()
                if (exp_type == 'ATAC-seq') {
                    meta_clone.id = exp_type
                } else if (antibody) {
                    meta_clone.id = exp_type + '_' + antibody
                } else {
                    meta_clone.id = exp_type + '_no_antibody_'
                }
                meta_clone.multiple_groups = groups.size() > 1
                meta_clone.replicates_exist = groups.max { it -> it.value }.value > 1
                [ meta_clone, peaks ]
        }
        .set { ch_antibody_peaks }


    //
    // Generate consensus peaks across samples
    //
    MACS3_CONSENSUS (
        ch_antibody_peaks,
        is_narrow_peak
    )

    //
    // Annotate consensus peaks
    //
    if (!skip_peak_annotation) {
        HOMER_ANNOTATEPEAKS (
            MACS3_CONSENSUS.out.bed,
            ch_fasta,
            ch_gtf
        )

        //
        // MODULE: Add boolean fields to annotated consensus peaks to aid filtering
        //
        ANNOTATE_BOOLEAN_PEAKS (
            MACS3_CONSENSUS.out.boolean_txt.join(HOMER_ANNOTATEPEAKS.out.txt, by: [0])
        )
    }

    // Create ID from BAM's meta to match consensus peaks ID
    ch_bam
        .map {
            meta, bam ->
                [ meta.antibody, meta.exp_type, meta.single_end, bam ]
        }
        .groupTuple(by: [0, 1, 2])
        .map {
            antibody, exp_type, se, bams ->
            def id
            if (exp_type == 'ATAC-seq') {
                id = exp_type
            } else if (antibody) {
                id = exp_type + '_' + antibody
            } else {
                id = exp_type + '_no_antibody_'
            }
            [ id, se, bams ]
        }
        .set { ch_antibody_bams }


    // Create channels: [ meta, [ ip_bams ], saf ]
    MACS3_CONSENSUS
        .out
        .saf
        .map {
            meta, saf ->
                [ meta.id, meta, saf ]
        }
        .set { ch_id_saf }

    // Split SE and PE samples that share consensus peaks, since featureCounts cannot handle mixed layouts
    ch_antibody_bams_se = ch_antibody_bams.filter { id, se, bams -> se }.map { id, se, bams -> [ id, bams ]}
    ch_antibody_bams_pe = ch_antibody_bams.filter { id, se, bams -> !se }.map { id, se, bams -> [ id, bams ]}

    ch_id_saf
        .join(ch_antibody_bams_se)
        .map {
            _id, meta, saf, bams ->
                def meta_clone = meta.clone()
                meta_clone.single_end = true
                [ meta_clone, bams.flatten().sort(), saf ]
        }
        .set { ch_bam_saf_se }

    ch_id_saf
        .join(ch_antibody_bams_pe)
        .map {
            _id, meta, saf, bams ->
                def meta_clone = meta.clone()
                meta_clone.single_end = false
                [ meta_clone, bams.flatten().sort(), saf ]
        }
        .set { ch_bam_saf_pe }

    //
    // Quantify peaks across samples with featureCounts
    //
    SUBREAD_FEATURECOUNTS (
        ch_bam_saf_se.mix(ch_bam_saf_pe)
    )
    ch_multiqc_files = ch_multiqc_files.mix(SUBREAD_FEATURECOUNTS.out.summary.collect { it -> it[1] })

    //
    // Generate QC plots with DESeq2
    //
    ch_deseq2_qc_pdf           = channel.empty()
    ch_deseq2_qc_rdata         = channel.empty()
    ch_deseq2_qc_rds           = channel.empty()
    ch_deseq2_qc_pca_txt       = channel.empty()
    ch_deseq2_qc_dists_txt     = channel.empty()
    ch_deseq2_qc_log           = channel.empty()
    ch_deseq2_qc_size_factors  = channel.empty()
    if (!skip_deseq2_qc) {
        DESEQ2_QC (
            SUBREAD_FEATURECOUNTS.out.counts,
            ch_deseq2_pca_header_multiqc,
            ch_deseq2_clustering_header_multiqc
        )
        ch_deseq2_qc_pdf           = DESEQ2_QC.out.pdf
        ch_deseq2_qc_rdata         = DESEQ2_QC.out.rdata
        ch_deseq2_qc_rds           = DESEQ2_QC.out.rds
        ch_deseq2_qc_pca_txt       = DESEQ2_QC.out.pca_txt
        ch_deseq2_qc_dists_txt     = DESEQ2_QC.out.dists_txt
        ch_deseq2_qc_log           = DESEQ2_QC.out.log
        ch_deseq2_qc_size_factors  = DESEQ2_QC.out.size_factors
        ch_multiqc_files = ch_multiqc_files.mix(DESEQ2_QC.out.pca_multiqc.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(DESEQ2_QC.out.dists_multiqc.collect { it -> it[1] })
    }

    if (!skip_consensus_plotprofile) {
        MACS3_CONSENSUS
            .out
            .bed
            .map { meta, peaks ->
                [ meta.antibody, meta.exp_type, peaks ]
            }
            .set { ch_cons_peaks }

        if (!input_cisrpm_in_plotprofile) {
            ch_bigwigs = ch_bigwigs.filter { meta, bws -> !(meta.is_input_control && meta.norm_factor_type == 'cisrpm') }
        }

        ch_bigwigs
            .map { meta, bw ->
                def antibody = meta.antibody ?: meta.input_control_of_antibody
                [ antibody, meta.exp_type, meta.norm_factor_type, meta.signal_vs_input_operation, meta.averaged_brep, meta.id, meta, bw ]
            }
            .groupTuple(by: [0, 1, 2, 3, 4])
            .map { antibody, exp_type, norm_factor_type, signal_vs_input_op, averaged_brep, ids, metas, bws ->
                // Sort ids, metas and bws by id to ensure consistent order in plots
                def sorted_ids = ids.sort()
                def sorted_metas = metas.sort { meta -> meta.id }
                def sorted_bws = bws.sort { bw -> bw.name }
                [ antibody, exp_type, norm_factor_type, signal_vs_input_op, averaged_brep, sorted_ids, sorted_metas, sorted_bws ]
            }
            .combine(ch_cons_peaks, by: [0, 1])
            .map {
                antibody, exp_type, norm_factor_type, signal_vs_input_op, averaged_brep, ids, metas, bws, cons_peaks ->
                    def meta_clone = metas[0].clone()
                    meta_clone.id = exp_type + '_' +
                        (antibody ? antibody : 'no_antibody') +
                        '_' + norm_factor_type +
                        (signal_vs_input_op ? '_' + signal_vs_input_op : '') +
                        (averaged_brep ? '_' + 'bRep_avg' : '')
                    meta_clone.antibody = antibody
                    meta_clone.ids = ids
                    [ meta_clone, bws.flatten(), cons_peaks ]
            }
            .set { ch_bigwigs_peaks }

        //
        // MODULE: deepTools computeMatrix for peaks
        //
        DEEPTOOLS_COMPUTEMATRIX_PEAKS (
            ch_bigwigs_peaks
        )

        //
        // MODULE: deepTools profile plots
        //
        DEEPTOOLS_PLOTPROFILE_PEAKS (
            DEEPTOOLS_COMPUTEMATRIX_PEAKS.out.matrix
        )
        ch_multiqc_files = ch_multiqc_files.mix(DEEPTOOLS_PLOTPROFILE_PEAKS.out.table.collect { it -> it[1] })

        //
        // MODULE: deepTools heatmaps
        //
        DEEPTOOLS_PLOTHEATMAP_PEAKS (
            DEEPTOOLS_COMPUTEMATRIX_PEAKS.out.matrix
        )
    }


    emit:
    consensus_bed           = MACS3_CONSENSUS.out.bed           // channel: [ bed ]
    consensus_saf           = MACS3_CONSENSUS.out.saf           // channel: [ saf ]
    consensus_pdf           = MACS3_CONSENSUS.out.pdf           // channel: [ pdf ]
    consensus_txt           = MACS3_CONSENSUS.out.txt           // channel: [ pdf ]
    consensus_boolean_txt   = MACS3_CONSENSUS.out.boolean_txt   // channel: [ txt ]
    consensus_intersect_txt = MACS3_CONSENSUS.out.intersect_txt // channel: [ txt ]

    featurecounts_txt       = SUBREAD_FEATURECOUNTS.out.counts  // channel: [ txt ]
    featurecounts_summary   = SUBREAD_FEATURECOUNTS.out.summary // channel: [ txt ]

    deseq2_qc_pdf           = ch_deseq2_qc_pdf                  // channel: [ pdf ]
    deseq2_qc_rdata         = ch_deseq2_qc_rdata                // channel: [ rdata ]
    deseq2_qc_rds           = ch_deseq2_qc_rds                  // channel: [ rds ]
    deseq2_qc_pca_txt       = ch_deseq2_qc_pca_txt              // channel: [ txt ]
    deseq2_qc_dists_txt     = ch_deseq2_qc_dists_txt            // channel: [ txt ]
    deseq2_qc_log           = ch_deseq2_qc_log                  // channel: [ txt ]
    deseq2_qc_size_factors  = ch_deseq2_qc_size_factors         // channel: [ txt ]

    multiqc_files           = ch_multiqc_files                   // channel: [ multiqc_files ]
}
