include { SAMTOOLS_SORT                                         } from '../../../modules/nf-core/samtools/sort/main'
include { BEDTOOLS_BAMTOBED                                     } from '../../../modules/nf-core/bedtools/bamtobed/main'
include { BED_TO_TAGALIGN                                       } from '../../../modules/local/bed_to_tagalign/main'
include { TAGALIGN_SELF_PSEUDOREPLICATES                        } from '../../../modules/local/tagalign_self_pseudoreplicates/main'
include { FIND_CONCATENATE as TAGALIGN_POOL                     } from '../../../modules/nf-core/find/concatenate/main'
include { PHANTOMPEAKQUALTOOLS as PHANTOMPEAKQUALTOOLS_SPP      } from '../../../modules/local/phantompeakqualtools/main'
include { BED_FILTER_BLACKLIST as PEAKS_FILTER_BLACKLIST        } from '../../../modules/local/bed_filter_blacklist/main'
include { IDR                                                   } from '../../../modules/nf-core/idr/main'
include { IDR_FILTER_THRESHOLD                                  } from '../../../modules/local/idr_filter_threshold/main'
include { BED_FILTER_BLACKLIST as CONSENSUS_FILTER_BLACKLIST    } from '../../../modules/local/bed_filter_blacklist/main'
include { PEAKS_NAIVE_OVERLAP                                   } from '../../../modules/local/peaks_naive_overlap/main'
include { TAGALIGN_FRIP_SCORE                                   } from '../../../modules/local/tagalign_frip_score/main'

workflow BAM_ENCODE_PIPELINE {
    take:
    ch_bam                            // channel: [ val(meta), [ ip_bam ], [ control_bam ] ]
    ch_fasta_fai                      // channel: [ val(meta), path(fasta), path(fai) ]
    ch_chromsizes                     // channel: [ val(meta), path(chromsizes) ]
    ctl_depth_ratio_threshold
    val_peak_type
    ch_blacklist
    idr_filtering_threshold
    encode_peak_max_score

    main:

    //
    // MODULE: Name-sorting BAM files
    //
    SAMTOOLS_SORT (
        ch_bam,
        ch_fasta_fai,
        ''
    )

    //
    // MODULE: Convert BAM to BED
    //
    BEDTOOLS_BAMTOBED (
        SAMTOOLS_SORT.out.bam
    )

    //
    // MODULE: Convert BED to TAGALIGN
    //
    BED_TO_TAGALIGN (
        BEDTOOLS_BAMTOBED.out.bed
    )


    BED_TO_TAGALIGN
        .out
        .tagalign
        .filter { meta, tagalign -> !meta.is_input_control } // We do not generate pseudoreplicates for input controls
        .map { meta, tagalign -> [ meta + [ is_pseudoreplicate: true ], tagalign ]}
        .set { ch_tagalign_ips_for_pseudoreps }

    //
    // MODULE: Generate self-pseudoreplicates of IP samples
    //
    TAGALIGN_SELF_PSEUDOREPLICATES (
        ch_tagalign_ips_for_pseudoreps
    )

    // Add pseudoreplicate to metadata
    TAGALIGN_SELF_PSEUDOREPLICATES.out.tagalign1
        .map { meta, tagalign -> [ meta + [ pseudoreplicate: '1' ], tagalign ]}
        .mix(
            TAGALIGN_SELF_PSEUDOREPLICATES.out.tagalign2
                .map { meta, tagalign -> [ meta + [ pseudoreplicate: '2' ], tagalign ] }
        )
        .set {ch_self_pseudoreps}


    // Create channel: [ meta, tagaligns ] to pool replicates and pseudoreplicates
    BED_TO_TAGALIGN
        .out
        .tagalign
        .mix(ch_self_pseudoreps)
        .set { ch_tas_reps_and_pseudoreps }

    ch_tas_reps_and_pseudoreps
        .map { meta, tagalign ->
            def meta_clone = meta.clone()
            meta_clone.id = meta_clone.id - ~/_bRep_.*$/
            meta_clone.input_control = meta_clone.input_control - ~/_bRep_.*$/
            def antibody = meta.input_control_of_antibody ?: meta.antibody
            def pseudoreplicate = meta.pseudoreplicate ?: false
            [  meta_clone.id, antibody, pseudoreplicate, meta_clone, tagalign ]
        }
        .groupTuple(by: [0, 1, 2])
        .map { id, antibody, pseudoreplicate, metas, tagaligns ->
            // Sort metas and tagaligns to ensure consistent ordering for caching and resuming
            def sorted_tagaligns = tagaligns.sort { it -> it.name }
            def sorted_metas = metas.sort { meta -> meta.brep }
            def meta_clone = sorted_metas[0].clone()
            meta_clone.is_pooled = true
            [ meta_clone, sorted_tagaligns ]
        }
        .set { ch_tas_reps_and_pseudoreps_to_pool }

    // TODO: save for debugging
    ch_tas_reps_and_pseudoreps_to_pool
        .map { meta, tagaligns ->
            "${meta}\t${tagaligns}"
        }
        .collectFile(name: 'ch_tas_reps_and_pseudoreps_to_pool.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_ENCODE_PIPELINE")

    //
    // MODULE: Pool replicates and pseudoreplicates with cat
    //
    TAGALIGN_POOL (
        ch_tas_reps_and_pseudoreps_to_pool
    )
    ch_tas_reps_and_pseudoreps_pooled = TAGALIGN_POOL.out.file_out

    // TODO: save for debugging
    ch_tas_reps_and_pseudoreps_pooled
        .map { meta, tagaligns ->
            "${meta}\t${tagaligns}"
        }
        .collectFile(name: 'ch_tas_reps_and_pseudoreps_pooled.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_ENCODE_PIPELINE")

    //
    // Create channel: [ meta, tagalign ] with metadata indicating whether to use pooled control or not for each sample
    //
    ch_tas_reps_and_pseudoreps
        .branch { meta, tagalign ->
            ips_with_ipcontrol: meta.input_control
                return [meta.input_control, meta.antibody, meta, tagalign]
            ips_wo_ipcontrol: !meta.input_control && !meta.is_input_control
                return [meta, tagalign]
            ipcontrols: !meta.input_control && meta.is_input_control
                return [meta.id, meta.input_control_of_antibody, meta, tagalign]
        }
        .set { ch_tas_reps_and_pseudoreps_by_type }

    ch_tas_reps_and_pseudoreps_by_type
        .ips_with_ipcontrol
        .combine(ch_tas_reps_and_pseudoreps_by_type.ipcontrols, by: [0, 1]) // combine by ipcontrol_id and antibody
        .map { ipcontrol_id, antibody, ip_meta, ip_tagalign, ipcontrol_meta, ipcontrol_tagalign ->
            def pooled_ipcontrol_id = ipcontrol_id - ~/_bRep_.*$/
            [ pooled_ipcontrol_id, antibody, ip_meta, ip_tagalign, ipcontrol_meta, ipcontrol_tagalign ]
        }
        .groupTuple(by: [0, 1]) // group samples
        .map { pooled_ipcontrol_id, antibody, ip_metas, ip_tagaligns, ipcontrol_metas, ipcontrol_tagaligns ->
            // if depth ratio between controls is higher than ctl_depth_ratio, then use pooled control
            def ipcontrol_depths = ipcontrol_metas.collect { meta ->
                meta[meta.ref_total_mapped_reads_for_rpm_key]
            }
            def ctl_depth_max = ipcontrol_depths.max()
            def ctl_depth_min = ipcontrol_depths.min()
            def ctl_depth_ratio = ctl_depth_max / ctl_depth_min
            def ctl_depth_ratio_threshold_exceeded = ctl_depth_ratio > ctl_depth_ratio_threshold
            [ pooled_ipcontrol_id, antibody, ctl_depth_max, ctl_depth_min, ctl_depth_ratio, ctl_depth_ratio_threshold_exceeded, ip_metas, ip_tagaligns, ipcontrol_metas, ipcontrol_tagaligns ]
        }
        .transpose()
        .map { pooled_ipcontrol_id, antibody, ctl_depth_max, ctl_depth_min, ctl_depth_ratio, ctl_depth_ratio_threshold_exceeded, ip_meta, ip_tagalign, ipcontrol_meta, ipcontrol_tagalign ->
            def meta_clone = ip_meta.clone()
            if (ctl_depth_ratio_threshold_exceeded) {
                meta_clone.input_control = pooled_ipcontrol_id
            }
            meta_clone.ctl_depth_max = ctl_depth_max
            meta_clone.ctl_depth_min = ctl_depth_min
            meta_clone.ctl_depth_ratio = ctl_depth_ratio
            meta_clone.ctl_depth_ratio_threshold_exceeded = ctl_depth_ratio_threshold_exceeded
            meta_clone.ipcontrol_is_pooled = ctl_depth_ratio_threshold_exceeded ?: false
            [ meta_clone, ip_tagalign]
        }
        .set { ch_tas_reps_and_pseudoreps_ips_with_ipcontrol }

    // TODO: save for debugging
    ch_tas_reps_and_pseudoreps_ips_with_ipcontrol
        .map { meta, ip_tagalign ->
            "${meta}\t${ip_tagalign}"
        }
        .collectFile(name: 'ch_tas_reps_and_pseudoreps_ips_with_ipcontrol.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_ENCODE_PIPELINE")

    // We remove id and antibody (used when branching above) to mix below
    ch_tas_reps_and_pseudoreps_by_type
        .ipcontrols
        .map { id, antibody, meta, tagalign ->
            [ meta, tagalign ]
        }
        .set { ch_tas_reps_and_pseudoreps_ipcontrols }

    // We mix back the rest with the ips now with updated pooled/non-pooled control metadata
    ch_tas_reps_and_pseudoreps_by_type.ips_wo_ipcontrol
        .mix(ch_tas_reps_and_pseudoreps_pooled)
        .mix(ch_tas_reps_and_pseudoreps_ips_with_ipcontrol)
        .mix(ch_tas_reps_and_pseudoreps_ipcontrols)
        .set { ch_tagalign }

    // Create channel: [ meta, ip_tagalign, ipcontrol_tagalign ]
    ch_tagalign
        .branch { meta, tagalign ->
            ips_with_ipcontrol: meta.input_control
                return [meta.input_control, meta.antibody, meta, tagalign]
            ips_wo_ipcontrol: !meta.input_control && !meta.is_input_control
                return [meta, tagalign]
            ipcontrols: !meta.input_control && meta.is_input_control
                return [meta.id, meta.input_control_of_antibody, meta, tagalign]
        }
        .set { ch_tagalign_by_type }

    ch_tagalign_by_type
        .ips_with_ipcontrol
        .combine(ch_tagalign_by_type.ipcontrols, by: [0, 1])
        .map { ipcontrol_id, antibody, ip_meta, ip_tagalign, ipcontrol_meta, ipcontrol_tagalign ->
            [ ip_meta, ip_tagalign, ipcontrol_tagalign ]
        }
        .mix(ch_tagalign_by_type.ips_wo_ipcontrol)
        .map { it -> [ it[0], it[1], it[2] ?: []] }
        .set { ch_tagalign_for_spp }

    // TODO: save for debugging
    ch_tagalign_for_spp
        .map { meta, ip_tagalign, ipcontrol_tagalign ->
            "${meta}\t${ip_tagalign}\t${ipcontrol_tagalign}"
        }
        .collectFile(name: 'ch_tagalign_for_spp.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_ENCODE_PIPELINE")

    //
    // MODULE: Call peaks with phantompeakqualtools SPP
    //
    PHANTOMPEAKQUALTOOLS_SPP (
        ch_tagalign_for_spp
    )

    //
    // MODULE: Filter peaks by blacklist, chromosomes, and max score
    //
    PEAKS_FILTER_BLACKLIST (
        PHANTOMPEAKQUALTOOLS_SPP.out.regionpeak,
        ch_blacklist,
        true, // filter_chr
        val_peak_type,
        encode_peak_max_score
    )

    // Create channel: [ meta, [peaks1, peaks2], pooled_peaks ]
    PEAKS_FILTER_BLACKLIST
        .out
        .peaks
        .branch { meta, peak ->
            true_replicates: !meta.is_pseudoreplicate && !meta.is_pooled
                return [ meta.id, meta, peak ]
            pooled_replicates: meta.is_pooled && !meta.is_pseudoreplicate
                return [ meta.id, peak ]
            pseudoreplicates: meta.is_pseudoreplicate && !meta.is_pooled
                return [ meta.id, meta, peak ]
            pooled_pseudoreplicates: meta.is_pseudoreplicate && meta.is_pooled
                return [ meta.id, meta, peak ]
        }
        .set { ch_spp_peaks_by_type }

    // True replicates
    ch_spp_peaks_by_type
        .true_replicates
        .map { id, meta, peak ->
            def pooled_id = id - ~/_bRep_.*$/
            [ pooled_id, meta, peak ]
        }
        .set { ch_spp_peaks_true_reps }

    ch_spp_peaks_true_reps
        .combine(ch_spp_peaks_true_reps, by: 0)
        // Ensure unordered combinations without repetition
        // See https://github.com/nextflow-io/nextflow/discussions/2109#discussioncomment-12501996
        .filter { id, meta1, peak1, meta2, peak2 -> meta1.id < meta2.id }
        .map{ id, meta1, peak1, meta2, peak2 ->
            def meta_clone = meta1.clone()
            meta_clone.id = id
            meta_clone.peak_consensus_pair_type = 'true_replicate'
            meta_clone.idr_pair_breps = [meta1.brep, meta2.brep]
            [ meta_clone.id, meta_clone, [ peak1, peak2 ] ]
        }
        .combine(ch_spp_peaks_by_type.pooled_replicates, by: 0)
        .map { id, meta, peaks, pooled_peak ->
            [ meta, peaks, pooled_peak ]
        }
        .set { ch_spp_peaks_true_reps_for_idr }

    // TODO: save for debugging
    ch_spp_peaks_true_reps_for_idr
        .map { meta, peaks, pooled_peak ->
            "${meta}\t${peaks}\t${pooled_peak}"
        }
        .collectFile(name: 'ch_spp_peaks_true_reps_for_idr.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_ENCODE_PIPELINE")

    // Pooled pseudoreplicates
    ch_spp_peaks_by_type
        .pooled_pseudoreplicates
        .groupTuple(by: 0)
        .combine(ch_spp_peaks_by_type.pooled_replicates, by: 0)
        .map { id, metas, peaks, pooled_peak ->
            // Sort metas and peaks to ensure consistent ordering for caching and resuming
            def sorted_peaks = peaks.sort { peak -> peak.name }
            def sorted_metas = metas.sort { meta -> meta.brep }
            def meta_clone = sorted_metas[0].clone()
            meta_clone.peak_consensus_pair_type = 'pooled_pseudoreplicate'
            [ meta_clone, sorted_peaks, pooled_peak ]
        }
        .set { ch_spp_peaks_pooled_pseudoreps_for_idr }

    // TODO: save for debugging
    ch_spp_peaks_pooled_pseudoreps_for_idr
        .map { meta, peaks, pooled_peak ->
            "${meta}\t${peaks}\t${pooled_peak}"
        }
        .collectFile(name: 'ch_spp_peaks_pooled_pseudoreps_for_idr.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_ENCODE_PIPELINE")

    // Self-pseudoreplicates
    ch_spp_peaks_by_type
        .pseudoreplicates
        .groupTuple(by: 0)
        // id, metas, peaks
        .combine(ch_spp_peaks_by_type.true_replicates, by: 0)
        .map { id, metas, peaks, true_rep_meta, true_rep_peak ->
            // Sort metas and peaks to ensure consistent ordering for caching and resuming
            def sorted_peaks = peaks.sort { peak -> peak.name }
            def sorted_metas = metas.sort { meta -> meta.pseudoreplicate }
            def meta_clone = sorted_metas[0].clone()
            meta_clone.remove('pseudoreplicate')
            meta_clone.id = id
            meta_clone.peak_consensus_pair_type = 'self_pseudoreplicate'
            [ meta_clone, sorted_peaks, true_rep_peak ]
        }
        .set { ch_spp_peaks_self_pseudoreps_for_idr }

    // TODO: save for debugging
    ch_spp_peaks_self_pseudoreps_for_idr
        .map { meta, peaks, true_rep_peak ->
            "${meta}\t${peaks}\t${true_rep_peak}"
        }
        .collectFile(name: 'ch_spp_peaks_self_pseudoreps_for_idr.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_ENCODE_PIPELINE")


    // Use IDR to compare all pairs of matched replicates
    // (1) True replicates narrowPeak files: ${REP1_PEAK_FILE} vs. ${REP2_PEAK_FILE} IDR results transferred to Pooled-replicates narrowPeak file  ${POOLED_PEAK_FILE}
    // (2) Pooled-pseudoreplicates: ${PPR1_PEAK_FILE} vs. ${PPR2_PEAK_FILE} IDR results transferred to Pooled-replicates narrowPeak file ${POOLED_PEAK_FILE}
    // (3) Rep1 self-pseudoreplicates: ${REP1_PR1_PEAK_FILE} vs. ${REP1_PR2_PEAK_FILE} IDR results transferred to Rep1 narrowPeak file ${REP1_PEAK_FILE}
    // (4) Rep2 self-pseudoreplicates: ${REP2_PR1_PEAK_FILE} vs. ${REP2_PR2_PEAK_FILE} IDR results transferred to Rep2 narrowPeak file ${REP2_PEAK_FILE}
    ch_spp_peaks_true_reps_for_idr
        .mix(ch_spp_peaks_pooled_pseudoreps_for_idr)
        .mix(ch_spp_peaks_self_pseudoreps_for_idr)
        .map { meta, peaks, pooled_peak ->
            [ meta + [ peak_consensus_type: 'idr' ], peaks, val_peak_type, pooled_peak ]
        }
        .set { ch_for_idr }


    // TODO: save for debugging
    ch_for_idr
        .map { meta, peaks, peak_type, pooled_peak ->
            "${meta}\t${peaks}\t${peak_type}\t${pooled_peak}"
        }
        .collectFile(name: 'ch_for_idr.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_ENCODE_PIPELINE")

    //
    // MODULE: IDR analysis
    //
    IDR (
        ch_for_idr
    )

    //
    // MODULE: Filter peaks by IDR threshold
    //
    IDR_FILTER_THRESHOLD (
        IDR.out.idr,
        val_peak_type,
        idr_filtering_threshold
    )


    ch_for_idr
        .map { meta, peaks, peak_type, pooled_peak ->
            [ meta + [ peak_consensus_type: 'naive_overlap' ], peaks, pooled_peak ]
        }
        .set { ch_for_naive_overlap }

    //
    // MODULE: Naive overlap thresholding as an alternative to IDR for histone marks
    // See Section 6) in https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit?tab=t.0#heading=h.9ecc41kilcvq
    //
    PEAKS_NAIVE_OVERLAP (
        ch_for_naive_overlap,
        val_peak_type
    )

    IDR_FILTER_THRESHOLD
        .out
        .peaks
        .mix(PEAKS_NAIVE_OVERLAP.out.peak_overlap)
        .set { ch_peaks_for_fltbl }

    //
    // MODULE: Filter peaks by blacklist
    //
    CONSENSUS_FILTER_BLACKLIST (
            ch_peaks_for_fltbl,
            ch_blacklist,
            true, // filter_chr
            val_peak_type,
            encode_peak_max_score
    )
    ch_peaks_fltbl = CONSENSUS_FILTER_BLACKLIST.out.peaks

    // Create channel: [ meta, idr_true_reps, idr_pseudo_reps, idr_pooled_pseudoreps ]
    //ch_peaks_fltbl
    //     .map { meta, peak ->
    //         // Remove the bRep from pseudoreplicates' ids
    //         def pooled_id = meta.id - ~/_bRep_.*$/
    //         [ pooled_id, meta, peak ]
    //     }
    //     .branch { pooled_id, meta, peak ->
    //         true_replicates: meta.peak_consensus_pair_type == 'true_replicate'
    //             return [ pooled_id, meta.antibody, meta, peak ]
    //         self_pseudoreplicates: meta.peak_consensus_pair_type == 'self_pseudoreplicate'
    //             return [ pooled_id, meta.antibody, meta, peak ]
    //         pooled_pseudoreplicates: meta.peak_consensus_pair_type == 'pooled_pseudoreplicate'
    //             return [ pooled_id, meta.antibody, peak ]
    //     }
    //     .set { ch_idr_peaks }

    // ch_idr_peaks
    //     .self_pseudoreplicates
    //     .groupTuple(by: [0, 1])
    //     .map { pooled_id, antibody, metas, peaks ->
    //         def sorted_peaks = peaks.sort { peak -> peak.name }
    //         [ pooled_id, antibody, peaks ]
    //     }
    //     .set { ch_idr_peaks_self_pseudoreps }

    // ch_idr_peaks
    //     .true_replicates
    //     .combine(ch_idr_peaks_self_pseudoreps, by: [0,1])
    //     .combine(ch_idr_peaks.pooled_pseudoreplicates, by: [0,1])
    //     .map { pooled_id, antibody, true_rep_meta, true_rep_peak, self_pseudo_rep_peaks, pooled_pseudo_rep_peaks ->
    //         [ true_rep_meta, true_rep_peak, self_pseudo_rep_peaks, pooled_pseudo_rep_peaks ]
    //     }
    //     .set { ch_idr_peaks_for_qc }

    // // TODO: save for debugging
    // ch_idr_peaks_for_qc
    //     .map { meta, true_rep_peaks, self_pseudo_rep_peaks, pooled_pseudo_rep_peaks ->
    //         "${meta}\t${true_rep_peaks}\t${self_pseudo_rep_peaks}\t${pooled_pseudo_rep_peaks}"
    //     }
    //     .collectFile(name: 'ch_idr_peaks_for_qc.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_ENCODE_PIPELINE")

    //
    // MODULE: Compute IDR QC scores
    //
    // PEAKS_IDR_QC (
    //     ch_idr_peaks_for_qc
    // )




    // Create channel: [ meta, tagalign, ccscores, idr_peaks ]

    // Join tagalign and ccscores
    PHANTOMPEAKQUALTOOLS_SPP
        .out
        .ccscores
        .filter { meta, ccscores -> !meta.is_pseudoreplicate }
        .map { meta, ccscores -> [ meta.id, meta.antibody, ccscores ] }
        .set { ch_ccscores}

    ch_tagalign
        .filter {meta, tagalign -> !meta.is_pseudoreplicate }
        .map { meta, tagalign -> [ meta.id, meta.antibody, meta, tagalign ] }
        .combine(ch_ccscores, by: [0, 1])
        // [ id, antibody, meta, tagalign, ccscores ]
        .set { ch_ta_ccscores }

    // We want to compute FRiP scores for the following comparisons:
    // Rep1 tagAlign vs. IDR/overlap peak from pseudo replicates of Rep1 with estimated fragment length of Rep1
    // Rep2 tagAlign vs. IDR/overlap peak from pseudo replicates of Rep2 with estimated fragment length of Rep2
    // Pooled tagAlign vs. IDR/overlap peak from true replicates (Nt) with mean estimated fragment length of Rep1 and Rep2
    // Pooled tagAlign vs. IDR/overlap peak from pooled pseudo replicates (Np) with mean estimated fragment length of Rep1 and Rep2
    ch_peaks_fltbl
        .map { meta, peak ->
            [ meta.id, meta.antibody, meta, peak ]
        }
        .set { ch_peaks_fltbl }

    ch_ta_ccscores
        .combine(ch_peaks_fltbl, by: [0, 1])
        .map { id, antibody, meta_tagalign, tagalign, ccscores, meta_peak, peak ->
            def meta_clone = meta_tagalign.clone()
            meta_clone.peak_consensus_type = meta_peak.peak_consensus_type
            meta_clone.peak_consensus_pair_type = meta_peak.peak_consensus_pair_type
            if (meta_clone.peak_consensus_pair_type == 'true_replicate') {
                meta_clone.idr_pair_breps = meta_peak.idr_pair_breps
            }
            [ meta_clone, tagalign, ccscores, peak ]
        }
        .set { ch_ta_ccscores_peaks }

    // TODO: save for debugging
    ch_ta_ccscores_peaks
        .map { meta, tagalign, ccscores, peak ->
            "${meta}\t${tagalign}\t${ccscores}\t${peak}"
        }
        .collectFile(name: 'ch_ta_ccscores_peaks.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_ENCODE_PIPELINE")


    //
    // MODULE: Compute FRiP scores for IDR and naive overlap peaks
    //
    TAGALIGN_FRIP_SCORE (
        ch_ta_ccscores_peaks,
        ch_chromsizes
    )

    emit:
    tagalign            = ch_tagalign                                // channel: [ val(meta), path(tagalign) ]
    ccscores            = PHANTOMPEAKQUALTOOLS_SPP.out.ccscores       // channel: [ val(meta), path(ccscores) ]
    spp_peaks           = PEAKS_FILTER_BLACKLIST.out.peaks            // channel: [ val(meta), path(peak) ]
    idr                 = IDR.out.idr                                // channel: [ val(meta), path(idr) ]
    idr_peaks           = IDR_FILTER_THRESHOLD.out.peaks              // channel: [ val(meta), path(peak) ]
    naive_overlap_peaks = PEAKS_NAIVE_OVERLAP.out.peak_overlap        // channel: [ val(meta), path(peak) ]
    consensus_peaks     = CONSENSUS_FILTER_BLACKLIST.out.peaks        // channel: [ val(meta), path(peak) ]
    frip_bed            = TAGALIGN_FRIP_SCORE.out.bed                 // channel: [ val(meta), path(bed) ]
    frip                = TAGALIGN_FRIP_SCORE.out.frip                // channel: [ val(meta), path(frip_txt) ]
}
