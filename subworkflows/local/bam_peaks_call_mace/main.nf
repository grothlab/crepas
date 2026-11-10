//
// Call peaks with MACE (for ChIP-exo samples)
//
include { MACE_PREPROCESSOR        } from '../../../modules/local/mace/preprocessor/main'
include { UCSC_WIGTOBIGWIG         } from '../../../modules/nf-core/ucsc/wigtobigwig/main'
include { MACE_MACE                } from '../../../modules/local/mace/mace/main'

workflow BAM_PEAKS_CALL_MACE {
    take:
    ch_bam_bai                        // channel: [ val(meta), bam, bai ]
    ch_chrom_sizes                    // channel: [ val(meta), path(chrom_sizes) ]

    main:

    // Create channel: [ meta, [bams_merged_reps] ]
    // MACE has no input control option, so input controls are left out rather than treated as samples
    ch_bam_bai
        .filter { meta, _bam, _bai -> !meta.is_input_control }
        .map { meta, bam, bai ->
            def meta_clone = meta.clone()
            def antibody = meta_clone.antibody
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
        .collectFile(name: 'ch_bam_merged_reps.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_PEAKS_CALL_MACE")

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

    //
    // MODULE: Convert the strand-specific wig files to bigWig
    //
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
    // Peaks are called on the pooled replicates, so join with the channel MACE was given
    ch_bam_merged_reps
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
        .set { ch_bam_peaks }

    emit:
    peaks                        = ch_mace_peaks                    // channel: [ val(meta), path(peaks) ]
    bam_peaks                    = ch_bam_peaks                     // channel: [ val(meta), path(bam), path(peaks) ]
    qc_peaks                     = ch_mace_peaks                    // channel: [ val(meta), path(peaks) ]
    border                       = MACE_MACE.out.border             // channel: [ val(meta), path(bed) ]
    border_cluster               = MACE_MACE.out.border_cluster     // channel: [ val(meta), path(bed) ]
    border_pair_elite            = MACE_MACE.out.border_pair_elite  // channel: [ val(meta), path(bed) ]
    border_pair                  = MACE_MACE.out.border_pair        // channel: [ val(meta), path(bed) ]
    border_pair_peak             = MACE_MACE.out.border_pair_peak   // channel: [ val(meta), path(bed) ]
}
