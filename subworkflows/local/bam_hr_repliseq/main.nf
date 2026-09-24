include { DEEPTOOLS_MULTIBAMSUMMARY as DEEPTOOLS_MULTIBAMSUMMARY_HR } from '../../../modules/nf-core/deeptools/multibamsummary/main'
include { HR_REPLISEQ_MAKE_ARRAY    } from '../../../modules/local/hr_repliseq_make_array/main'
include { HR_REPLISEQ_CALL_FEATURES } from '../../../modules/local/hr_repliseq_call_features/main'
include { HR_REPLISEQ_PLOT          } from '../../../modules/local/hr_repliseq_plot/main'
include { FIND_CONCATENATE as HR_REPLISEQ_FEATURES_MULTIQC } from '../../../modules/nf-core/find/concatenate/main'
include { rtFractionOrder           } from '../utils_grothlab_crepas_pipeline'

workflow BAM_HR_REPLISEQ {

    take:
    ch_bam_bai      // channel: [ val(meta), path(bam), path(bai) ]
    ch_blacklist    // channel: [ val(meta), path(blacklist) ]
    ch_hr_header    // channel: path(hr_repliseq_features_header.txt)
    ch_el_track     // channel: [ val(condition_id), path(bedgraph) ], may be empty

    main:

    //
    // Group every fraction of the same condition together, in replication order (G1 first, then
    // S1..S16), so that the count columns and the fraction labels stay aligned.
    //
    ch_bam_bai
        .map { meta, bam, bai ->
            def meta_clone = meta.clone()
            meta_clone.id = "${meta_clone.exp_type}_${meta_clone.rt_condition}"
            [ meta_clone.id, meta_clone, bam, bai ]
        }
        .groupTuple(by: 0)
        .map { _id, metas, bams, bais ->
            def replicates = [metas, bams, bais].transpose().sort { a, b ->
                rtFractionOrder(a[0].rt_fraction) <=> rtFractionOrder(b[0].rt_fraction) ?: a[0].brep <=> b[0].brep
            }
            def sorted_metas = replicates.collect { it -> it[0] }
            def meta_clone = sorted_metas[0].clone()
            meta_clone.remove('brep')
            meta_clone.remove('rt_fraction')
            [ meta_clone,
              replicates.collect { it -> it[1] },
              replicates.collect { it -> it[2] },
              sorted_metas.collect { meta -> meta.rt_fraction },
              sorted_metas.collect { meta -> meta.brep } ]
        }
        .set { ch_condition_bam_bai }

    //
    // MODULE: Count reads per genomic window across every fraction of a condition.
    //
    ch_condition_bam_bai
        .map { meta, bams, bais, fractions, breps ->
            def labels = [fractions, breps].transpose().collect { fraction, brep -> "${fraction}_${brep}" }
            [ meta, bams, bais, labels ]
        }
        .set { ch_for_multibamsummary }

    DEEPTOOLS_MULTIBAMSUMMARY_HR (
        ch_for_multibamsummary,
        ch_blacklist
    )

    ch_fractions = ch_condition_bam_bai.map { meta, _bams, _bais, fractions, _breps -> [ meta, fractions ] }

    //
    // MODULE: Build the Gaussian-smoothed, column-scaled Repli-seq array. When a G1 control is
    // among the fractions it is used to correct the S-phase fractions for mappability
    //
    DEEPTOOLS_MULTIBAMSUMMARY_HR.out.raw_counts
        .join(ch_fractions, by: 0)
        .set { ch_for_make_array }

    HR_REPLISEQ_MAKE_ARRAY (
        ch_for_make_array
    )

    //
    // MODULE: Call initiation zones, timing transition regions, breakages, termination sites and
    // late constant-timing regions from the array.
    //
    HR_REPLISEQ_CALL_FEATURES (
        HR_REPLISEQ_MAKE_ARRAY.out.array
    )

    //
    // MODULE: Prepend the MultiQC custom-content header to each sample's feature counts.
    //
    HR_REPLISEQ_CALL_FEATURES.out.summary
        .combine(ch_hr_header)
        .map { meta, summary, header -> [ meta, [ header, summary ] ] }
        .set { ch_hr_mqc }

    HR_REPLISEQ_FEATURES_MULTIQC (
        ch_hr_mqc
    )

    ch_plots = channel.empty()
    ch_mqc_heatmap = channel.empty()
    ch_mqc_sizes = channel.empty()

    //
    // MODULE: Draw the Repli-seq heatmap per chromosome and summarize the features called
    //
    if (!params.skip_hr_repliseq_plots) {
        HR_REPLISEQ_MAKE_ARRAY.out.array
            .join(HR_REPLISEQ_CALL_FEATURES.out.partition, by: 0)
            .join(HR_REPLISEQ_CALL_FEATURES.out.partition_disjoint, by: 0)
            .join(HR_REPLISEQ_CALL_FEATURES.out.partition_flanked, by: 0)
            .join(HR_REPLISEQ_CALL_FEATURES.out.ttr_speed, by: 0)
            .join(HR_REPLISEQ_CALL_FEATURES.out.cluster_rank, by: 0)
            .join(HR_REPLISEQ_CALL_FEATURES.out.iz, by: 0)
            .join(HR_REPLISEQ_CALL_FEATURES.out.ttr, by: 0)
            .join(HR_REPLISEQ_CALL_FEATURES.out.breakage, by: 0)
            .join(HR_REPLISEQ_CALL_FEATURES.out.termination, by: 0)
            .join(HR_REPLISEQ_CALL_FEATURES.out.ctr, by: 0)
            .map { meta, array, partition, partition_disjoint, partition_flanked, speeds, cluster_rank, iz, ttr, breakage, termination, ctr ->
                [ meta.id, meta, array, partition, partition_disjoint, partition_flanked, speeds,
                  cluster_rank, [ iz, ttr, breakage, termination, ctr ] ]
            }
            // Left-join the early/late track of the same condition, if there is one: a sample may
            // carry both designs, but most carry only the high-resolution one.
            .join(ch_el_track, by: 0, remainder: true)
            .filter { it -> it[1] != null }
            .map { _id, meta, array, partition, partition_disjoint, partition_flanked, speeds, cluster_rank, calls, el_track ->
                [ meta, array, partition, partition_disjoint, partition_flanked, speeds,
                  cluster_rank, el_track ?: [], calls ]
            }
            .set { ch_for_plot }

        HR_REPLISEQ_PLOT (
            ch_for_plot
        )
        ch_plots = HR_REPLISEQ_PLOT.out.plots
        ch_mqc_heatmap = HR_REPLISEQ_PLOT.out.mqc_heatmap
        ch_mqc_sizes = HR_REPLISEQ_PLOT.out.mqc_sizes
    }

    emit:
    array       = HR_REPLISEQ_MAKE_ARRAY.out.array         // channel: [ val(meta), path(hr_array.csv) ]
    array_qc    = HR_REPLISEQ_MAKE_ARRAY.out.qc            // channel: [ val(meta), path(hr_array.qc.txt) ]
    iz          = HR_REPLISEQ_CALL_FEATURES.out.iz         // channel: [ val(meta), path(hr_IZ.bed) ]
    ttr         = HR_REPLISEQ_CALL_FEATURES.out.ttr        // channel: [ val(meta), path(hr_TTR.bed) ]
    breakage    = HR_REPLISEQ_CALL_FEATURES.out.breakage   // channel: [ val(meta), path(hr_breakage.bed) ]
    termination = HR_REPLISEQ_CALL_FEATURES.out.termination// channel: [ val(meta), path(hr_termination.bed) ]
    ctr         = HR_REPLISEQ_CALL_FEATURES.out.ctr        // channel: [ val(meta), path(hr_CTR.bed) ]
    partition   = HR_REPLISEQ_CALL_FEATURES.out.partition  // channel: [ val(meta), path(hr_partition.bed) ]
    partition_disjoint = HR_REPLISEQ_CALL_FEATURES.out.partition_disjoint // channel: [ val(meta), path(hr_partition.disjoint.bed) ]
    breakage_disjoint  = HR_REPLISEQ_CALL_FEATURES.out.breakage_disjoint  // channel: [ val(meta), path(hr_breakage.disjoint.bed) ]
    partition_flanked  = HR_REPLISEQ_CALL_FEATURES.out.partition_flanked  // channel: [ val(meta), path(hr_partition.flanked.bed) ]
    breakage_flanked   = HR_REPLISEQ_CALL_FEATURES.out.breakage_flanked   // channel: [ val(meta), path(hr_breakage.flanked.bed) ]
    ttr_speed          = HR_REPLISEQ_CALL_FEATURES.out.ttr_speed          // channel: [ val(meta), path(hr_TTR.speed.tsv) ]
    cluster_rank       = HR_REPLISEQ_CALL_FEATURES.out.cluster_rank       // channel: [ val(meta), path(hr_cluster_rank.tsv) ]
    qc          = HR_REPLISEQ_CALL_FEATURES.out.qc         // channel: [ val(meta), path(hr_features.qc.txt) ]
    plots       = ch_plots                                 // channel: [ val(meta), path(hr_repliseq.plots.pdf) ]
    mqc_heatmap = ch_mqc_heatmap                           // channel: [ val(meta), path(*_mqc.png) ]
    mqc_sizes   = ch_mqc_sizes                             // channel: [ val(meta), path(*_mqc.json) ]
    mqc         = HR_REPLISEQ_FEATURES_MULTIQC.out.file_out// channel: [ val(meta), path(hr_features_mqc.tsv) ]
}
