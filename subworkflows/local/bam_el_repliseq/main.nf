include { DEEPTOOLS_MULTIBAMSUMMARY } from '../../../modules/nf-core/deeptools/multibamsummary/main'
include { REPLISEQ_RTNORMALIZE      } from '../../../modules/local/repliseq_rtnormalize/main'
include { FILE_SORT                 } from '../../../modules/local/file_sort/main'
include { UCSC_BEDGRAPHTOBIGWIG     } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'
include { FIND_CONCATENATE as REPLISEQ_RT_MULTIQC } from '../../../modules/nf-core/find/concatenate/main'
include { FIND_CONCATENATE as REPLISEQ_GENE_CLASS_MULTIQC } from '../../../modules/nf-core/find/concatenate/main'
include { REPLISEQ_RT_DOMAINS      } from '../../../modules/local/repliseq_rt_domains/main'
include { BED_TO_SAF               } from '../../../modules/local/bed_to_saf/main'
include { SUBREAD_FEATURECOUNTS as SUBREAD_FEATURECOUNTS_GENES } from '../../../modules/nf-core/subread/featurecounts/main'
include { REPLISEQ_CLASSIFY_GENES  } from '../../../modules/local/repliseq_classify_genes/main'
include { rtFractionOrder             } from '../utils_grothlab_crepas_pipeline'

workflow BAM_EL_REPLISEQ {

    take:
    ch_bam_bai          // channel: [ val(meta), path(bam), path(bai) ]
    ch_chrom_sizes_endo // channel: [ val(meta), path(chrom_sizes) ]
    ch_blacklist        // channel: [ val(meta), path(blacklist) ]
    ch_rt_header        // channel: path(repliseq_rt_header.txt)
    ch_gene_class_header // channel: path(repliseq_gene_class_header.txt)
    ch_gene_bed         // channel: [ val(meta), path(gene.bed) ]

    main:

    //
    // Group all replicates (every fraction) of the same condition together, sorted
    // by (rt_fraction, brep) so bams/bais/fractions/breps stay aligned.
    //
    ch_bam_bai
        .map { meta, bam, bai ->
            def meta_clone = meta.clone()
            // The id at this point is "<exp_type>_<sample>_<rt_fraction>_bRep_<n>"
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
            def sorted_bams = replicates.collect { it -> it[1] }
            def sorted_bais = replicates.collect { it -> it[2] }
            def sorted_fractions = sorted_metas.collect { meta -> meta.rt_fraction }
            def sorted_breps = sorted_metas.collect { meta -> meta.brep }
            meta_clone.remove('brep')
            meta_clone.remove('rt_fraction')
            [ meta_clone, sorted_bams, sorted_bais, sorted_fractions, sorted_breps ]
        }
        .set { ch_condition_bam_bai }

    ch_condition_bam_bai
        .map { meta, bams, bais, fractions, breps ->
            def labels = [fractions, breps].transpose().collect { fraction, brep -> "${fraction}_${brep}" }
            [ meta, bams, bais, labels ]
        }
        .set { ch_for_multibamsummary }

    //
    // MODULE: Count reads per genomic window (bins mode)
    // across all replicates of a condition
    //
    DEEPTOOLS_MULTIBAMSUMMARY (
        ch_for_multibamsummary,
        ch_blacklist
    )

    ch_fractions_breps = ch_condition_bam_bai.map { meta, _bams, _bais, fractions, breps -> [ meta, fractions, breps ] }


    DEEPTOOLS_MULTIBAMSUMMARY.out.raw_counts
        .join(ch_fractions_breps, by: 0)
        .map { meta, counts, fractions, breps ->
            [ meta, counts, fractions, breps ]
        }
        .set { ch_counts_for_rtnormalize }

    //
    // MODULE: Calculate and normalize the RT track: CPM, paired/pooled/auto replicate
    // combination, optional quantile normalization between replicates, loess/roll smoothing
    //
    REPLISEQ_RTNORMALIZE (
        ch_counts_for_rtnormalize
    )

    //
    // MODULE: Prepend the MultiQC custom-content header to each sample's summary row.
    //
    REPLISEQ_RTNORMALIZE.out.summary
        .combine(ch_rt_header)
        .map { meta, summary, header -> [ meta, [ header, summary ] ] }
        .set { ch_rt_mqc }

    REPLISEQ_RT_MULTIQC (
        ch_rt_mqc
    )

    ch_rt_raw = REPLISEQ_RTNORMALIZE.out.raw.map { meta, bdg ->
        def meta_clone = meta.clone()
        meta_clone.rt_measure = 'ratio'
        meta_clone.rt_track_type = 'raw'
        [ meta_clone, bdg ]
    }

    ch_rt_smooth = REPLISEQ_RTNORMALIZE.out.smooth.map { meta, bdg ->
        def meta_clone = meta.clone()
        meta_clone.rt_measure = 'ratio'
        meta_clone.rt_track_type = 'smooth'
        [ meta_clone, bdg ]
    }

    // Only emitted when a third fraction was supplied, so these channels are empty when only early/late are inputted
    ch_rt_index_raw = REPLISEQ_RTNORMALIZE.out.rt_index_raw.map { meta, bdg ->
        def meta_clone = meta.clone()
        meta_clone.rt_measure = 'rt_index'
        meta_clone.rt_track_type = 'raw'
        [ meta_clone, bdg ]
    }

    ch_rt_index_smooth = REPLISEQ_RTNORMALIZE.out.rt_index_smooth.map { meta, bdg ->
        def meta_clone = meta.clone()
        meta_clone.rt_measure = 'rt_index'
        meta_clone.rt_track_type = 'smooth'
        [ meta_clone, bdg ]
    }

    ch_rt_raw_covered = REPLISEQ_RTNORMALIZE.out.raw_covered.map { meta, bdg ->
        def meta_clone = meta.clone()
        meta_clone.rt_measure = 'ratio'
        meta_clone.rt_track_type = 'raw'
        meta_clone.rt_covered = true
        [ meta_clone, bdg ]
    }

    ch_rt_smooth_covered = REPLISEQ_RTNORMALIZE.out.smooth_covered.map { meta, bdg ->
        def meta_clone = meta.clone()
        meta_clone.rt_measure = 'ratio'
        meta_clone.rt_track_type = 'smooth'
        meta_clone.rt_covered = true
        [ meta_clone, bdg ]
    }

    // One smoothed coverage bedGraph per input replicate (only emitted under the default --smooth_stage replicate)
    ch_rt_replicate_smooth = REPLISEQ_RTNORMALIZE.out.replicate_smooth.flatMap { meta, bdg ->
        bdg.collect { f ->
            def replicate_id = (f.name =~ /\.([^.]+)\.coverage\.smooth\.bedGraph$/)[0][1]
            def meta_clone = meta.clone()
            meta_clone.id = "${meta.id}.${replicate_id}"
            meta_clone.rt_measure = 'coverage'
            meta_clone.rt_track_type = 'smooth'
            meta_clone.rt_replicate = replicate_id
            [ meta_clone, f ]
        }
    }

    ch_domains = channel.empty()
    ch_domain_qc = channel.empty()
    ch_domain_box = channel.empty()
    ch_gene_classes = channel.empty()
    ch_gene_class_qc = channel.empty()
    ch_featurecounts_summary = channel.empty()
    ch_gene_class_mqc = channel.empty()
    ch_gene_class_box = channel.empty()
    ch_gene_class_bed = channel.empty()

    //
    // MODULE: Call domains of constant replication timing. Which track they are called on is a
    // parameter, since smoothing changes where the boundaries fall.
    //
    if (!params.skip_repliseq_domains) {
        ch_rt_raw_covered
            .mix(ch_rt_smooth_covered)
            .filter { meta, _bdg ->
                params.repliseq_domain_track == 'both' || meta.rt_track_type == params.repliseq_domain_track
            }
            .set { ch_for_domains }

        REPLISEQ_RT_DOMAINS (
            ch_for_domains
        )
        ch_domains = REPLISEQ_RT_DOMAINS.out.domains
        ch_domain_qc = REPLISEQ_RT_DOMAINS.out.qc
        ch_domain_box = REPLISEQ_RT_DOMAINS.out.mqc_box
    }

    //
    // MODULE: Classify each gene by the S-phase fraction with the highest read density over its
    // gene body.
    //
    if (!params.skip_repliseq_gene_classification) {
        BED_TO_SAF (
            ch_gene_bed
        )

        ch_condition_bam_bai
            .combine(BED_TO_SAF.out.saf.map { it -> it[1] })
            .map { meta, bams, _bais, _fractions, _breps, saf -> [ meta, bams, saf ] }
            .set { ch_genes_for_featurecounts }

        SUBREAD_FEATURECOUNTS_GENES (
            ch_genes_for_featurecounts
        )

        SUBREAD_FEATURECOUNTS_GENES.out.counts
            .join(ch_fractions_breps, by: 0)
            .set { ch_gene_counts_for_classify }

        REPLISEQ_CLASSIFY_GENES (
            ch_gene_counts_for_classify
        )
        ch_gene_classes = REPLISEQ_CLASSIFY_GENES.out.classes
        ch_gene_class_qc = REPLISEQ_CLASSIFY_GENES.out.qc
        ch_featurecounts_summary = SUBREAD_FEATURECOUNTS_GENES.out.summary
        ch_gene_class_box = REPLISEQ_CLASSIFY_GENES.out.mqc_box

        ch_gene_class_bed = REPLISEQ_CLASSIFY_GENES.out.early.map { meta, bed -> [ meta, bed, 'early' ] }
            .mix(REPLISEQ_CLASSIFY_GENES.out.mid.map { meta, bed -> [ meta, bed, 'mid' ] })
            .mix(REPLISEQ_CLASSIFY_GENES.out.late.map { meta, bed -> [ meta, bed, 'late' ] })
            .mix(REPLISEQ_CLASSIFY_GENES.out.unclassified.map { meta, bed -> [ meta, bed, 'unclassified' ] })
            .map { meta, bed, cls ->
                def meta_clone = meta.clone()
                meta_clone.rt_gene_class = cls
                [ meta_clone, bed ]
            }

        // Same filename-order constraint as the track summary above:
        // "repliseq_gene_class_header.txt" sorts before "rt_gene_class_counts.tsv".
        REPLISEQ_CLASSIFY_GENES.out.summary
            .combine(ch_gene_class_header)
            .map { meta, summary, header -> [ meta, [ header, summary ] ] }
            .set { ch_gene_class_for_mqc }

        REPLISEQ_GENE_CLASS_MULTIQC (
            ch_gene_class_for_mqc
        )
        ch_gene_class_mqc = REPLISEQ_GENE_CLASS_MULTIQC.out.file_out
    }

    //
    // MODULE: Lexicographically sort the RT tracks (required by UCSC_BEDGRAPHTOBIGWIG)
    //
    FILE_SORT (
        ch_rt_raw.mix(ch_rt_smooth).mix(ch_rt_raw_covered).mix(ch_rt_smooth_covered).mix(ch_rt_index_raw).mix(ch_rt_index_smooth).mix(ch_rt_replicate_smooth),
        'bedGraph'
    )

    //
    // MODULE: Convert RT tracks (raw + smoothed, plus per-replicate smoothed coverage) to bigWig
    //
    UCSC_BEDGRAPHTOBIGWIG (
        FILE_SORT.out.sorted,
        ch_chrom_sizes_endo.map { it -> it[1] }
    )

    emit:
    bedgraph = FILE_SORT.out.sorted             // channel: [ val(meta), path(bedgraph) ]
    bigwig   = UCSC_BEDGRAPHTOBIGWIG.out.bigwig // channel: [ val(meta), path(bigwig) ]
    qc       = REPLISEQ_RTNORMALIZE.out.qc      // channel: [ val(meta), path(qc.txt) ]
    domains  = ch_domains                       // channel: [ val(meta), path(RT_domains.bed) ]
    domain_qc = ch_domain_qc                    // channel: [ val(meta), path(RT_domains.qc.txt) ]
    domain_box = ch_domain_box                  // channel: [ val(meta), path(*_mqc.json) ]
    gene_classes = ch_gene_classes              // channel: [ val(meta), path(gene_RT_class.tsv) ]
    gene_class_bed = ch_gene_class_bed          // channel: [ val(meta), path(gene_RT_class.<class>.bed) ]
    gene_class_qc = ch_gene_class_qc            // channel: [ val(meta), path(gene_RT_class.qc.txt) ]
    featurecounts_summary = ch_featurecounts_summary // channel: [ val(meta), path(summary) ]
    gene_class_mqc = ch_gene_class_mqc          // channel: [ val(meta), path(gene_class_mqc.tsv) ]
    gene_class_box = ch_gene_class_box          // channel: [ val(meta), path(*_mqc.json) ]
    mqc      = REPLISEQ_RT_MULTIQC.out.file_out // channel: [ val(meta), path(rt_mqc.tsv) ]
}
