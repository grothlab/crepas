//
// UMI-tools dedup, index BAM file and run samtools stats, flagstat and idxstats
//
include { BAM_SPLIT_BY_CHROMOSOME } from '../../../modules/local/bam_split_by_chromosome/main'
include { UMITOOLS_DEDUP     } from '../../../modules/nf-core/umitools/dedup/main'
include { SAMTOOLS_MERGE      } from '../../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS } from '../../../subworkflows/nf-core/bam_stats_samtools/main'

workflow BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS {
    take:
    ch_bam_index          // channel: [ val(meta), path(bam), path(index/csi) ]
    ch_chrom_sizes       // channel: [ val(meta), path(chrom_sizes) ]
    val_get_dedup_stats // boolean: true/false
    skip_split_by_chrom // boolean: true/false

    main:

    if (!skip_split_by_chrom) {
        // creating a channel with each chromosome to iterate over
        ch_chrom_sizes
            .map {
                meta, bed ->
                    bed.splitCsv(header:false, sep:'\t')
            }
            .flatMap { chrom_list ->
            chrom_list.collect { chroms -> chroms[0] }
            }
            .set { ch_chroms }

            // get only the chromosomes that do not contain "_" or "."
            ch_endo_chroms = ch_chroms
                .filter { chrom -> !chrom.contains('_') && !chrom.contains('.') }
                .map { chrom -> [ chrom, [ chrom ] ] }

            ch_exo_and_scaff = ch_chroms
                .filter { chrom -> chrom.contains('_') || chrom.contains('.') }
                .toList()
                .map { exo_and_scaff_list -> [ 'exo_and_scaff', exo_and_scaff_list ] }

            ch_chroms = ch_endo_chroms.mix(ch_exo_and_scaff)
            ch_split_count = ch_chroms.collect { it -> it[0]}.size()

        // print ch_chroms to file for debugging
        // ch_chroms
        //     .map {
        //         split_id, chrom ->
        //             "${split_id}\t${chrom}"
        //     }
        //     .collectFile( name: 'ch_chroms.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS" )

        // print ch_split_count to file for debugging
        // ch_split_count
        //     .map { split_count -> "split_count\t${split_count}" }
        //     .collectFile( name: 'ch_split_count.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS" )


        ch_bam_index
            .combine(ch_chroms)
            .combine(ch_split_count)
            .map {
                meta, bam, index, split_id, chrom_list, split_count ->
                    def meta_clone = meta.clone()
                    meta_clone.split_id = split_id
                    meta_clone.chrom_list = chrom_list
                    meta_clone.split_count = split_count
                    [ meta_clone, bam, index, chrom_list ]
            }
            .set { ch_bam_index_chroms }

        // print for debugging
        ch_bam_index_chroms
            .map {
                meta, bam, index, chrom ->
                    "${meta}\t${bam}\t${index}\t${chrom}"
            }
            .collectFile( name: 'ch_bam_index_chroms.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS" )

        // Split BAMs by chromosome
        BAM_SPLIT_BY_CHROMOSOME (
            ch_bam_index_chroms,
            'bai' // index_format
        )
        ch_bam_index = BAM_SPLIT_BY_CHROMOSOME.out.bam.join(BAM_SPLIT_BY_CHROMOSOME.out.bai, by: [0])
    }

    //
    // UMI-tools dedup
    //
    UMITOOLS_DEDUP ( ch_bam_index, val_get_dedup_stats )
    ch_dedup_bam = UMITOOLS_DEDUP.out.bam

    if (!skip_split_by_chrom) {

        ch_dedup_bam
            .map { meta, bam ->
                def key = groupKey(meta.id, meta.split_count)
                [ key, meta, bam ]
            }
            .tap { ch_dedup_bams_for_merge0 }
            .groupTuple(by: 0)
            .map {
                key, metas, bams ->
                def sorted_metas = metas.sort { meta -> meta.split_id }
                def meta_clone = sorted_metas[0].clone()
                meta_clone.remove('split_id')
                meta_clone.remove('chrom_list')
                def sorted_bams = bams.sort { bam -> bam.name }
                    [ meta_clone, sorted_bams, [] ]
            }
            .set { ch_dedup_bams_for_merge }

        // print for debugging
        ch_dedup_bams_for_merge0
            .map {
                key, meta, bam ->
                    "${key}\t${meta}\t${bam}"
            }
            .collectFile( name: 'ch_dedup_bams_for_merge0.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS" )

        // print for debugging
        ch_dedup_bams_for_merge
            .map {
                meta, bams, indices ->
                    "${meta}\t${bams}\t${indices}"
            }
            .collectFile( name: 'ch_dedup_bams_for_merge.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS" )


        //
        // MODULE: Merge deduplicated BAM files of each chromosome back together
        //
        SAMTOOLS_MERGE (
            ch_dedup_bams_for_merge,
            channel.value([[:], [], [], []])
        )
        ch_dedup_bam = SAMTOOLS_MERGE.out.bam
        ch_dedup_index = SAMTOOLS_MERGE.out.index

    } else {
        //
        // Index BAM file and run samtools stats, flagstat and idxstats
        //
        SAMTOOLS_INDEX ( UMITOOLS_DEDUP.out.bam )
        ch_dedup_index = SAMTOOLS_INDEX.out.index
    }

    ch_bam_index_dedup = ch_dedup_bam.join(ch_dedup_index, by: [0])

    BAM_STATS_SAMTOOLS ( ch_bam_index_dedup, [ [:], [], [] ] )

    emit:
    bam                  = ch_dedup_bam                             // channel: [ val(meta), path(bam) ]
    deduplog             = UMITOOLS_DEDUP.out.log                  // channel: [ val(meta), path(log) ]
    tsv_edit_distance    = UMITOOLS_DEDUP.out.tsv_edit_distance    // channel: [ val(meta), path(tsv) ]
    tsv_per_umi          = UMITOOLS_DEDUP.out.tsv_per_umi          // channel: [ val(meta), path(tsv) ]
    tsv_umi_per_position = UMITOOLS_DEDUP.out.tsv_umi_per_position // channel: [ val(meta), path(tsv) ]

    index                = ch_dedup_index               // channel: [ val(meta), path(index) ]
    stats                = BAM_STATS_SAMTOOLS.out.stats            // channel: [ val(meta), path(stats) ]
    flagstat             = BAM_STATS_SAMTOOLS.out.flagstat         // channel: [ val(meta), path(flagstat) ]
    idxstats             = BAM_STATS_SAMTOOLS.out.idxstats         // channel: [ val(meta), path(idxstats) ]
}
