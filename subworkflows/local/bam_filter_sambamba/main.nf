/*
 * Filter BAM file
 */

include { SAMBAMBA_VIEW }            from '../../../modules/local/sambamba/view/main'
include { BAM_SORT_STATS_SAMTOOLS  } from '../../../subworkflows/nf-core/bam_sort_stats_samtools'

workflow BAM_FILTER_SAMBAMBA {
    take:
    ch_bam_index              // channel: [ val(meta), [ bam ], [ index ]]
    ch_bed                    // channel: [ val(meta2), [ bed ] ]
    ch_fasta                  // channel: [ fasta ]

    main:
    ch_versions = Channel.empty()

    SAMBAMBA_VIEW (
        ch_bam_index,
        ch_bed
    )
    ch_versions = ch_versions.mix(SAMBAMBA_VIEW.out.versions)

    BAM_SORT_STATS_SAMTOOLS (
        SAMBAMBA_VIEW.out.bam,
        ch_fasta
    )
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    emit:
    bam      = BAM_SORT_STATS_SAMTOOLS.out.bam      // channel: [ val(meta), [ bam ] ]
    bai      = BAM_SORT_STATS_SAMTOOLS.out.bai      // channel: [ val(meta), [ bai ] ]
    stats    = BAM_SORT_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]

    versions = ch_versions                    // channel: [ versions.yml ]
}
