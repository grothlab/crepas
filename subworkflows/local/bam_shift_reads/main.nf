// TODO: copied from https://github.com/nf-core/atacseq/blob/347d98e820d59074b58719dce5f0deca5af85bb2/subworkflows/local/bam_shift_reads.nf
// add appropiate credits in meta.yaml

include { DEEPTOOLS_ALIGNMENTSIEVE } from '../../../modules/nf-core/deeptools/alignmentsieve'
include { BAM_SORT_STATS_SAMTOOLS } from '../../../subworkflows/nf-core/bam_sort_stats_samtools'

workflow BAM_SHIFT_READS {
    take:
    ch_bam_bai                   // channel: [ val(meta), [ bam ], [bai] ]
    ch_fasta                     // channel: [ val(meta), [ fasta ] ]

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Shift ATAC-seq reads
    //
    DEEPTOOLS_ALIGNMENTSIEVE (
        ch_bam_bai
    )
    ch_versions = ch_versions.mix(DEEPTOOLS_ALIGNMENTSIEVE.out.versions)

    //
    // SUBWORKFLOW: Sort, index and generate stats for the shifted BAM
    //
    BAM_SORT_STATS_SAMTOOLS (
        DEEPTOOLS_ALIGNMENTSIEVE.out.bam,
        ch_fasta
    )
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    emit:
    bam      = BAM_SORT_STATS_SAMTOOLS.out.bam      // channel: [ val(meta), [ bam ] ]
    bai    = BAM_SORT_STATS_SAMTOOLS.out.bai    // channel: [ val(meta), [ bai ] ]
    stats    = BAM_SORT_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
    versions = ch_versions                          // channel: [ versions.yml ]
}