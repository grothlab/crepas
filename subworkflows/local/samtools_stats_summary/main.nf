//
// Create a summary of samtools stat tables
//

include { SAMTOOLS_STATS_TRANSPOSE             } from '../../../modules/local/samtools/stats_transpose/main'
include { SAMTOOLS_STATS_CAT                  } from '../../../modules/local/samtools/stats_cat/main'
include { FINAL_STAT_SUMMARY } from '../../../modules/local/samtools/final_stat_summary/main'

workflow SAMTOOLS_STATS_SUMMARY {

    take:
    ch_stats
    genome
    spikein_genome

    main:

    ch_versions = Channel.empty()

    SAMTOOLS_STATS_TRANSPOSE (
        ch_stats
    )
    ch_col_stats = SAMTOOLS_STATS_TRANSPOSE.out.t_stats.collect{ it[1] }
    ch_versions = ch_versions.mix(SAMTOOLS_STATS_TRANSPOSE.out.versions.first())

    SAMTOOLS_STATS_CAT (
        ch_col_stats
    )
    ch_versions = ch_versions.mix(SAMTOOLS_STATS_CAT.out.versions.first())

    FINAL_STAT_SUMMARY (
        SAMTOOLS_STATS_CAT.out.cat,
        genome,
        spikein_genome
    )
    ch_versions = ch_versions.mix(FINAL_STAT_SUMMARY.out.versions.first())

    emit:

    summary     = FINAL_STAT_SUMMARY.out.tsv    // channel: [ val(meta), [ tsv ] ]
    
    versions    = ch_versions               // channel: [ versions.yml ]
}