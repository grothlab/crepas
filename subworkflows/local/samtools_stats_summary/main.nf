//
// Create a summary of samtools stat tables
//

include { STATS_TRANSPOSE     } from '../../../modules/local/stats_transpose/main'
include { STATS_CAT           } from '../../../modules/local/stats_cat/main'
include { STATS_SUMMARY       } from '../../../modules/local/stats_summary/main'

workflow SAMTOOLS_STATS_SUMMARY {

    take:
    ch_stats
    genome
    spikein_genome

    main:

    STATS_TRANSPOSE (
        ch_stats
    )
    ch_col_stats = STATS_TRANSPOSE.out.t_stats.collect{ it -> it[1] }

    STATS_CAT (
        ch_col_stats
    )

    STATS_SUMMARY (
        STATS_CAT.out.cat,
        genome,
        spikein_genome
    )

    emit:

    summary     = STATS_SUMMARY.out.tsv    // channel: [ val(meta), [ tsv ] ]

}
