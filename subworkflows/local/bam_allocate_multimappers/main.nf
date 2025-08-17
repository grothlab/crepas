//
// Allocate multmappers from a BAM file using one of a number of methods
//

include { SAMTOOLS_SORT             } from '../../../modules/nf-core/samtools/sort/main'
include { ALLO                      } from '../../../modules/local/allo/main'
include { MMR                       } from '../../../modules/local/mmr/main'
include { BAM_SORT_STATS_SAMTOOLS   } from '../../../subworkflows/nf-core/bam_sort_stats_samtools/main'

workflow BAM_ALLOCATE_MULTIMAPPERS {

    take:
    ch_bam     // channel: [ val(meta), [bam] ]
    ch_fasta            // channel: [ val(meta), fasta ]
    allocation_method   // string: e.g. 'allo'

    main:

    ch_versions = Channel.empty()

    SAMTOOLS_SORT (
        ch_bam,
        ch_fasta
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    if (allocation_method == 'allo') {
        ALLO (
            SAMTOOLS_SORT.out.bam
        )
        ch_allocated_bam = ALLO.out.bam
        ch_versions = ch_versions.mix(ALLO.out.versions.first())
    }
   
    if (allocation_method == 'mmr') {
        MMR (
            SAMTOOLS_SORT.out.bam
        )
        ch_allocated_bam = MMR.out.bam
        ch_versions = ch_versions.mix(MMR.out.versions.first())
    }

    BAM_SORT_STATS_SAMTOOLS (
        ch_allocated_bam,
        ch_fasta
    )
    ch_allocated_bam = BAM_SORT_STATS_SAMTOOLS.out.bam
    ch_allocated_index = BAM_SORT_STATS_SAMTOOLS.out.bai
    ch_allocated_flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat
    ch_allocated_stats = BAM_SORT_STATS_SAMTOOLS.out.stats
    ch_allocated_idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    emit:
    bam         = ch_allocated_bam          // channel: [ val(meta), [ bam ] ]
    bai         = ch_allocated_index        // channel: [ val(meta), [ index ] ]
    flagstat    = ch_allocated_flagstat     // channel: [ val(meta), [ flagstat ] ]
    stats       = ch_allocated_stats        // channel: [ val(meta), [ stats ] ]
    idxstats    = ch_allocated_idxstats     // channel: [ val(meta), [ idxstats ] ]

    versions    = ch_versions               // channel: [ versions.yml ]
}