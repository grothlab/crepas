//
// Splitting input BAMs by spike-in genome (into endogenous and exogenous BAMs)
//

include { BAM_SPLIT_BY_GENOME as BAM_SPLIT_BY_GENOME_ENDO } from '../../../modules/local/bam_split_by_genome/main'
include { BAM_SPLIT_BY_GENOME as BAM_SPLIT_BY_GENOME_EXO } from '../../../modules/local/bam_split_by_genome/main'

// TODO: simplify the following (maybe by concatenating ENDO and EXO channel and then splitting them in the end)
include { SAMTOOLS_INDEX } from '../../../modules/nf-core/samtools/index/main'
include { SAMBAMBA_VIEW } from '../../../modules/local/sambamba/view/main'
include { BAM_SORT_STATS_SAMTOOLS  } from '../../../subworkflows/nf-core/bam_sort_stats_samtools/main'

workflow BAM_SPIKEIN_SPLIT {
    take:
    ch_bam               // channel: [ val(meta), [ bam ], [bai] ]
    ch_fasta             // channel: [ val(meta), path(fasta) ]
    ch_bed               // channel: [ val(meta), path(bed) ]
    genome               // val
    spikein_genome       // val

    main:
    ch_versions = Channel.empty()

    // split BAMs by spike-in genome
    BAM_SPLIT_BY_GENOME_ENDO(ch_bam, spikein_genome, genome, true)
    BAM_SPLIT_BY_GENOME_EXO(ch_bam, spikein_genome, spikein_genome, false)

    // add genome as meta field
    ch_bam_endo = BAM_SPLIT_BY_GENOME_ENDO.out.bam.map { [ it[0] + [ genome: genome ], it[1] ] }
    ch_bam_exo = BAM_SPLIT_BY_GENOME_EXO.out.bam.map { [ it[0] + [ genome: spikein_genome ], it[1] ] }

    ch_bam = ch_bam_endo.mix(ch_bam_exo)

    SAMTOOLS_INDEX(ch_bam)
    ch_bam_bai = ch_bam.join(SAMTOOLS_INDEX.out.bai, by: [0])
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    SAMBAMBA_VIEW(ch_bam_bai, ch_bed)

    // TODO: I would need to separate the genome fasta, right now both endo and exo stats
    // are analyzed with the same (main) genome
    // either way, this has no effect on BAM output (https://bioinformatics.stackexchange.com/a/4218)
    BAM_SORT_STATS_SAMTOOLS(SAMBAMBA_VIEW.out.bam, ch_fasta)

    ch_versions = ch_versions.mix(SAMBAMBA_VIEW.out.versions,
                    BAM_SORT_STATS_SAMTOOLS.out.versions)

    emit:

    bam           = BAM_SORT_STATS_SAMTOOLS.out.bam.filter { it[0].genome == genome }            // channel: [ val(meta), [ bam ] ]
    exo_bam       = BAM_SORT_STATS_SAMTOOLS.out.bam.filter { it[0].genome == spikein_genome }    // channel: [ val(meta), [ bam ] ]

    bai           = BAM_SORT_STATS_SAMTOOLS.out.bai.filter { it[0].genome == genome }            // channel: [ val(meta), [ bai ] ]
    exo_bai       = BAM_SORT_STATS_SAMTOOLS.out.bai.filter { it[0].genome == spikein_genome }

    stats         = BAM_SORT_STATS_SAMTOOLS.out.stats.filter { it[0].genome == genome }            // channel: [ val(meta), [ stats ] ]
    exo_stats     = BAM_SORT_STATS_SAMTOOLS.out.stats.filter { it[0].genome == spikein_genome }    // channel: [ val(meta), [ stats ] ]

    flagstat      = BAM_SORT_STATS_SAMTOOLS.out.flagstat.filter { it[0].genome == genome }            // channel: [ val(meta), [ flagstat ] ]
    exo_flagstat  = BAM_SORT_STATS_SAMTOOLS.out.flagstat.filter { it[0].genome == spikein_genome }    // channel: [ val(meta), [ flagstat ] ]

    idxstats      = BAM_SORT_STATS_SAMTOOLS.out.idxstats.filter { it[0].genome == genome }            // channel: [ val(meta), [ idxstats ] ]
    exo_idxstats  = BAM_SORT_STATS_SAMTOOLS.out.idxstats.filter { it[0].genome == spikein_genome }    // channel: [ val(meta), [ idxstats ] ]

    versions = ch_versions                    // channel: [ versions.yml ]
}
