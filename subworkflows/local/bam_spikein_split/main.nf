//
// Splitting input BAMs by spike-in genome (into endogenous and exogenous BAMs)
//
include { BAM_SPLIT_BY_GENOME as BAM_SPLIT_BY_GENOME_ENDO   } from '../../../modules/local/bam_split_by_genome/main'
include { BAM_SPLIT_BY_GENOME as BAM_SPLIT_BY_GENOME_EXO    } from '../../../modules/local/bam_split_by_genome/main'
include { SAMTOOLS_SORT as SAMTOOLS_NSORT                   } from '../../../modules/nf-core/samtools/sort/main'
include { BAM_REMOVE_ORPHANS                                } from '../../../modules/local/bam_remove_orphans/main'
include { BAM_SORT_STATS_SAMTOOLS                           } from '../../../subworkflows/nf-core/bam_sort_stats_samtools/main'
include { BAM_FLAGSTAT_MAPPED                               } from '../../../modules/local/bam_flagstat_mapped/main'

workflow BAM_SPIKEIN_SPLIT {
    take:
    ch_bam                  // channel: [ val(meta), [ bam ], [index] ]
    ch_fasta_fai            // channel: [ val(meta), path(fasta), path(fai) ]
    genome                  // String
    spikein_genome          // String
    total_mapped_reads_key  // String

    main:

    ch_multiqc_files = channel.empty()

    //
    // MODULE: split BAMs by spike-in genome (keep endogenous)
    //
    BAM_SPLIT_BY_GENOME_ENDO (
        ch_bam,
        genome,
        spikein_genome,
        'endo'
    )

    //
    // MODULE: split BAMs by spike-in genome (keep exogenous)
    //
    BAM_SPLIT_BY_GENOME_EXO (
        ch_bam,
        genome,
        spikein_genome,
        'exo'
    )

    // add genome as meta field
    ch_bam_endo = BAM_SPLIT_BY_GENOME_ENDO.out.bam.map { it -> [ it[0] + [ genome: genome ], it[1] ] }
    ch_bam_exo = BAM_SPLIT_BY_GENOME_EXO.out.bam.map { it -> [ it[0] + [ genome: spikein_genome ], it[1] ] }
    ch_bam = ch_bam_endo.mix(ch_bam_exo)

    // Separate single-end and paired-end BAM files (SE do not have orphans)
    ch_bam
        .branch { meta, bam ->
            se: meta.single_end
            pe: !meta.single_end
        }
        .set { ch_bam }

    //
    // MODULE: Sort BAM files by query name
    //
    SAMTOOLS_NSORT (
        ch_bam.pe,
        ch_fasta_fai,
        ''
    )

    //
    // MODULE: Remove orphan reads left by splitting by genome (only for paired-end BAMs)
    //
    BAM_REMOVE_ORPHANS (
        SAMTOOLS_NSORT.out.bam,
        true
    )

    // Mixing SE and PE (removed orphans) files
    ch_bam = ch_bam.se.mix(BAM_REMOVE_ORPHANS.out.bam)

    //
    // MODULE: Sort BAM by coordinate and generate stats
    //
    // TODO: The genome fasta would need to be separated, right now both endo and exo stats
    // are analyzed with the same (main) genome
    // either way, this has no effect on BAM output (https://bioinformatics.stackexchange.com/a/4218)
    BAM_SORT_STATS_SAMTOOLS (
        ch_bam,
        ch_fasta_fai
    )
    ch_bam = BAM_SORT_STATS_SAMTOOLS.out.bam
    ch_index = BAM_SORT_STATS_SAMTOOLS.out.index
    ch_flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat
    ch_multiqc_files = ch_multiqc_files.mix(BAM_SORT_STATS_SAMTOOLS.out.stats.collect { it -> it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(BAM_SORT_STATS_SAMTOOLS.out.flagstat.collect { it -> it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(BAM_SORT_STATS_SAMTOOLS.out.idxstats.collect { it -> it[1] })

    //
    // MODULE: Extract total mapped reads from flagstats
    //
    BAM_FLAGSTAT_MAPPED (
        ch_flagstat
    )

    // Extract the total mapped reads from the text file
    BAM_FLAGSTAT_MAPPED.out.txt
        .map { meta, total ->
            [meta, total.splitCsv(header: false)[0][0]]
        }
        .set { ch_total_reads }

    // Add the total_mapped_reads both endo and exo bams' and indexes' metas
    ch_bam
        .join(ch_index, by: 0)
        .combine(ch_total_reads, by: 0)
        .map { meta, bam, index, total ->
            def meta_clone = meta.clone()
            meta_clone[total_mapped_reads_key] = total.toDouble()
            meta_clone.ref_total_mapped_reads_key = total_mapped_reads_key
            [meta_clone, bam, index]
        }
        .multiMap { meta, bam, index ->
            bam: [ meta, bam ]
            index: [ meta, index ]
        }
        .set { ch_bam_index }


    emit:

    bam           = ch_bam_index.bam
    endo_bam      = ch_bam_index.bam.filter { it -> it[0].genome == genome }                 // channel: [ val(meta), [ bam ] ]
    exo_bam       = ch_bam_index.bam.filter { it -> it[0].genome == spikein_genome }         // channel: [ val(meta), [ bam ] ]

    index         = ch_bam_index.index
    endo_index    = ch_bam_index.index.filter { it -> it[0].genome == genome }               // channel: [ val(meta), [ index ] ]
    exo_index     = ch_bam_index.index.filter { it -> it[0].genome == spikein_genome }       // channel: [ val(meta), [ index ] ]

    stats         = BAM_SORT_STATS_SAMTOOLS.out.stats
    endo_stats    = BAM_SORT_STATS_SAMTOOLS.out.stats.filter { it -> it[0].genome == genome }               // channel: [ val(meta), [ stats ] ]
    exo_stats     = BAM_SORT_STATS_SAMTOOLS.out.stats.filter { it -> it[0].genome == spikein_genome }       // channel: [ val(meta), [ stats ] ]

    flagstat      = ch_flagstat
    endo_flagstat = ch_flagstat.filter { it -> it[0].genome == genome }            // channel: [ val(meta), [ flagstat ] ]
    exo_flagstat  = ch_flagstat.filter { it -> it[0].genome == spikein_genome }    // channel: [ val(meta), [ flagstat ] ]

    idxstats      = BAM_SORT_STATS_SAMTOOLS.out.idxstats
    endo_idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats.filter { it -> it[0].genome == genome }            // channel: [ val(meta), [ idxstats ] ]
    exo_idxstats  = BAM_SORT_STATS_SAMTOOLS.out.idxstats.filter { it -> it[0].genome == spikein_genome }    // channel: [ val(meta), [ idxstats ] ]

    multiqc_files = ch_multiqc_files    // channel: [ multiqc_files ]

}
