/*
 * Filter out blacklisted regions from BAM file using a whitelist BED file.
 */

include { BAM_FILTER_SAMBAMBA     } from '../../../subworkflows/local/bam_filter_sambamba/main'
include { SAMTOOLS_SORT           } from '../../../modules/nf-core/samtools/sort/main'
include { BAM_REMOVE_ORPHANS      } from '../../../modules/local/bam_remove_orphans/main'
include { BAM_SORT_STATS_SAMTOOLS } from '../../../subworkflows/nf-core/bam_sort_stats_samtools/main'
include { BAM_FLAGSTAT_MAPPED     } from '../../../modules/local/bam_flagstat_mapped/main'

workflow BAM_FILTER_BLACKLIST {
    take:
    ch_bam_index     // channel: [ val(meta), [ bam ], [ index ]]
    ch_whitelist_bed // channel: [ val(meta2), [ bed ] ]
    ch_fasta         // channel: [ fasta ]

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Filter BAM file with SAMBAMBA using blacklist
    //
    BAM_FILTER_SAMBAMBA(
        ch_bam_index,
        ch_whitelist_bed,
        ch_fasta
    )
    ch_filtered_bam = BAM_FILTER_SAMBAMBA.out.bam
    ch_filtered_index = BAM_FILTER_SAMBAMBA.out.bai
    ch_flagstat = BAM_FILTER_SAMBAMBA.out.flagstat
    ch_multiqc_files = BAM_FILTER_SAMBAMBA.out.stats.collect { it[1] }
    ch_multiqc_files = ch_multiqc_files.mix(BAM_FILTER_SAMBAMBA.out.flagstat.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(BAM_FILTER_SAMBAMBA.out.idxstats.collect { it[1] })
    ch_versions = ch_versions.mix(BAM_FILTER_SAMBAMBA.out.versions)

    // Separate single-end and paired-end BAM files (SE do not have orphans)
    ch_filtered_bam
        .branch { meta, bam ->
            se: meta.single_end
            pe: !meta.single_end
        }
        .set { ch_filtered_bam }

    ch_filtered_index
        .branch { meta, bam ->
            se: meta.single_end
            pe: !meta.single_end
        }
        .set { ch_filtered_index }

    ch_flagstat
        .branch { meta, flagstat ->
            se: meta.single_end
            pe: !meta.single_end
        }
        .set { ch_flagstat }

    //
    // MODULE: Sort BAM files by query name
    //
    SAMTOOLS_SORT(
        ch_filtered_bam.pe,
        ch_fasta
    )
    ch_filtered_bam_pe = SAMTOOLS_SORT.out.bam
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    //
    // MODULE: Remove orphan reads left by SAMBAMBA's blacklist filtering
    //
    BAM_REMOVE_ORPHANS(
        ch_filtered_bam_pe,
        true
    )
    ch_versions = ch_versions.mix(BAM_REMOVE_ORPHANS.out.versions.first())

    //
    // MODULE: Sort BAM by coordinate and generate stats
    //
    BAM_SORT_STATS_SAMTOOLS(
        BAM_REMOVE_ORPHANS.out.bam,
        ch_fasta
    )
    ch_multiqc_files = ch_multiqc_files.mix(BAM_SORT_STATS_SAMTOOLS.out.stats.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(BAM_SORT_STATS_SAMTOOLS.out.flagstat.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(BAM_SORT_STATS_SAMTOOLS.out.idxstats.collect { it[1] })
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)


    // Mixing SE and PE (removed orphans) files
    ch_filtered_bam = ch_filtered_bam.se.mix(BAM_SORT_STATS_SAMTOOLS.out.bam)
    ch_filtered_index = ch_filtered_index.se.mix(BAM_SORT_STATS_SAMTOOLS.out.bai)
    ch_flagstat = ch_flagstat.se.mix(BAM_SORT_STATS_SAMTOOLS.out.flagstat)

    //
    // MODULE: Extract total mapped reads from flagstats
    //
    BAM_FLAGSTAT_MAPPED(
        ch_flagstat
    )
    ch_versions = ch_versions.mix(BAM_FLAGSTAT_MAPPED.out.versions)

    // Extract the total mapped reads from the text file
    BAM_FLAGSTAT_MAPPED.out.txt
        .map { meta, total ->
            [meta, total.splitCsv(header: false)[0][0]]
        }
        .set { ch_flTbl_total }

    // Add the total_mapped_reads to the bams' and bais' metas
    ch_filtered_bam
        .combine(ch_filtered_index, by: 0)
        .map { meta, bam, bai ->
            [meta, bam, bai]
        }
        .combine(ch_flTbl_total, by: 0)
        .map { meta, bam, bai, total ->
            def meta_clone = meta.clone()
            meta_clone.flTbl_total_mapped_reads = total.toDouble()
            [meta_clone, bam, bai]
        }
        .set { ch_filtered_bam_bai }

    ch_filtered_bam_bai
        .map { meta, bam, bai ->
            [meta, bam]
        }
        .set { ch_filtered_bam }

    ch_filtered_bam_bai
        .map { meta, bam, bai ->
            [meta, bai]
        }
        .set { ch_filtered_index }


    // TODO: save for debugging
    ch_filtered_bam_bai
        .map { meta, bam, bai ->
            "${meta}\t${bam}\t${bai}"
        }
        .collectFile(name: 'ch_filtered_bam_bai.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_FILTER_BLACKLIST")

    emit:
    bam           = ch_filtered_bam // channel: [ val(meta), [ bam ] ]
    bai           = ch_filtered_index // channel: [ val(meta), [ bai ] ]
    stats         = BAM_FILTER_SAMBAMBA.out.stats.mix(BAM_SORT_STATS_SAMTOOLS.out.stats) // channel: [ val(meta), [ stats ] ]
    flagstat      = BAM_FILTER_SAMBAMBA.out.flagstat.mix(BAM_SORT_STATS_SAMTOOLS.out.flagstat) // channel: [ val(meta), [ flagstat ] ]
    idxstats      = BAM_FILTER_SAMBAMBA.out.idxstats.mix(BAM_SORT_STATS_SAMTOOLS.out.idxstats) // channel: [ val(meta), [ idxstats ] ]
    multiqc_files = ch_multiqc_files // channel: [ val(meta), [ multiqc_files ] ]
    versions      = ch_versions // channel: [ versions.yml ]
}
