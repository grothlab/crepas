//
// Alignment with minibwa
//

include { MINIBWA_MAP             } from '../../../modules/nf-core/minibwa/map/main'
include { BAM_SORT_STATS_SAMTOOLS } from '../../../subworkflows/nf-core/bam_sort_stats_samtools/main'

workflow FASTQ_ALIGN_MINIBWA {
    take:
    ch_reads        // channel (mandatory): [ val(meta), [ path(reads) ] ]
    ch_index        // channel (mandatory): [ val(meta2), path(index) ]
    val_sort_bam    // boolean (mandatory): true or false
    ch_fasta_fai    // channel (optional) : [ val(meta3), path(fasta), path(fai) ]

    main:

    //
    // MODULE: Map reads with MINIBWA_MAP (unsorted BAM when val_sort_bam is false)
    //
    ch_fasta = ch_fasta_fai.map { meta, fasta, _fai -> [ meta, fasta ] }

    MINIBWA_MAP ( ch_reads, ch_index, ch_fasta, val_sort_bam )

    //
    // SUBWORKFLOW: Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_STATS_SAMTOOLS ( MINIBWA_MAP.out.aligned, ch_fasta_fai )

    emit:
    bam_orig = MINIBWA_MAP.out.aligned              // channel: [ val(meta), path(bam) ]
    bam      = BAM_SORT_STATS_SAMTOOLS.out.bam      // channel: [ val(meta), path(bam) ]
    index    = BAM_SORT_STATS_SAMTOOLS.out.index    // channel: [ val(meta), path(index) ]
    stats    = BAM_SORT_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), path(stats) ]
    flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), path(flagstat) ]
    idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), path(idxstats) ]

}
