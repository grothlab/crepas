//
// Alignment with strobealign
//

include { STROBEALIGN             } from "../../../modules/nf-core/strobealign/main"
include { BAM_SORT_STATS_SAMTOOLS } from '../../nf-core/bam_sort_stats_samtools/main'

workflow FASTQ_ALIGN_STROBEALIGN {
    take:
    ch_reads          // channel: [ val(meta), [ reads ] ]
    ch_strobealign_index
    ch_fasta_fai      // channel: [ val(meta), path(fasta), path(fai) ]
    val_sort_bam    // boolean (mandatory): true or false

    main:

    ch_fasta = ch_fasta_fai.map { meta, fasta, _fai -> [ meta, fasta ] }

    //
    // MODULE: Map reads with strobealign
    //
    STROBEALIGN  (
        ch_reads,
        ch_fasta,
        ch_strobealign_index,
        val_sort_bam
    )

    //
    // SUBWORKFLOW: Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_STATS_SAMTOOLS ( STROBEALIGN.out.bam, ch_fasta_fai )

    emit:
    bam_orig         = STROBEALIGN.out.bam                  // channel: [ val(meta), [ bam ] ]
    bam              = BAM_SORT_STATS_SAMTOOLS.out.bam      // channel: [ val(meta), [ bam ] ]
    index            = BAM_SORT_STATS_SAMTOOLS.out.index    // channel: [ val(meta), [ index ] ]
    stats            = BAM_SORT_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat         = BAM_SORT_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats         = BAM_SORT_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
}
