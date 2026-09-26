//
// Alignment with Bowtie
//

include { BOWTIE_ALIGN           } from '../../../modules/nf-core/bowtie/align/main'
include { BAM_SORT_STATS_SAMTOOLS } from '../../nf-core/bam_sort_stats_samtools/main'

workflow FASTQ_ALIGN_BOWTIE {
    take:
    ch_reads          // channel: [ val(meta), [ reads ] ]
    ch_index          // channel: /path/to/bowtie2/index/
    save_unaligned    // val
    ch_fasta_fai      // channel: [ val(meta), path(fasta), path(fai) ]

    main:

    //
    // Map reads with Bowtie2
    //
    BOWTIE_ALIGN ( ch_reads, ch_index, save_unaligned )

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_STATS_SAMTOOLS ( BOWTIE_ALIGN.out.bam, ch_fasta_fai )

    emit:
    bam_orig         = BOWTIE_ALIGN.out.bam          // channel: [ val(meta), aligned ]
    log_out          = BOWTIE_ALIGN.out.log          // channel: [ val(meta), log     ]
    fastq            = BOWTIE_ALIGN.out.fastq        // channel: [ val(meta), fastq   ]

    bam              = BAM_SORT_STATS_SAMTOOLS.out.bam      // channel: [ val(meta), [ bam ] ]
    index            = BAM_SORT_STATS_SAMTOOLS.out.index    // channel: [ val(meta), [ index ] ]
    stats            = BAM_SORT_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat         = BAM_SORT_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats         = BAM_SORT_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
}
