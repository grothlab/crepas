//
// Counting reads in transposable elements (TEs)
//
include { SAMTOOLS_SORT             } from '../../../modules/nf-core/samtools/sort/main'
include { TECOUNT }                 from '../../../modules/local/tecount/main'
include { TELOCAL }                 from '../../../modules/local/telocal/telocal/main'

workflow TE_COUNTING {

    take:
    ch_bam
    ch_fasta_fai
    skip_name_sort // boolean: skip name sorting of BAM files
    ch_tecount_gene_index
    ch_tecount_te_index
    ch_telocal_gene_index
    ch_telocal_te_index
    skip_telocal
    skip_tecount_gz
    skip_telocal_gz

    main:

    if (!skip_name_sort) {
        //
        // MODULE: Sort BAM files (name sort order)
        //
        SAMTOOLS_SORT (
            ch_bam,
            ch_fasta_fai,
            ''
        )
        ch_bam = SAMTOOLS_SORT.out.bam
    }

    ch_te_counting_no_split = ch_bam.filter { it -> !(it[0].exp_type in ['SCAR-seq', 'OK-seq', 'eSPAN']) }
    ch_te_counting_split = ch_bam.filter { it -> it[0].exp_type in ['SCAR-seq', 'OK-seq', 'eSPAN'] }

    ch_te_counting_split
        .map { meta, bam -> [ meta + [ te_counting_strandedness: 'forward' ], bam ] }
        .set { ch_te_counting_split_fwd }

    ch_te_counting_split
        .map { meta, bam -> [ meta + [ te_counting_strandedness: 'reverse' ], bam ] }
        .set { ch_te_counting_split_rev }

    ch_te_counting_no_split
        .mix(ch_te_counting_split_fwd)
        .mix(ch_te_counting_split_rev)
        .set { ch_bam }

    //
    // MODULE: Count reads in transposable elements (TEs) at the subfamily level
    //
    TECOUNT (
        ch_bam,
        ch_tecount_gene_index,
        ch_tecount_te_index,
        skip_tecount_gz
    )

    //
    // MODULE: Count reads in transposable elements (TEs) at the instance (location) level
    //
    ch_telocal_counts = channel.empty()
    if (!skip_telocal) {
        TELOCAL (
            ch_bam,
            ch_telocal_gene_index,
            ch_telocal_te_index,
            skip_telocal_gz
        )
        ch_telocal_counts = TELOCAL.out.counts
    }

    emit:

    tecount_counts     = TECOUNT.out.counts    // channel: [ te_counts.tsv ]
    telocal_counts     = ch_telocal_counts    // channel: [ te_local_counts.tsv ]

}
