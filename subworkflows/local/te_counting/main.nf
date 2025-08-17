//
// Counting reads in transposable elements (TEs)
//
include { SAMTOOLS_SORT             } from '../../../modules/nf-core/samtools/sort/main'
include { TECOUNT }                 from '../../../modules/local/tecount/main'
include { TELOCAL }                 from '../../../modules/local/telocal/telocal/main'

workflow TE_COUNTING {

    take:
    ch_bam
    ch_fasta
    skip_name_sort // boolean: skip name sorting of BAM files
    ch_tecount_gene_index
    ch_tecount_te_index
    ch_telocal_gene_index
    ch_telocal_te_index
    skip_telocal

    main:

    ch_versions = Channel.empty()


    if (!skip_name_sort) {
        //
        // MODULE: Sort BAM files (name sort order)
        //
        SAMTOOLS_SORT (
            ch_bam,
            ch_fasta
        )
        ch_bam = SAMTOOLS_SORT.out.bam
        ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())
    }

    //
    // MODULE: Count reads in transposable elements (TEs) at the subfamily level
    //
    TECOUNT (
        ch_bam,
        ch_tecount_gene_index,
        ch_tecount_te_index
    )
    ch_versions = ch_versions.mix(TECOUNT.out.versions.first())

    //
    // MODULE: Count reads in transposable elements (TEs) at the instance (location) level
    //
    ch_telocal_counts = Channel.empty()
    if (!skip_telocal) {
        TELOCAL (
            ch_bam,
            ch_telocal_gene_index,
            ch_telocal_te_index
        )
        ch_telocal_counts = TELOCAL.out.counts
        ch_versions = ch_versions.mix(TELOCAL.out.versions.first())
    }


    emit:

    tecount_counts     = TECOUNT.out.counts    // channel: [ te_counts.tsv ]
    telocal_counts     = ch_telocal_counts    // channel: [ te_local_counts.tsv ]
    
    versions    = ch_versions               // channel: [ versions.yml ]
}