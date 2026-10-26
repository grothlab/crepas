//
// Counting reads in transposable elements (TEs)
//
include { SAMTOOLS_SORT             } from '../../../modules/nf-core/samtools/sort/main'
include { TECOUNT }                 from '../../../modules/local/tecount/main'
include { TELOCAL }                 from '../../../modules/local/telocal/telocal/main'
include { TE_COUNTS_QC as TE_COUNTS_QC_TECOUNT } from '../../../modules/local/te_counts_qc/main'
include { TE_COUNTS_QC as TE_COUNTS_QC_TELOCAL } from '../../../modules/local/te_counts_qc/main'

//
// Group per-sample count tables into one QC run per experiment type, filtering stage,
// downsampling method and strand, so that a PCA only compares samples processed the same way
//
def groupCountsForQc(ch_counts) {
    return ch_counts
        .map { meta, counts ->
            // SCAR-seq/OK-seq tables are labelled by strand
            def strand = meta.te_counting_strandedness
            if (strand && meta.exp_type in ['SCAR-seq', 'OK-seq'] && meta.strandedness != 'forward') {
                strand = strand == 'forward' ? 'reverse' : 'forward'
            }
            def key = [
                id                       : meta.exp_type + (strand ? '.' + strand : ''),
                exp_type                 : meta.exp_type,
                te_counting_stage        : meta.te_counting_stage,
                downsampling_method      : meta.downsampling_method,
                te_counting_strandedness : strand
            ]
            [ key, counts ]
        }
        .groupTuple()
        .filter { _meta, counts -> counts.size() > 1 }
        .map { meta, counts -> [ meta, counts.sort { count -> count.name } ] }
}

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
    ch_te_counts_pca_header        // channel: [ header_file ]
    ch_te_counts_clustering_header // channel: [ header_file ]
    skip_te_counting_qc            // boolean: skip the PCA and sample distance QC of the count tables

    main:

    ch_multiqc_files = channel.empty()

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

    if (!skip_te_counting_qc) {
        //
        // MODULE: Normalise TEcount tables, PCA and sample distance heatmap
        //
        TE_COUNTS_QC_TECOUNT (
            groupCountsForQc(TECOUNT.out.counts),
            ch_te_counts_pca_header,
            ch_te_counts_clustering_header
        )
        ch_multiqc_files = ch_multiqc_files.mix(TE_COUNTS_QC_TECOUNT.out.pca_multiqc.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(TE_COUNTS_QC_TECOUNT.out.dists_multiqc.collect { it -> it[1] })

        if (!skip_telocal) {
            //
            // MODULE: Normalise TElocal tables, PCA and sample distance heatmap
            //
            TE_COUNTS_QC_TELOCAL (
                groupCountsForQc(TELOCAL.out.counts),
                ch_te_counts_pca_header,
                ch_te_counts_clustering_header
            )
            ch_multiqc_files = ch_multiqc_files.mix(TE_COUNTS_QC_TELOCAL.out.pca_multiqc.collect { it -> it[1] })
            ch_multiqc_files = ch_multiqc_files.mix(TE_COUNTS_QC_TELOCAL.out.dists_multiqc.collect { it -> it[1] })
        }
    }

    emit:

    tecount_counts     = TECOUNT.out.counts    // channel: [ te_counts.tsv ]
    telocal_counts     = ch_telocal_counts    // channel: [ te_local_counts.tsv ]
    multiqc_files      = ch_multiqc_files     // channel: [ multiqc_files ]

}
