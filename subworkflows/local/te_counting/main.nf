//
// Counting reads in transposable elements (TEs)
//
include { SAMTOOLS_SORT             } from '../../../modules/nf-core/samtools/sort/main'
include { TECOUNT }                 from '../../../modules/local/tecount/main'
include { TELOCAL }                 from '../../../modules/local/telocal/telocal/main'
include { TE_COUNTING_QC as TECOUNT_QC } from '../../../modules/local/te_counting_qc/main'
include { TE_COUNTING_QC as TELOCAL_QC } from '../../../modules/local/te_counting_qc/main'


// TODO: Improve readability
def groupCountsForQc(ch_counts) {
    // Antibodies of the IPs that use each input control
    def ch_input_antibodies = ch_counts
        // IPs that have an input control
        .filter { meta, _counts -> meta.input_control }
        .map { meta, _counts -> [ meta.input_control, meta.antibody ] }
        .unique()

    // pre_flT3 and pre_flTbl tables of input controls that are the reference for multiple antibodies need to be copied
    return ch_counts
        // Input controls not yet assigned to an antibody (pre_flT3 and pre_flTbl stages, and final stage without downsampling)
        .filter { meta, _counts -> meta.is_input_control && !meta.input_control_of_antibody }
        .map { meta, counts -> [ meta.id, meta, counts ] }
        .combine(ch_input_antibodies, by: 0)
        .map { _id, meta, counts, antibody -> [ meta + [ input_control_of_antibody: antibody ], counts ] }
        // IPs, and input controls already duplicated for each antibody (final stage with downsampling)
        .mix(ch_counts.filter { meta, _counts -> !meta.is_input_control || meta.input_control_of_antibody })
        .map { meta, counts ->
            // SCAR-seq/OK-seq tables are labelled by strand
            def strand = meta.te_counting_strandedness
            if (strand && meta.exp_type in ['SCAR-seq', 'OK-seq'] && meta.strandedness != 'forward') {
                strand = strand == 'forward' ? 'reverse' : 'forward'
            }
            // IPs have meta.antibody, input controls have meta.input_control_of_antibody
            def antibody = meta.antibody ?: meta.input_control_of_antibody
            def key = [
                id                       : meta.exp_type + (antibody ? '_' + antibody : '') + (strand ? '.' + strand : ''),
                exp_type                 : meta.exp_type,
                antibody                 : antibody,
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
    ch_te_counting_qc_pca_header        // channel: [ header_file ]
    ch_te_counting_qc_clustering_header // channel: [ header_file ]
    skip_te_counting_qc                 // boolean: skip the PCA and sample distance QC of the count tables

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
        TECOUNT_QC (
            groupCountsForQc(TECOUNT.out.counts),
            ch_te_counting_qc_pca_header,
            ch_te_counting_qc_clustering_header
        )
        ch_multiqc_files = ch_multiqc_files.mix(TECOUNT_QC.out.pca_multiqc.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(TECOUNT_QC.out.dists_multiqc.collect { it -> it[1] })

        if (!skip_telocal) {
            //
            // MODULE: Normalise TElocal tables, PCA and sample distance heatmap
            //
            TELOCAL_QC (
                groupCountsForQc(TELOCAL.out.counts),
                ch_te_counting_qc_pca_header,
                ch_te_counting_qc_clustering_header
            )
            ch_multiqc_files = ch_multiqc_files.mix(TELOCAL_QC.out.pca_multiqc.collect { it -> it[1] })
            ch_multiqc_files = ch_multiqc_files.mix(TELOCAL_QC.out.dists_multiqc.collect { it -> it[1] })
        }
    }

    emit:

    tecount_counts     = TECOUNT.out.counts    // channel: [ te_counts.tsv ]
    telocal_counts     = ch_telocal_counts    // channel: [ te_local_counts.tsv ]
    multiqc_files      = ch_multiqc_files     // channel: [ multiqc_files ]

}
