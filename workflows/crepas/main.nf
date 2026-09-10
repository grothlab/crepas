/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { IGV                                                         } from '../../modules/local/igv/main'
include { MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS                         } from '../../modules/local/multiqc_custom_phantompeakqualtools/main'

//
// SUBWORKFLOWS: Consisting of a mix of local and nf-core/modules
//
include { samplesheetToList                } from 'plugin/nf-schema'
include { paramsSummaryMap                                                  } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                                              } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                                            } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                                            } from '../../subworkflows/local/utils_grothlab_crepas_pipeline'
include { rtElFractions                                                    } from '../../subworkflows/local/utils_grothlab_crepas_pipeline'
include { rtHrFractions                                                    } from '../../subworkflows/local/utils_grothlab_crepas_pipeline'
include { INPUT_CHECK                                                       } from '../../subworkflows/local/utils_grothlab_crepas_pipeline'

include {
    BAM_FILTER_SAMBAMBA_RMO_STATS as BAM_FILTER_SAMBAMBA_FLT1
    BAM_FILTER_SAMBAMBA_RMO_STATS as BAM_FILTER_SAMBAMBA_FLT3
    BAM_FILTER_SAMBAMBA_RMO_STATS as BAM_FILTER_SAMBAMBA_BLACKLIST
    } from '../../subworkflows/local/bam_filter_sambamba_rmo_stats/main'

include { BAM_SPIKEIN_SPLIT                                                 } from '../../subworkflows/local/bam_spikein_split/main'
include { FASTQ_FASTQC_UMITOOLS_UMITRANSFER_TRIMGALORE                      } from '../../subworkflows/local/fastq_fastqc_umitools_umitransfer_trimgalore/main'
include { BAM_ENCODE_PIPELINE                                               } from '../../subworkflows/local/bam_encode_pipeline/main'
include { BAM_CREATE_PARTITIONS                                             } from '../../subworkflows/local/bam_create_partitions/main'
include { BAM_EL_REPLISEQ                                                   } from '../../subworkflows/local/bam_el_repliseq/main'
include { BAM_HR_REPLISEQ                                                   } from '../../subworkflows/local/bam_hr_repliseq/main'
include { BAM_ALLOCATE_MULTIMAPPERS                                         } from '../../subworkflows/local/bam_allocate_multimappers/main'
include { BAM_SHIFT_READS                                                   } from '../../subworkflows/local/bam_shift_reads/main'
include { SAMTOOLS_STATS_SUMMARY                                            } from '../../subworkflows/local/samtools_stats_summary/main'
include { BAM_NORMALIZE_BIGWIG_DEEPTOOLS                                    } from '../../subworkflows/local/bam_normalize_bigwig_deeptools/main'
include { BAM_DOWNSAMPLE                                                    } from '../../subworkflows/local/bam_downsample/main'
include { TE_COUNTING                                                       } from '../../subworkflows/local/te_counting/main'
include { FASTQ_ALIGN                                                       } from '../../subworkflows/local/fastq_align/main'
include { SPIKEIN_BARCODES                                                  } from '../../subworkflows/local/spikein_barcodes/main'
include { CALL_PEAKS                                                        } from '../../subworkflows/local/call_peaks/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULES: Installed directly from nf-core/modules
//
include { SAMTOOLS_INDEX                                                    } from '../../modules/nf-core/samtools/index/main'
include { PICARD_MERGESAMFILES                                              } from '../../modules/nf-core/picard/mergesamfiles/main'
include { PICARD_COLLECTMULTIPLEMETRICS                                     } from '../../modules/nf-core/picard/collectmultiplemetrics/main'
include { PRESEQ_LCEXTRAP                                                   } from '../../modules/nf-core/preseq/lcextrap/main'
include { PHANTOMPEAKQUALTOOLS                                              } from '../../modules/local/phantompeakqualtools/main'
include { DEEPTOOLS_COMPUTEMATRIX as DEEPTOOLS_COMPUTEMATRIX_GENES          } from '../../modules/nf-core/deeptools/computematrix/main'
include { DEEPTOOLS_PLOTPROFILE as DEEPTOOLS_PLOTPROFILE_GENES              } from '../../modules/nf-core/deeptools/plotprofile/main'
include { DEEPTOOLS_PLOTHEATMAP as DEEPTOOLS_PLOTHEATMAP_GENES              } from '../../modules/nf-core/deeptools/plotheatmap/main'
include { DEEPTOOLS_PLOTFINGERPRINT                                         } from '../../modules/nf-core/deeptools/plotfingerprint/main'
include { MULTIQC                                                           } from '../../modules/nf-core/multiqc/main'

//
// SUBWORKFLOWS: Installed directly from nf-core/modules
//
include { BAM_MARKDUPLICATES_PICARD                                         } from '../../subworkflows/nf-core/bam_markduplicates_picard'
include { BAM_DEDUP_UMI                                                     } from '../../subworkflows/local/bam_dedup_umi'
include { BAM_STATS_SAMTOOLS                                                } from '../../subworkflows/nf-core/bam_stats_samtools'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CREPAS {
    take:
    ch_versions               // channel: [ path(versions.yml) ]
    ch_fasta                  // channel: path(genome.fa)
    ch_fai                    // channel: path(genome.fai)
    ch_chrom_sizes
    ch_endo_chromsizes       // path(chrom.sizes.endo)
    ch_exo_chromsizes
    ch_gtf                    // channel: path(genome.gtf)
    ch_gene_bed               // channel: path(gene.beds)
    ch_effective_gsize
    ch_effective_gfraction
    ch_whitelist           // channel: path(filtered.bed)
    ch_blacklist              // channel: path(blacklist.bed)
    ch_spikein_barcode_table  // channel: [ val(meta), path(spikein_barcode_table.tsv) ]
    ch_sparsebed              // channel: path(sparse.bed)
    ch_active_regions         // channel: path(active_regions.bed)
    ch_rocco_params           // channel: path(params.csv)
    ch_okseq_rfd_file       // channel: [ val(meta), [ bed ] ]
    ch_initiation_zones       // channel: path(initiation_zones)
    ch_bwa_index              // channel: path(bwa/index/)
    ch_bwamem2_index        // channel: path(bwamem2/index/)
    ch_minibwa_index        // channel: path(minibwa/index/)
    ch_bowtie_index            // channel: path(bowtie/index/)
    ch_bowtie2_index          // channel: path(bowtie2/index)
    ch_chromap_index          // channel: path(chromap.index)
    ch_star_index             // channel: path(star/index/)
    ch_hisat2_index           // channel: path(hisat2/index)
    ch_minimap2_index         // channel: path(minimap2/index/)
    ch_strobealign_index      // channel: path(strobealign/index)
    ch_splicesites            // channel: path(splicesites)
    ch_tecount_gene_index // channel: val(meta), path(tecount_gene_index.Ind)
    ch_telocal_gene_index // channel: val(meta), path(telocal_gene_index.Ind)
    ch_tecount_te_index       // channel: val(meta), path(tecount_te_index.Ind)
    ch_telocal_te_index       // channel: val(meta), path(telocal_te_index.locInd)

    main:
    ch_multiqc_files = channel.empty()
    ch_samtools_stats_summary = channel.empty()

    // TODO: organize these:

    // Header files for MultiQC
    ch_spp_nsc_header = file("${projectDir}/assets/multiqc/spp_nsc_header.txt", checkIfExists: true)
    ch_spp_rsc_header = file("${projectDir}/assets/multiqc/spp_rsc_header.txt", checkIfExists: true)
    ch_spp_correlation_header = file("${projectDir}/assets/multiqc/spp_correlation_header.txt", checkIfExists: true)
    ch_macs3_peak_count_header = file("${projectDir}/assets/multiqc/peak_count_header.txt", checkIfExists: true)
    ch_gr_peak_count_header = file("${projectDir}/assets/multiqc/gr_peak_count_header.txt", checkIfExists: true)
    ch_mace_peak_count_header = file("${projectDir}/assets/multiqc/mace_peak_count_header.txt", checkIfExists: true)
    ch_epic2_peak_count_header = file("${projectDir}/assets/multiqc/epic2_peak_count_header.txt", checkIfExists: true)
    ch_seacr_peak_count_header = file("${projectDir}/assets/multiqc/seacr_peak_count_header.txt", checkIfExists: true)
    ch_macs3_frip_score_header = file("${projectDir}/assets/multiqc/frip_score_header.txt", checkIfExists: true)
    ch_gr_frip_score_header = file("${projectDir}/assets/multiqc/gr_frip_score_header.txt", checkIfExists: true)
    ch_mace_frip_score_header = file("${projectDir}/assets/multiqc/mace_frip_score_header.txt", checkIfExists: true)
    ch_epic2_frip_score_header = file("${projectDir}/assets/multiqc/epic2_frip_score_header.txt", checkIfExists: true)
    ch_seacr_frip_score_header = file("${projectDir}/assets/multiqc/seacr_frip_score_header.txt", checkIfExists: true)
    ch_macs3_peak_annotation_header = file("${projectDir}/assets/multiqc/peak_annotation_header.txt", checkIfExists: true)
    ch_gr_peak_annotation_header = file("${projectDir}/assets/multiqc/gr_peak_annotation_header.txt", checkIfExists: true)
    ch_mace_peak_annotation_header = file("${projectDir}/assets/multiqc/mace_peak_annotation_header.txt", checkIfExists: true)
    ch_epic2_peak_annotation_header = file("${projectDir}/assets/multiqc/epic2_peak_annotation_header.txt", checkIfExists: true)
    ch_seacr_peak_annotation_header = file("${projectDir}/assets/multiqc/seacr_peak_annotation_header.txt", checkIfExists: true)
    ch_deseq2_pca_header = channel.value(file("${projectDir}/assets/multiqc/deseq2_pca_header.txt", checkIfExists: true))
    ch_repliseq_rt_header = channel.value(file("${projectDir}/assets/multiqc/repliseq_rt_header.txt", checkIfExists: true))
    ch_repliseq_gene_class_header = channel.value(file("${projectDir}/assets/multiqc/repliseq_gene_class_header.txt", checkIfExists: true))
    ch_hr_repliseq_features_header = channel.value(file("${projectDir}/assets/multiqc/hr_repliseq_features_header.txt", checkIfExists: true))
    ch_deseq2_clustering_header = channel.value(file("${projectDir}/assets/multiqc/deseq2_clustering_header.txt", checkIfExists: true))

    //
    // Create channel from input file provided through params.input
    //
    channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .map {
            meta, fastq_1, fastq_2, fastq_umi ->
                def meta_clone = meta.clone()
                if (!fastq_2) {
                    meta_clone.single_end = true
                    if (!fastq_umi) {
                        meta_clone.sep_umi_fq = false
                        return [ meta_clone, [ fastq_1 ] ]
                    } else {
                        meta_clone.sep_umi_fq = true
                        return [ meta_clone, [ fastq_1, fastq_umi ] ]
                    }
                } else {
                    meta_clone.single_end = false
                    if (!fastq_umi) {
                        meta_clone.sep_umi_fq = false
                        return [ meta_clone, [ fastq_1, fastq_2 ] ]
                    } else {
                        meta_clone.sep_umi_fq = true
                        return [ meta_clone, [ fastq_1, fastq_2, fastq_umi ] ]
                    }
                }
        }
        .set { ch_fastq }

    //
    // SUBWORKFLOW: Extra validation of input samplesheet
    //
    INPUT_CHECK (
        ch_fastq,
        params.seq_center
    )

    //
    // SUBWORKFLOW: Read QC and trim adapters
    //
    FASTQ_FASTQC_UMITOOLS_UMITRANSFER_TRIMGALORE (
        INPUT_CHECK.out.fastq,
        params.skip_fastqc || params.skip_qc,
        params.skip_spikein_barcode_extract,
        ch_spikein_barcode_table,
        params.skip_trimming,
        params.umi_discard_read,
        params.min_trimmed_reads,
        params.hardtrim3_length,
        params.hardtrim5_length
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_UMITRANSFER_TRIMGALORE.out.fastqc_zip.collect { it -> it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_UMITRANSFER_TRIMGALORE.out.trim_zip.collect { it -> it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_UMITRANSFER_TRIMGALORE.out.trim_log.collect { it -> it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_UMITRANSFER_TRIMGALORE.out.no_sep_umi_fq_log.collect { it -> it[1] })


    ch_fasta_fai = ch_fasta.combine(ch_fai.map { meta, fai -> fai }).first()

    //
    // SUBWORKFLOW: Alignment
    //
    ch_genome_bam = channel.empty()
    FASTQ_ALIGN (
        FASTQ_FASTQC_UMITOOLS_UMITRANSFER_TRIMGALORE.out.reads,
        ch_fasta_fai,
        params.aligner,
        ch_bwa_index,
        ch_bwamem2_index,
        ch_minibwa_index,
        ch_bowtie_index,
        ch_bowtie2_index,
        ch_chromap_index,
        ch_star_index,
        ch_hisat2_index,
        ch_minimap2_index,
        ch_strobealign_index,
        ch_gtf,
        ch_splicesites,
        params.save_unaligned,
        false

    )
    ch_genome_bam = FASTQ_ALIGN.out.bam
    ch_samtools_stats_summary = ch_samtools_stats_summary.mix(FASTQ_ALIGN.out.samtools_stats_summary)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN.out.multiqc_files)

    //
    // MODULE: Merge resequenced BAM files
    //
    ch_genome_bam
        .map { meta, bam ->
            def meta_clone = meta.clone()
            meta_clone.remove('read_group')
            meta_clone.id = meta_clone.id.split('_')[0..-3].join('_')
            if (meta_clone.input_control) {
                meta_clone.input_control = meta_clone.input_control.split('_')[0..-3].join('_')
            }
            def key = groupKey(meta_clone.id, meta_clone.trep_count) // trep_count defined in INPUT_CHECK subworkflow
            [key, meta_clone, bam]
        }
        .groupTuple(by: 0)
        .map { key, metas, bams ->
            def sorted_metas = metas.sort { meta -> meta.trep }
            def meta_clone = sorted_metas[0].clone()
            meta_clone.remove('trep')
            def sorted_bams = bams.sort { bam -> bam.name }
            [meta_clone, sorted_bams]
        }
        .set { ch_sort_bam }

    // TODO: print for debugging
    ch_sort_bam
        .map {
            meta, bam ->
                "${meta}\t${bam}"
        }
        .collectFile( name: 'ch_sort_bam.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/" )


    PICARD_MERGESAMFILES (
        ch_sort_bam
    )
    ch_merged_bam = PICARD_MERGESAMFILES.out.bam

    SAMTOOLS_INDEX (
        ch_merged_bam
    )
    ch_merged_index = SAMTOOLS_INDEX.out.index
    ch_merged_bam_index = ch_merged_bam.join(ch_merged_index, by: 0)

    BAM_STATS_SAMTOOLS (
        ch_merged_bam_index,
        ch_fasta_fai
    )
    ch_samtools_stats_summary = ch_samtools_stats_summary.mix(BAM_STATS_SAMTOOLS.out.stats)
    ch_multiqc_files = ch_multiqc_files.mix(BAM_STATS_SAMTOOLS.out.stats.collect { it -> it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(BAM_STATS_SAMTOOLS.out.flagstat.collect { it -> it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(BAM_STATS_SAMTOOLS.out.idxstats.collect { it -> it[1] })


    //
    // SUBWORKFLOW: Allocation of multimappers
    //
    if (params.multimap_allocation_method && params.multimap_allocation_method != 'chromap') {
        BAM_ALLOCATE_MULTIMAPPERS (
            ch_merged_bam,
            ch_fasta_fai,
            params.multimap_allocation_method
        )
        ch_merged_bam = BAM_ALLOCATE_MULTIMAPPERS.out.bam
        ch_merged_index = BAM_ALLOCATE_MULTIMAPPERS.out.index
        ch_merged_bam_index = ch_merged_bam.join(ch_merged_index, by: 0)
        ch_samtools_stats_summary = ch_samtools_stats_summary.mix(BAM_ALLOCATE_MULTIMAPPERS.out.stats)
        ch_multiqc_files = ch_multiqc_files.mix(BAM_ALLOCATE_MULTIMAPPERS.out.stats.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(BAM_ALLOCATE_MULTIMAPPERS.out.flagstat.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(BAM_ALLOCATE_MULTIMAPPERS.out.idxstats.collect { it -> it[1] })

    }

    //
    // MODULE: Preseq coverage analysis
    //
    // TODO: this is done on the bams with spike-in included
    if (!params.skip_preseq) {
        PRESEQ_LCEXTRAP(
            ch_merged_bam
        )
        ch_multiqc_files = ch_multiqc_files.mix(PRESEQ_LCEXTRAP.out.lc_extrap.collect { it -> it[1] })
    }

    ch_merged_bam_index_with_umi = ch_merged_bam_index.filter { meta, bam, index -> meta.with_umi }
    ch_merged_index_without_umi = ch_merged_index.filter { meta, index -> !meta.with_umi }
    ch_merged_bam_without_umi = ch_merged_bam.filter { meta, bam -> !meta.with_umi }

    //
    // SUBWORKFLOW: Deduplicate BAM files
    //
    ch_transcriptome_bam = channel.empty()
    ch_transcriptome_fasta = channel.empty()
    ch_umidedup_bam = channel.empty()
    ch_umidedup_index = channel.empty()
    BAM_DEDUP_UMI (
        ch_merged_bam_index_with_umi,
        ch_chrom_sizes,
        [],
        params.umi_dedup_tool,
        params.get_dedup_stats,
        ch_transcriptome_bam,
        ch_transcriptome_fasta,
        params.skip_split_by_chrom
    )
    ch_umidedup_bam = BAM_DEDUP_UMI.out.bam
    ch_umidedup_index = BAM_DEDUP_UMI.out.index
    ch_samtools_stats_summary = ch_samtools_stats_summary.mix(BAM_DEDUP_UMI.out.stats)
    ch_multiqc_files = ch_multiqc_files.mix(BAM_DEDUP_UMI.out.multiqc_files)


    if (!params.skip_markduplicates ) {
        //
        // SUBWORKFLOW: Mark duplicates & filter BAM files
        //
        ch_mkdup_bam = channel.empty()
        ch_mkdup_index = channel.empty()
        BAM_MARKDUPLICATES_PICARD (
            ch_merged_bam_without_umi,
            ch_fasta_fai
        )
        ch_mkdup_bam = BAM_MARKDUPLICATES_PICARD.out.bam
        ch_mkdup_index = BAM_MARKDUPLICATES_PICARD.out.index
        ch_samtools_stats_summary = ch_samtools_stats_summary.mix(BAM_MARKDUPLICATES_PICARD.out.stats)
        ch_multiqc_files = ch_multiqc_files.mix(BAM_MARKDUPLICATES_PICARD.out.stats.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(BAM_MARKDUPLICATES_PICARD.out.flagstat.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(BAM_MARKDUPLICATES_PICARD.out.idxstats.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(BAM_MARKDUPLICATES_PICARD.out.metrics.collect { it -> it[1] })

        ch_dedup_bam = ch_umidedup_bam.mix(ch_mkdup_bam)
        ch_dedup_index = ch_umidedup_index.mix(ch_mkdup_index)

    } else {

        ch_dedup_bam = ch_umidedup_bam.mix(ch_merged_bam_without_umi)
        ch_dedup_index = ch_umidedup_index.mix(ch_merged_index_without_umi)
    }

    //
    // SUBWORKFLOW: Filter BAM file with SAMBAMBA
    //
    BAM_FILTER_SAMBAMBA_FLT1 (
        ch_dedup_bam.join(ch_dedup_index, by: 0),
        channel.value([[:], []]),
        ch_fasta_fai,
        true, // skip orphan removal
        'flT1_total_mapped_reads'
    )
    ch_filtered_bam = BAM_FILTER_SAMBAMBA_FLT1.out.bam
    ch_filtered_index = BAM_FILTER_SAMBAMBA_FLT1.out.index
    ch_flT1_total = BAM_FILTER_SAMBAMBA_FLT1.out.total_reads
    ch_samtools_stats_summary = ch_samtools_stats_summary.mix(BAM_FILTER_SAMBAMBA_FLT1.out.stats)
    ch_multiqc_files = ch_multiqc_files.mix(BAM_FILTER_SAMBAMBA_FLT1.out.stats.collect { it -> it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(BAM_FILTER_SAMBAMBA_FLT1.out.flagstat.collect { it -> it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(BAM_FILTER_SAMBAMBA_FLT1.out.idxstats.collect { it -> it[1] })

    if (!params.skip_spikein_barcode_extract) {

        //
        // MODULE: Merge spikein barcode counts of resequenced samples
        //
        SPIKEIN_BARCODES (
            FASTQ_FASTQC_UMITOOLS_UMITRANSFER_TRIMGALORE.out.barcode_counts,
            ch_flT1_total
        )

    }

    //
    // SUBWORKFLOW: Spike-in splitting
    //
    // TODO: if fasta and gtf are specified but not genome, val keep_genome_string in
    // BAM_SPLIT_BY_GENOME will fail
    ch_filtered_exo_bam = channel.empty()
    ch_filtered_exo_index = channel.empty()
    if (params.spikein_genome) {
        BAM_SPIKEIN_SPLIT (
            ch_filtered_bam,
            ch_fasta_fai,
            params.genome,
            params.spikein_genome,
            'flT2_total_mapped_reads'
        )
        ch_filtered_bam = BAM_SPIKEIN_SPLIT.out.endo_bam
        ch_filtered_exo_bam = BAM_SPIKEIN_SPLIT.out.exo_bam
        ch_filtered_index = BAM_SPIKEIN_SPLIT.out.endo_index
        ch_filtered_exo_index = BAM_SPIKEIN_SPLIT.out.exo_index
        ch_samtools_stats_summary = ch_samtools_stats_summary.mix(BAM_SPIKEIN_SPLIT.out.stats)
        ch_multiqc_files = ch_multiqc_files.mix(BAM_SPIKEIN_SPLIT.out.multiqc_files)

    } else {
        // If no spike-in genome add genome to metas
        ch_filtered_bam
            .join(ch_filtered_index, by: 0)
            .map { meta, bam, index ->
                def meta_clone = meta.clone()
                meta_clone.genome = params.genome
                [meta_clone, bam, index]
            }
            .set { ch_filtered_bam_index }

        ch_filtered_bam = ch_filtered_bam_index.map { meta, bam, index -> [meta, bam] }
        ch_filtered_index = ch_filtered_bam_index.map { meta, bam, index -> [meta, index] }
    }


    //
    // MODULE: Picard post alignment QC
    //
    if (!params.skip_picard_metrics) {
        PICARD_COLLECTMULTIPLEMETRICS (
            ch_filtered_bam.join(ch_filtered_index, by: 0),
            ch_fasta,
            ch_fai
        )
        ch_multiqc_files = ch_multiqc_files.mix(PICARD_COLLECTMULTIPLEMETRICS.out.metrics.collect { it -> it[1] })
    }

    // Mix the exogenous and endogenous BAM and index files
    ch_filtered_bam
        .mix(ch_filtered_exo_bam)
        .set { ch_filtered_bam }

    ch_filtered_index
        .mix(ch_filtered_exo_index)
        .set { ch_filtered_index }

    // Split the BAM and indexes into 'ATAC-seq' and other (for shifting)
    ch_filtered_bam
        .branch { meta, bam ->
            atacseq: meta.exp_type == 'ATAC-seq'
            other: meta.exp_type != 'ATAC-seq'
        }
        .set { ch_filtered_bam }

    ch_filtered_index
        .branch { meta, index ->
            atacseq: meta.exp_type == 'ATAC-seq'
            other: meta.exp_type != 'ATAC-seq'
        }
        .set { ch_filtered_index }

    //
    // SUBWORKFLOW: Shift ATAC-seq reads
    //
    BAM_SHIFT_READS (
        ch_filtered_bam.atacseq.join(ch_filtered_index.atacseq, by: 0),
        ch_fasta_fai
    )
    ch_filtered_bam = ch_filtered_bam.other.mix(BAM_SHIFT_READS.out.bam)
    ch_filtered_index = ch_filtered_index.other.mix(BAM_SHIFT_READS.out.index)
    ch_filtered_bam_index = ch_filtered_bam.join(ch_filtered_index, by: 0)

    ch_pre_flT3_endo_bam = channel.empty()
    if (!params.skip_flT3) {

        // These are to later run TE counting on both pre- and post-blacklist-filtering BAM files
        ch_pre_flT3_endo_bam = ch_filtered_bam.filter { it -> it[0].genome == params.genome }

        //
        // MODULE: Final filtering of BAM file with SAMBAMBA (quality filtering)
        //
        BAM_FILTER_SAMBAMBA_FLT3 (
            ch_filtered_bam.join(ch_filtered_index, by: 0),
            channel.value([[:], []]),
            ch_fasta_fai,
            params.skip_te_counting, // do orphan removal if TE counting downstream
            'flT3_total_mapped_reads'

        )
        ch_filtered_bam = BAM_FILTER_SAMBAMBA_FLT3.out.bam
        ch_filtered_index = BAM_FILTER_SAMBAMBA_FLT3.out.index
        ch_filtered_bam_index = ch_filtered_bam.join(ch_filtered_index, by: 0)
        ch_samtools_stats_summary = ch_samtools_stats_summary.mix(BAM_FILTER_SAMBAMBA_FLT3.out.stats)
        ch_multiqc_files = ch_multiqc_files.mix(BAM_FILTER_SAMBAMBA_FLT3.out.stats.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(BAM_FILTER_SAMBAMBA_FLT3.out.flagstat.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(BAM_FILTER_SAMBAMBA_FLT3.out.idxstats.collect { it -> it[1] })

    }

    ch_pre_flTbl_endo_bam = channel.empty()
    if (!params.skip_flTbl) {

        // These are to later run TE counting on both pre- and post-blacklist-filtering BAM files
        ch_pre_flTbl_endo_bam = ch_filtered_bam.filter { it -> it[0].genome == params.genome }

        // Separating endogenous and exogenous samples
        // TODO: could add a param "exo_blacklist" to use a different blacklist
        // for exogenous samples instead of skipping them
        ch_filtered_bam_index
            .branch { meta, bam, index ->
                endo: meta.genome == params.genome
                exo: meta.genome == params.spikein_genome
            }
            .set { ch_flt_bam_index_by_genome }

        ch_flt_bam_index_by_genome
            .exo
            // .map { meta, bam, index ->
            //     def meta_clone = meta.clone()
            //     meta_clone.flTbl_total_mapped_reads = meta[meta.ref_total_mapped_reads_key]
            //     meta_clone.ref_total_mapped_reads_key = 'flTbl_total_mapped_reads'
            //     [meta_clone, bam, index]
            // }
            .set { ch_flt_bam_index_by_genome_exo }

        // TODO: print for debugging
        ch_flt_bam_index_by_genome.endo
            .map { it -> "${it}"}
            .collectFile(name: 'ch_flt_bam_index_by_genome_endo_flTbl.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/")

        //
        // SUBWORKFLOW: Filter BAM file with SAMBAMBA using blacklist (whitelist)
        //
        BAM_FILTER_SAMBAMBA_BLACKLIST (
            ch_flt_bam_index_by_genome.endo,
            ch_whitelist,
            ch_fasta_fai,
            false, // do not skip orphan removal
            'flTbl_total_mapped_reads'
        )
        ch_filtered_bam = BAM_FILTER_SAMBAMBA_BLACKLIST.out.bam.mix(ch_flt_bam_index_by_genome_exo.map { meta, bam, index -> [meta, bam] })
        ch_filtered_index = BAM_FILTER_SAMBAMBA_BLACKLIST.out.index.mix(ch_flt_bam_index_by_genome_exo.map { meta, bam, index -> [meta, index] })
        ch_filtered_bam_index = ch_filtered_bam.join(ch_filtered_index, by: 0)
        ch_samtools_stats_summary = ch_samtools_stats_summary.mix(BAM_FILTER_SAMBAMBA_BLACKLIST.out.stats)
        ch_multiqc_files = ch_multiqc_files.mix(BAM_FILTER_SAMBAMBA_BLACKLIST.out.multiqc_files)

    }

    // TODO: print for debugging
    ch_filtered_bam_index
        .map { meta, bam, index ->
            "${meta}\t${bam}\t${index}"
        }
        .collectFile(name: 'ch_filtered_bam_index_flTbl.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/")


    //
    // MODULE: Phantompeaktools strand cross-correlation and QC metrics
    //
    if (!params.skip_spp_qc) {

        // Remove exogenous samples from SPP analysis
        ch_filtered_bam
            .filter { meta, bam ->
                meta.genome == params.genome
            }
            .map { meta, bam -> [meta, bam, []] }
            .set { ch_filtered_bam_for_spp }

        PHANTOMPEAKQUALTOOLS(
            ch_filtered_bam_for_spp
        )
        ch_multiqc_files = ch_multiqc_files.mix(PHANTOMPEAKQUALTOOLS.out.ccscores.collect { it -> it[1] })

        //
        // MODULE: MultiQC custom content for Phantompeaktools
        //
        MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS (
            PHANTOMPEAKQUALTOOLS.out.ccscores.join(PHANTOMPEAKQUALTOOLS.out.rdata, by: 0),
            ch_spp_nsc_header,
            ch_spp_rsc_header,
            ch_spp_correlation_header
        )
        ch_multiqc_files = ch_multiqc_files.mix(MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS.out.nsc.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS.out.rsc.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS.out.correlation.collect { it -> it[1] })
    }

    //
    // Duplicate input controls for each antibody
    // This is done for the cases where one input control is used for multiple IPs
    //
    ch_filtered_bam_index
        .branch { meta, bam, index ->
            ips_with_ipcontrol: meta.input_control
                return [meta.input_control, meta.antibody, meta, bam, index]
            ips_wo_ipcontrol: !meta.input_control && !meta.is_input_control
                    return [meta, bam, index]
            ipcontrols: !meta.input_control && meta.is_input_control
                return [meta.id, meta, bam, index]
        }
        .set { ch_bam_index_by_type }

    ch_bam_index_by_type.ipcontrols
        .combine(ch_bam_index_by_type.ips_with_ipcontrol, by: 0) // combine by control id only
        .map { ipcontrol_id, ipcontrol_meta, ipcontrol_bam, ipcontrol_index, ip_antibody, ip_meta, ip_bam, ip_index ->
            def meta_clone = ipcontrol_meta.clone()
            meta_clone.input_control_of_antibody = ip_antibody
            [ meta_clone, ipcontrol_bam, ipcontrol_index ]
        }
        .unique()
        .set { ch_bam_index_ipcontrols }

    ch_bam_index_by_type
        .ips_with_ipcontrol
        .map { ipcontrol_id, antibody, meta, bam, index ->
            [ meta, bam, index ]
        }
        .mix(ch_bam_index_by_type.ips_wo_ipcontrol)
        .mix(ch_bam_index_ipcontrols)
        .set { ch_filtered_bam_index }

    //
    // Define the reference total mapped reads key to be used for downsampling, normalization, etc.
    // Note: this is the place to define more complex rules if needed
    // For example: if one wants to prefer flT2 for downsampling, but not for normalization,
    // then change definition of meta.ref_total_mapped_reads_for_dSp_key here:
    //
    ch_filtered_bam_index
        .map { meta, bam, index ->
            def meta_clone = meta.clone()
            // samples have meta.antibody, while input controls have meta.input_control_of_antibody
            def antibody = meta.antibody ?: meta.input_control_of_antibody
            def total_key = meta.ref_total_mapped_reads_key
            // Save the last filtering total mapped reads
            meta_clone.last_total_mapped_reads_key = total_key
            if (antibody in params.use_flT2_as_total_ref.split(',').collect { it -> it.trim() }) {
                if (meta.flT2_total_mapped_reads) {
                    total_key = 'flT2_total_mapped_reads'
                } else {
                    total_key = 'flT1_total_mapped_reads'
                }
            }
            meta_clone.ref_total_mapped_reads_key = total_key
            meta_clone.ref_total_mapped_reads_for_dSp_key = total_key
            def norm_key = params.bam_downsampling_method ? 'dSp_total_mapped_reads' : total_key
            meta_clone.ref_total_mapped_reads_for_rpm_key = norm_key
            meta_clone.ref_total_mapped_reads_for_srpm_key = norm_key
            meta_clone.ref_total_mapped_reads_for_cisrpm_key = norm_key
            [meta_clone, bam, index]
        }
        .set { ch_filtered_bam_index }

    ch_filtered_bam = ch_filtered_bam_index.map { meta, bam, index -> [meta, bam] }
    ch_filtered_index = ch_filtered_bam_index.map { meta, bam, index -> [meta, index] }

    // TODO: print for debugging
    ch_filtered_bam_index
        .map { meta, bam, index ->
            "${meta}\t${bam}\t${index}"
        }
        .collectFile(name: 'ch_filtered_bam_index_before_dSp.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/")

    if (params.bam_downsampling_method) {
        //
        // SUBWORKFLOW: Downsample IP and input control BAM files
        //
        BAM_DOWNSAMPLE (
            ch_filtered_bam_index,
            ch_fasta_fai,
            params.genome,
            params.spikein_genome,
            params.bam_downsampling_method,
            params.downsampling_endo_threshold,
            params.downsampling_exo_threshold
        )
        ch_filtered_bam = BAM_DOWNSAMPLE.out.bam
        ch_filtered_index = BAM_DOWNSAMPLE.out.index
        ch_filtered_bam_index = ch_filtered_bam.join(ch_filtered_index, by: 0)
        ch_samtools_stats_summary = ch_samtools_stats_summary.mix(BAM_DOWNSAMPLE.out.stats)
        ch_multiqc_files = ch_multiqc_files.mix(BAM_DOWNSAMPLE.out.stats.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(BAM_DOWNSAMPLE.out.flagstat.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(BAM_DOWNSAMPLE.out.idxstats.collect { it -> it[1] })

    }

    // TODO: print for debugging
    ch_filtered_bam_index
        .map { meta, bam, index ->
            "${meta}\t${bam}\t${index}"
        }
        .collectFile(name: 'ch_filtered_bam_index_after_dSp.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/")

    //
    // SUBWORKFLOW: Generate normalized bigWig coverage tracks
    //
    BAM_NORMALIZE_BIGWIG_DEEPTOOLS (
        ch_filtered_bam_index,
        ch_endo_chromsizes,
        ch_exo_chromsizes,
        params.genome,
        params.spikein_genome,
        params.min_reads_for_norm,
        params.skip_srpm,
        params.skip_cisrpm,
        params.skip_signal_vs_input,
        params.signal_vs_input_operation,
        params.skip_bw_average,
        params.skip_exo_bw,
        params.skip_multibigwigsummary,
        params.skip_plotcorrelation,
        params.skip_plotpca
    )
    ch_multiqc_files = ch_multiqc_files.mix(BAM_NORMALIZE_BIGWIG_DEEPTOOLS.out.multiqc_files)


    if (!params.skip_genes_plotprofile) {

        ch_bigwigs_genes = BAM_NORMALIZE_BIGWIG_DEEPTOOLS.out.bigwig_all_endo

        if (!params.input_cisrpm_in_plotprofile) {
            ch_bigwigs_genes
                .filter { meta, bws ->
                    !(meta.is_input_control && meta.norm_factor_type == 'cisrpm')
                }
                .set { ch_bigwigs_genes }
        }

       ch_bigwigs_genes
            .map { meta, bw ->
                def antibody = meta.antibody ?: meta.input_control_of_antibody
                [ antibody, meta.exp_type, meta.norm_factor_type, meta.signal_vs_input_operation, meta.averaged_brep, meta.id, meta, bw ]
            }
            .groupTuple(by: [0, 1, 2, 3, 4])
            .map { antibody, exp_type, norm_factor_type, signal_vs_input_op, averaged_brep, ids, metas, bws ->
                // Sort ids, metas and bws by id to ensure consistent order in plots
                def sorted_ids = ids.sort()
                def sorted_metas = metas.sort { meta -> meta.id }
                def sorted_bws = bws.sort { it -> it.name }
                [ antibody, exp_type, norm_factor_type, signal_vs_input_op, averaged_brep, sorted_ids, sorted_metas, sorted_bws ]
            }
            .combine(ch_gene_bed.map { it -> it[1] })
            .map {
                antibody, exp_type, norm_factor_type, signal_vs_input_op, averaged_brep, ids, metas, bws, gene_bed ->
                    def meta_new = metas[0].clone()
                    meta_new.id = exp_type + '_' +
                        (antibody ? antibody : 'no_antibody') +
                        '_' + norm_factor_type +
                        (signal_vs_input_op ? '_' + signal_vs_input_op : '') +
                        (averaged_brep ? '_bRep_avg' : '')
                    meta_new.antibody = antibody
                    meta_new.ids = ids
                    [ meta_new, bws.flatten(), gene_bed ]
            }
            .set { ch_bigwigs_genes }


        //
        // MODULE: deepTools matrix generation for plotting
        //
        DEEPTOOLS_COMPUTEMATRIX_GENES (
            ch_bigwigs_genes
        )

        //
        // MODULE: deepTools profile plots
        //
        DEEPTOOLS_PLOTPROFILE_GENES (
            DEEPTOOLS_COMPUTEMATRIX_GENES.out.matrix
        )
        ch_multiqc_files = ch_multiqc_files.mix(DEEPTOOLS_PLOTPROFILE_GENES.out.table.collect { it -> it[1] })

        //
        // MODULE: deepTools heatmaps
        //
        DEEPTOOLS_PLOTHEATMAP_GENES (
            DEEPTOOLS_COMPUTEMATRIX_GENES.out.matrix
        )
    }

    // Removing the exogenous samples from the filtered_bam_index channel
    ch_filtered_bam = ch_filtered_bam.filter { it -> it[0].genome == params.genome }
    ch_filtered_index = ch_filtered_index.filter { it -> it[0].genome == params.genome }
    ch_filtered_bam_index = ch_filtered_bam_index.filter { it -> it[0].genome == params.genome }

    //
    // SUBWORKFLOW: Counting reads in transposable elements
    //
    if (!params.skip_te_counting) {

        TE_COUNTING (
            // Here we run TE counting on both pre- and post flTbl and flT3 BAM files
            ch_filtered_bam.mix(ch_pre_flTbl_endo_bam).mix(ch_pre_flT3_endo_bam),
            ch_fasta_fai,
            false,
            ch_tecount_gene_index,
            ch_tecount_te_index,
            ch_telocal_gene_index,
            ch_telocal_te_index,
            params.skip_telocal,
            params.skip_tecount_gz,
            params.skip_telocal_gz
        )
    }

    BAM_NORMALIZE_BIGWIG_DEEPTOOLS
        .out
        .bedgraph_endo
        .filter { it -> !(it[0].signal_vs_input)}
        .set { ch_bedgraph_endo_for_seacr }


    ch_edd_peaks = channel.empty()
    ch_macs3_peaks = channel.empty()
    ch_genrich_peaks = channel.empty()
    ch_epic2_peaks = channel.empty()
    ch_mace_peaks = channel.empty()
    ch_consenrich_tracks = channel.empty()
    ch_rocco_peaks = channel.empty()
    ch_consensus_bed = channel.empty()
    ch_consensus_txt = channel.empty()

    if (!params.skip_peak_calling) {

        ch_peak_callers = channel.of( "${params.peak_callers}" ).flatMap { callers -> callers.tokenize(',') }

        //
        // SUBWORKFLOW: Call peaks
        //
        CALL_PEAKS (
            ch_filtered_bam_index,
            BAM_NORMALIZE_BIGWIG_DEEPTOOLS.out.bigwig_all_endo,
            ch_bedgraph_endo_for_seacr,
            ch_peak_callers,
            ch_fasta_fai,
            ch_gtf,
            ch_effective_gfraction,
            ch_endo_chromsizes,
            ch_blacklist.ifEmpty([[:], []]),
            ch_sparsebed.ifEmpty([[:], []]),
            ch_active_regions.ifEmpty([[:], []]),
            ch_rocco_params,
            ch_effective_gsize,
            ch_epic2_peak_count_header,
            ch_epic2_frip_score_header,
            ch_epic2_peak_annotation_header,
            ch_gr_peak_count_header,
            ch_gr_frip_score_header,
            ch_gr_peak_annotation_header,
            ch_mace_peak_count_header,
            ch_mace_frip_score_header,
            ch_mace_peak_annotation_header,
            ch_macs3_peak_count_header,
            ch_macs3_frip_score_header,
            ch_macs3_peak_annotation_header,
            ch_seacr_peak_count_header,
            ch_seacr_frip_score_header,
            ch_seacr_peak_annotation_header,
            ch_deseq2_pca_header,
            ch_deseq2_clustering_header,
            params.narrow_peak,
            params.skip_peak_annotation,
            params.skip_peak_qc,
            params.skip_bdgcmp,
            params.skip_consensus_peaks,
            params.skip_deseq2_qc,
            params.skip_consensus_plotprofile,
            params.input_cisrpm_in_plotprofile,
            params.seacr_peak_threshold
        )

        ch_edd_peaks = CALL_PEAKS.out.edd_peaks
        ch_macs3_peaks = CALL_PEAKS.out.macs3_peaks
        ch_genrich_peaks = CALL_PEAKS.out.genrich_peaks
        ch_epic2_peaks = CALL_PEAKS.out.epic2_peaks
        ch_mace_peaks = CALL_PEAKS.out.mace_peaks
        ch_consenrich_tracks = CALL_PEAKS.out.consenrich_tracks
        ch_rocco_peaks = CALL_PEAKS.out.rocco_peaks
        ch_consensus_bed = CALL_PEAKS.out.consensus_bed
        ch_consensus_txt = CALL_PEAKS.out.consensus_txt

        ch_multiqc_files = ch_multiqc_files.mix(CALL_PEAKS.out.multiqc_files)

    }

    //
    // Create channel for downstream processes: [ meta, [ ip_bam, ipcontrol_bam ] [ ip_index, ipcontrol_index ] ]
    // (Excluding ips_wo_ipcontrol as they don't need to be compared to anything)
    //
    ch_filtered_bam_index
        .branch { meta, bam, index ->
            ips_with_ipcontrol: meta.input_control
                return [meta.input_control, meta.antibody, meta, bam, index]
            ips_wo_ipcontrol: !meta.input_control && !meta.is_input_control
                return [meta, bam, index]
            ipcontrols: !meta.input_control && meta.is_input_control
                return [meta.id, meta.input_control_of_antibody, meta, bam, index]
        }
        .set { ch_bam_index_by_type }

    ch_bam_index_by_type
        .ips_with_ipcontrol
        .combine(ch_bam_index_by_type.ipcontrols, by: [0,1])
        .map { ipcontrol_id, antibody, ip_meta, ip_bam, ip_index, ipcontrol_meta, ipcontrol_bam, ipcontrol_index ->
            [ ip_meta, [ip_bam] + [ipcontrol_bam], [ip_index] + [ipcontrol_index] ]
        }
        .set { ch_ip_and_ipcontrols_bam_index }

    //
    // MODULE: deepTools plotFingerprint joint QC for IP and control
    //
    if (!params.skip_plot_fingerprint) {
        DEEPTOOLS_PLOTFINGERPRINT(
            ch_ip_and_ipcontrols_bam_index
        )
        ch_multiqc_files = ch_multiqc_files.mix(DEEPTOOLS_PLOTFINGERPRINT.out.matrix.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(DEEPTOOLS_PLOTFINGERPRINT.out.metrics.collect { it -> it[1] })
    }

    //
    // SUBWORKFLOW: Run ENCODE3 ChIP-seq pipeline
    //
    if (!params.skip_encode_pipeline) {
        BAM_ENCODE_PIPELINE (
            ch_filtered_bam,
            ch_fasta_fai,
            ch_endo_chromsizes,
            params.ctl_depth_ratio_threshold,
            params.narrow_peak ? 'narrowPeak' : 'broadPeak',
            ch_blacklist.ifEmpty([[:], []]),
            params.idr_filtering_threshold,
            params.encode_peak_max_score
        )

    }

    ch_filtered_bam_ss = channel.empty()
    ch_filtered_bam_ss = ch_filtered_bam.filter { it -> it[0].exp_type in ['SCAR-seq', 'OK-seq', 'eSPAN'] }

    // TODO: remove when optional inputs to subworkflows are implemented
    // Make ch_endo_chromsizes empty if there are no SCAR-seq samples
    // This is to avoid unnecessarily running modules in the BAM_CREATE_PARTITIONS
    ch_endo_chromsizes
        .combine(ch_filtered_bam_ss)
        .first()
        .map { sizes_meta, sizes, ss_meta, ss_bam ->
            [sizes_meta, sizes]
        }
        .set { ch_endo_chromsizes_ss }

    //
    // SUBWORKFLOW: SCAR-seq and OK-seq analysis: partitioning of reads
    //
    ch_partition_smooth = channel.empty()
    BAM_CREATE_PARTITIONS (
        ch_filtered_bam_ss,
        ch_fasta_fai,
        ch_endo_chromsizes_ss,
        ch_blacklist.ifEmpty([[:], []]),
        ch_okseq_rfd_file.ifEmpty([[:], [[]]]),
        ch_initiation_zones.ifEmpty([[:], [[]]]),
        params.smooth_radius,
        params.derivative_radius,
        params.zero_crossing_radius,
        params.skip_partition_group_plot
    )
    ch_partition_smooth = BAM_CREATE_PARTITIONS.out.tab

    //
    // SUBWORKFLOW: Repli-seq analysis: E/L ratio replication-timing (RT) tracks
    //

    ch_filtered_bam_index_repliseq = ch_filtered_bam_index.filter { it ->
        it[0].exp_type == 'Repli-seq' && it[0].rt_fraction in rtElFractions()
    }

    BAM_EL_REPLISEQ (
        ch_filtered_bam_index_repliseq,
        ch_endo_chromsizes,
        ch_blacklist.ifEmpty([[:], []]),
        ch_repliseq_rt_header,
        ch_repliseq_gene_class_header,
        ch_gene_bed
    )
    ch_multiqc_files = ch_multiqc_files.mix(BAM_EL_REPLISEQ.out.mqc.collect { it -> it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(BAM_EL_REPLISEQ.out.featurecounts_summary.collect { it -> it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(BAM_EL_REPLISEQ.out.gene_class_mqc.collect { it -> it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(BAM_EL_REPLISEQ.out.gene_class_box.collect { it -> it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(BAM_EL_REPLISEQ.out.domain_box.collect { it -> it[1] })

    //
    // SUBWORKFLOW: High-resolution (16-fraction) Repli-seq: array and replication features
    //
    ch_filtered_bam_index_hr_repliseq = ch_filtered_bam_index.filter { it ->
        it[0].exp_type == 'Repli-seq' && it[0].rt_fraction in rtHrFractions()
    }

    ch_el_track_for_hr = BAM_EL_REPLISEQ.out.bedgraph
        .filter { meta, _bedgraph ->
            meta.rt_measure == 'ratio' && meta.rt_track_type == 'smooth' && !meta.rt_covered
        }
        .map { meta, bedgraph -> [ meta.id, bedgraph ] }

    BAM_HR_REPLISEQ (
        ch_filtered_bam_index_hr_repliseq,
        ch_blacklist.ifEmpty([[:], []]),
        ch_hr_repliseq_features_header,
        ch_el_track_for_hr
    )
    ch_multiqc_files = ch_multiqc_files.mix(BAM_HR_REPLISEQ.out.mqc.collect { it -> it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(BAM_HR_REPLISEQ.out.mqc_heatmap.collect { it -> it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(BAM_HR_REPLISEQ.out.mqc_sizes.collect { it -> it[1] })

    //
    // SUBWORKFLOW: Create SAMtools summary table
    //
    SAMTOOLS_STATS_SUMMARY (
        ch_samtools_stats_summary,
        params.genome,
        params.spikein_genome ?: channel.value([])
    )

    //
    // MODULE: Create IGV session
    //
    ch_files_and_outpaths = channel.empty()
    if (!params.skip_igv) {

        BAM_NORMALIZE_BIGWIG_DEEPTOOLS.out.bigwig_endo
        .mix(BAM_NORMALIZE_BIGWIG_DEEPTOOLS.out.bigwig_cmp_endo)
        .mix(BAM_NORMALIZE_BIGWIG_DEEPTOOLS.out.bigwig_avg_endo)
            .map { meta, bw ->
                def outpath = "${params.outdir}/${params.aligner}/mergedLibrary/" +
                "${params.multimap_allocation_method ? params.multimap_allocation_method == 'chromap' ? 'cm_allo' : params.multimap_allocation_method : ''}" +
                "/${meta.exp_type}" +
                "${meta.downsampling_method ? '/downsampled' : ''}" +
                "/coverage/" +
                "/${meta.norm_factor_type}" +
                "/${meta.signal_vs_input_operation ? 'signal_vs_input/' + meta.signal_vs_input_operation : ''}" +
                "${meta.averaged_brep ? '/bRep_avg' : ''}" +
                "/${bw.getName()}"
                [meta, bw, outpath, "0,0,178"] // dark blue
            }
            .mix(ch_files_and_outpaths)
            .set { ch_files_and_outpaths }

        ch_edd_peaks
            .map { meta, peak ->
                def outpath = "${params.outdir}/${params.aligner}/mergedLibrary/" +
                    "${params.multimap_allocation_method ? params.multimap_allocation_method == 'chromap' ? 'cm_allo' : params.multimap_allocation_method : ''}" +
                    "/${meta.exp_type}" +
                    "${meta.downsampling_method ? '/downsampled' : ''}" +
                    "/edd/" +
                    "${peak.getName()}"
                [meta, peak, outpath, "255,140,0"] // dark orange
            }
            .mix(ch_files_and_outpaths)
            .set { ch_files_and_outpaths }

        ch_macs3_peaks
            .map { meta, peak ->
                def outpath = "${params.outdir}/${params.aligner}/mergedLibrary/" +
                    "${params.multimap_allocation_method ? params.multimap_allocation_method == 'chromap' ? 'cm_allo' : params.multimap_allocation_method : ''}" +
                    "/${meta.exp_type}" +
                    "${meta.downsampling_method ? '/downsampled' : ''}" +
                    "/macs3/" +
                    "${params.narrow_peak ? 'narrow_peak' : 'broad_peak'}" +
                    "/${peak.getName()}"
                [meta, peak, outpath, "178,34,34"] // firebrick
            }
            .mix(ch_files_and_outpaths)
            .set { ch_files_and_outpaths }

        ch_consensus_bed
            .mix(ch_consensus_txt)
            .map { meta, bed ->
                def outpath = "${params.outdir}/${params.aligner}/mergedLibrary/" +
                    "${params.multimap_allocation_method ? params.multimap_allocation_method == 'chromap' ? 'cm_allo' : params.multimap_allocation_method : ''}" +
                    "/${meta.exp_type}" +
                    "${meta.downsampling_method ? '/downsampled' : ''}" +
                    "/macs3/" +
                    "${params.narrow_peak ? 'narrow_peak' : 'broad_peak'}" +
                    "/consensus/${meta.id}" +
                    "${bed.getName()}"
                [meta, bed, outpath, "255,0,255"] // magenta
            }
            .mix(ch_files_and_outpaths)
            .set { ch_files_and_outpaths }

        ch_genrich_peaks
            .map { meta, peak ->
                def outpath = "${params.outdir}/${params.aligner}/mergedLibrary/" +
                    "${params.multimap_allocation_method ? params.multimap_allocation_method == 'chromap' ? 'cm_allo' : params.multimap_allocation_method : ''}" +
                    "/${meta.exp_type}" +
                    "${meta.downsampling_method ? '/downsampled' : ''}" +
                    "/genrich/" +
                    "${params.narrow_peak ? 'narrow_peak' : 'broad_peak'}" +
                    "/${peak.getName()}"
                [meta, peak, outpath, "34,139,34"] // forest green
            }
            .mix(ch_files_and_outpaths)
            .set { ch_files_and_outpaths }

        ch_mace_peaks
            .map { meta, peak ->
                def outpath = "${params.outdir}/${params.aligner}/mergedLibrary/" +
                    "${params.multimap_allocation_method ? params.multimap_allocation_method == 'chromap' ? 'cm_allo' : params.multimap_allocation_method : ''}" +
                    "/${meta.exp_type}" +
                    "${meta.downsampling_method ? '/downsampled' : ''}" +
                    "/mace/" +
                    "${peak.getName()}"
                [meta, peak, outpath, "255,140,0"] // dark orange
            }
            .mix(ch_files_and_outpaths)
            .set { ch_files_and_outpaths }

        ch_epic2_peaks
            .map { meta, peak ->
                def outpath = "${params.outdir}/${params.aligner}/mergedLibrary/" +
                    "${params.multimap_allocation_method ? params.multimap_allocation_method == 'chromap' ? 'cm_allo' : params.multimap_allocation_method : ''}" +
                    "/${meta.exp_type}" +
                    "${meta.downsampling_method ? '/downsampled' : ''}" +
                    "/epic2/" +
                    "${peak.getName()}"
                [meta, peak, outpath, "138,43,226"] // blue violet
            }
            .mix(ch_files_and_outpaths)
            .set { ch_files_and_outpaths }

        ch_consenrich_tracks
            .map { meta, signal ->
                def outpath = "${params.outdir}/${params.aligner}/mergedLibrary/" +
                    "${params.multimap_allocation_method ? params.multimap_allocation_method == 'chromap' ? 'cm_allo' : params.multimap_allocation_method : ''}" +
                    "/${meta.exp_type}" +
                    "${meta.downsampling_method ? '/downsampled' : ''}" +
                    "/consenrich/" +
                    "${signal.getName()}"
                [meta, signal, outpath, "255,20,147"] // deep pink
            }
            .mix(ch_files_and_outpaths)
            .set { ch_files_and_outpaths }

        ch_rocco_peaks
            .map { meta, peak ->
                def outpath = "${params.outdir}/${params.aligner}/mergedLibrary/" +
                    "${params.multimap_allocation_method ? params.multimap_allocation_method == 'chromap' ? 'cm_allo' : params.multimap_allocation_method : ''}" +
                    "/${meta.exp_type}" +
                    "${meta.downsampling_method ? '/downsampled' : ''}" +
                    "/consenrich/rocco/" +
                    "${peak.getName()}"
                [meta, peak, outpath, "255,105,180"] // hot pink
            }
            .mix(ch_files_and_outpaths)
            .set { ch_files_and_outpaths }

        ch_gtf
            .map { meta, gtf ->
                def outpath = "${params.outdir}/genome/" +
                    "${gtf.getName()}"
                [meta, gtf, outpath, "0,128,0"] // green
            }
            .mix(ch_files_and_outpaths)
            .set { ch_files_and_outpaths }

        BAM_EL_REPLISEQ.out.bigwig
            .map { meta, bw ->
                def measure_label = meta.rt_measure == 'rt_index' ? 'RT_index' :
                    (meta.rt_measure == 'coverage' ? 'coverage' :
                    (params.repliseq_ratio_direction == 'late_over_early' ? 'LE_ratio' : 'EL_ratio'))
                def outpath = "${params.outdir}/${params.aligner}/mergedLibrary/" +
                    "${params.multimap_allocation_method ? (params.multimap_allocation_method != 'chromap' ? params.multimap_allocation_method : '') : ''}" +
                    "/${meta.exp_type}" +
                    "/EL_repliseq/bigwig/${measure_label}/${meta.rt_track_type}/" +
                    "${bw.getName()}"
                def color = meta.rt_measure == 'rt_index' ?
                    (meta.rt_track_type == 'smooth' ? "128,0,128" : "216,191,216") : // purple / thistle
                    meta.rt_measure == 'coverage' ? "105,105,105" : // dim grey
                    (meta.rt_track_type == 'smooth' ? "0,100,0" : "144,238,144") // dark green / light green
                [meta, bw, outpath, color]
            }
            .mix(ch_files_and_outpaths)
            .set { ch_files_and_outpaths }

        BAM_EL_REPLISEQ.out.domains
            .map { meta, bed ->
                def outpath = "${params.outdir}/${params.aligner}/mergedLibrary/" +
                    "${params.multimap_allocation_method ? (params.multimap_allocation_method != 'chromap' ? params.multimap_allocation_method : '') : ''}" +
                    "/${meta.exp_type}" +
                    "/EL_repliseq/domains/" +
                    "${bed.getName()}"
                [meta, bed, outpath, "47,79,79"] // dark slate grey
            }
            .mix(ch_files_and_outpaths)
            .set { ch_files_and_outpaths }

        BAM_EL_REPLISEQ.out.gene_class_bed
            .map { meta, bed ->
                def outpath = "${params.outdir}/${params.aligner}/mergedLibrary/" +
                    "${params.multimap_allocation_method ? (params.multimap_allocation_method != 'chromap' ? params.multimap_allocation_method : '') : ''}" +
                    "/${meta.exp_type}" +
                    "/EL_repliseq/gene_classification/" +
                    "${bed.getName()}"
                def color = [
                    early: "33,102,172", // #2166AC
                    mid: "178,171,210",  // #B2ABD2
                    late: "178,24,43",   // #B2182B
                    unclassified: "128,128,128"
                ][meta.rt_gene_class]
                [meta, bed, outpath, color]
            }
            .mix(ch_files_and_outpaths)
            .set { ch_files_and_outpaths }

        BAM_HR_REPLISEQ.out.iz.map { meta, bed -> [meta, bed, "220,20,60"] } // crimson
            .mix(BAM_HR_REPLISEQ.out.ttr.map { meta, bed -> [meta, bed, "255,165,0"] }) // orange
            .mix(BAM_HR_REPLISEQ.out.breakage.map { meta, bed -> [meta, bed, "139,0,0"] }) // dark red
            .mix(BAM_HR_REPLISEQ.out.termination.map { meta, bed -> [meta, bed, "0,0,139"] }) // dark blue
            .mix(BAM_HR_REPLISEQ.out.ctr.map { meta, bed -> [meta, bed, "0,139,139"] }) // dark cyan
            .mix(BAM_HR_REPLISEQ.out.partition.map { meta, bed -> [meta, bed, "85,107,47"] }) // dark olive green
            .map { meta, bed, color ->
                def outpath = "${params.outdir}/${params.aligner}/mergedLibrary/" +
                    "${params.multimap_allocation_method ? (params.multimap_allocation_method != 'chromap' ? params.multimap_allocation_method : '') : ''}" +
                    "/${meta.exp_type}" +
                    "/hr_repliseq/features/" +
                    "${bed.getName()}"
                [meta, bed, outpath, color]
            }
            .mix(ch_files_and_outpaths)
            .set { ch_files_and_outpaths }

        // create channel: [ list_of_files, list_of_outpaths ]
        ch_files_and_outpaths
            .map { meta, file, outpath, color -> [1, file, outpath, color]}
            .groupTuple()
            .map { id, files, outpaths, colors -> [files, outpaths, colors]}
            .set { ch_files_and_outpaths }

        ch_fasta
            .map { meta, fasta ->
                def outpath = "${params.outdir}/genome/" +
                    "${fasta.getName()}"
                [fasta, outpath]
            }
            .set { ch_fasta_outpath }

        IGV (
            ch_files_and_outpaths,
            ch_fasta_outpath
        )
    }

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'grothlab' + 'crepas_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    if (!params.skip_multiqc) {

        // Prepare the workflow summary
        ch_workflow_summary = channel.value(
                paramsSummaryMultiqc(
                    paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
                )
            )
            .collectFile(name: 'workflow_summary_mqc.yaml')

        // Prepare the methods section
        ch_methods_description = channel.value(
                methodsDescriptionText(
                    params.multiqc_methods_description
                        ? file(params.multiqc_methods_description, checkIfExists: true)
                        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
                )
            )
            .collectFile(name: 'methods_description_mqc.yaml', sort: true)

        // Add summary, versions, and methods to the MultiQC input file list
        ch_multiqc_files = ch_multiqc_files
            .mix(ch_workflow_summary)
            .mix(ch_collated_versions)
            .mix(ch_methods_description)

        //
        // MODULE: MultiQC
        //
        MULTIQC (
            ch_multiqc_files.flatten().collect().map { files ->
                [
                    [id: 'crepas'],
                    files,
                    params.multiqc_config
                        ? file(params.multiqc_config, checkIfExists: true)
                        : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                    params.multiqc_logo ? file(params.multiqc_logo, checkIfExists: true) : [],
                    [],
                    [],
                ]
            }
        )

        ch_multiqc_report  = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
        ch_multiqc_publish = MULTIQC.out.data.mix(MULTIQC.out.plots, MULTIQC.out.report)
        ch_versions = ch_versions.mix(MULTIQC.out.versions)

    }

    emit:
    multiqc_report  = ch_multiqc_report
    multiqc_publish = ch_multiqc_publish
    versions        = ch_versions // channel: [ path(versions.yml) ]
}
