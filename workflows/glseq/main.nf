/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { IGV                                                         } from '../../modules/local/igv/main'
include { MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS                         } from '../../modules/local/multiqc_custom_phantompeakqualtools/main'
include {
    BAM_FLAGSTAT_MAPPED as BAM_FLAGSTAT_MAPPED_FLT1 ;
    BAM_FLAGSTAT_MAPPED as BAM_FLAGSTAT_MAPPED_FLT2 ;
    BAM_FLAGSTAT_MAPPED as BAM_FLAGSTAT_MAPPED_FLT3
} from '../../modules/local/bam_flagstat_mapped/main'
include { EDD } from '../../modules/local/edd/main'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { paramsSummaryMap                                            } from 'plugin/nf-validation'
include { paramsSummaryMultiqc                                        } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                                      } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                                      } from '../../subworkflows/local/utils_grothlab_glseq_pipeline'
include { INPUT_CHECK                                                 } from '../../subworkflows/local/input_check/main'
include {
    BAM_FILTER_SAMBAMBA as BAM_FILTER_SAMBAMBA_FLT1 ;
    BAM_FILTER_SAMBAMBA as BAM_FILTER_SAMBAMBA_FLT3
} from '../../subworkflows/local/bam_filter_sambamba/main'
include { BAM_SPIKEIN_SPLIT                                           } from '../../subworkflows/local/bam_spikein_split/main'
include { FASTQ_FASTQC_UMITOOLS_UMITRANSFER_TRIMGALORE                } from '../../subworkflows/local/fastq_fastqc_umitools_umitransfer_trimgalore/main'
include { BAM_PEAKS_CALL_QC_ANNOTATE_EPIC2_HOMER                      } from '../../subworkflows/local/bam_peaks_call_qc_annotate_epic2_homer/main'
include { BAM_PEAKS_CALL_QC_ANNOTATE_MACS3_HOMER                      } from '../../subworkflows/local/bam_peaks_call_qc_annotate_macs3_homer/main'
include { BAM_PEAKS_CALL_QC_ANNOTATE_GENRICH_HOMER                    } from '../../subworkflows/local/bam_peaks_call_qc_annotate_genrich_homer/main'
include { BAM_PEAKS_CALL_QC_ANNOTATE_MACE_HOMER                       } from '../../subworkflows/local/bam_peaks_call_qc_annotate_mace_homer/main'
include { BED_CONSENSUS_QUANTIFY_QC_BEDTOOLS_FEATURECOUNTS_DESEQ2     } from '../../subworkflows/local/bed_consensus_quantify_qc_bedtools_featurecounts_deseq2/main'
include { BAM_CREATE_PARTITIONS                                       } from '../../subworkflows/local/bam_create_partitions/main'
include { BAM_ALLOCATE_MULTIMAPPERS as BAM_ALLOCATE_MULTIMAPPERS_ENDO } from '../../subworkflows/local/bam_allocate_multimappers/main'
include { BAM_ALLOCATE_MULTIMAPPERS as BAM_ALLOCATE_MULTIMAPPERS_EXO  } from '../../subworkflows/local/bam_allocate_multimappers/main'
include { BAM_PEAKS_CALL_QC_ANNOTATE_CONSENRICH_HOMER                 } from '../../subworkflows/local/bam_peaks_call_qc_annotate_consenrich_homer/main'
include { BAM_SHIFT_READS                                             } from '../../subworkflows/local/bam_shift_reads/main'
include { SAMTOOLS_STATS_SUMMARY                                      } from '../../subworkflows/local/samtools_stats_summary/main'
include { BAM_FILTER_BLACKLIST                                        } from '../../subworkflows/local/bam_filter_blacklist/main'
include { BAM_NORMALIZE_BIGWIG_DEEPTOOLS                              } from '../../subworkflows/local/bam_normalize_bigwig_deeptools/main'
include { BAM_DOWNSAMPLE                                              } from '../../subworkflows/local/bam_downsample/main'
include { TE_COUNTING                                                 } from '../../subworkflows/local/te_counting/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { SAMTOOLS_INDEX                                              } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_COLLATE                                            } from '../../modules/nf-core/samtools/collate/main'
include { PICARD_MERGESAMFILES                                        } from '../../modules/nf-core/picard/mergesamfiles/main'
include { PICARD_COLLECTMULTIPLEMETRICS                               } from '../../modules/nf-core/picard/collectmultiplemetrics/main'
include { PRESEQ_LCEXTRAP                                             } from '../../modules/nf-core/preseq/lcextrap/main'
include { PHANTOMPEAKQUALTOOLS                                        } from '../../modules/nf-core/phantompeakqualtools/main'
include { DEEPTOOLS_COMPUTEMATRIX                                     } from '../../modules/nf-core/deeptools/computematrix/main'
include { DEEPTOOLS_PLOTPROFILE                                       } from '../../modules/nf-core/deeptools/plotprofile/main'
include { DEEPTOOLS_PLOTHEATMAP                                       } from '../../modules/nf-core/deeptools/plotheatmap/main'
include { DEEPTOOLS_PLOTFINGERPRINT                                   } from '../../modules/nf-core/deeptools/plotfingerprint/main'
include { MULTIQC                                                     } from '../../modules/nf-core/multiqc/main'

//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//

// include { FASTQ_FASTQC_UMITOOLS_TRIMGALORE      } from '../../subworkflows/nf-core/fastq_fastqc_umitools_trimgalore'
include { FASTQ_ALIGN_BWA                                             } from '../../subworkflows/nf-core/fastq_align_bwa'
include { FASTQ_ALIGN_BOWTIE2                                         } from '../../subworkflows/nf-core/fastq_align_bowtie2'
include { FASTQ_ALIGN_CHROMAP                                         } from '../../subworkflows/nf-core/fastq_align_chromap'
include { FASTQ_ALIGN_STAR                                            } from '../../subworkflows/nf-core/fastq_align_star'
include { FASTQ_ALIGN_HISAT2                                          } from '../../subworkflows/nf-core/fastq_align_hisat2'
include { BAM_MARKDUPLICATES_PICARD                                   } from '../../subworkflows/nf-core/bam_markduplicates_picard'
include { BAM_DEDUP_UMI                                               } from '../../subworkflows/nf-core/bam_dedup_umi'
include { BAM_STATS_SAMTOOLS                                          } from '../../subworkflows/nf-core/bam_stats_samtools'
include { BAM_SORT_STATS_SAMTOOLS                                     } from '../../subworkflows/nf-core/bam_sort_stats_samtools'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GLSEQ {
    take:
    ch_input                  // channel: path(sample_sheet.csv)
    ch_versions               // channel: [ path(versions.yml) ]
    ch_fasta                  // channel: path(genome.fa)
    ch_fai                    // channel: path(genome.fai)
    ch_gtf                    // channel: path(genome.gtf)
    ch_gene_bed               // channel: path(gene.beds)
    ch_chrom_sizes_endo       // path(chrom.sizes.endo)
    ch_chrom_sizes_exo
    ch_effective_gsize        
    ch_effective_gfraction
    ch_whitelist           // channel: path(filtered.bed)
    ch_blacklist              // channel: path(blacklist.bed)
    ch_sparsebed              // channel: path(sparse.bed)
    ch_active_regions         // channel: path(active_regions.bed)
    ch_rocco_params           // channel: path(params.csv)
    ch_initiation_zones       // channel: path(initiation_zones)
    ch_bwa_index              // channel: path(bwa/index/)
    ch_bowtie2_index          // channel: path(bowtie2/index)
    ch_chromap_index          // channel: path(chromap.index)
    ch_star_index             // channel: path(star/index/)
    ch_hisat2_index           // channel: path(hisat2/index)
    ch_splicesites            // channel: path(splicesites)
    ch_tecount_gene_index // channel: val(meta), path(tecount_gene_index.Ind)
    ch_telocal_gene_index // channel: val(meta), path(telocal_gene_index.Ind)
    ch_tecount_te_index       // channel: val(meta), path(tecount_te_index.Ind)
    ch_telocal_te_index       // channel: val(meta), path(telocal_te_index.locInd)

    main:
    ch_multiqc_files = Channel.empty()
    ch_samtools_stats_summary = Channel.empty()

    // TODO: organize these:

    // Header files for MultiQC
    ch_spp_nsc_header = file("${projectDir}/assets/multiqc/spp_nsc_header.txt", checkIfExists: true)
    ch_spp_rsc_header = file("${projectDir}/assets/multiqc/spp_rsc_header.txt", checkIfExists: true)
    ch_spp_correlation_header = file("${projectDir}/assets/multiqc/spp_correlation_header.txt", checkIfExists: true)
    ch_peak_count_header = file("${projectDir}/assets/multiqc/peak_count_header.txt", checkIfExists: true)
    ch_gr_peak_count_header = file("${projectDir}/assets/multiqc/gr_peak_count_header.txt", checkIfExists: true)
    ch_mace_peak_count_header = file("${projectDir}/assets/multiqc/mace_peak_count_header.txt", checkIfExists: true)
    ch_epic2_peak_count_header = file("${projectDir}/assets/multiqc/epic2_peak_count_header.txt", checkIfExists: true)
    ch_frip_score_header = file("${projectDir}/assets/multiqc/frip_score_header.txt", checkIfExists: true)
    ch_gr_frip_score_header = file("${projectDir}/assets/multiqc/gr_frip_score_header.txt", checkIfExists: true)
    ch_mace_frip_score_header = file("${projectDir}/assets/multiqc/mace_frip_score_header.txt", checkIfExists: true)
    ch_epic2_frip_score_header = file("${projectDir}/assets/multiqc/epic2_frip_score_header.txt", checkIfExists: true)
    ch_peak_annotation_header = file("${projectDir}/assets/multiqc/peak_annotation_header.txt", checkIfExists: true)
    ch_gr_peak_annotation_header = file("${projectDir}/assets/multiqc/gr_peak_annotation_header.txt", checkIfExists: true)
    ch_mace_peak_annotation_header = file("${projectDir}/assets/multiqc/mace_peak_annotation_header.txt", checkIfExists: true)
    ch_epic2_peak_annotation_header = file("${projectDir}/assets/multiqc/epic2_peak_annotation_header.txt", checkIfExists: true)
    ch_deseq2_pca_header = Channel.value(file("${projectDir}/assets/multiqc/deseq2_pca_header.txt", checkIfExists: true))
    ch_deseq2_clustering_header = Channel.value(file("${projectDir}/assets/multiqc/deseq2_clustering_header.txt", checkIfExists: true))

    // Save AWS IGenomes file containing annotation version
    def anno_readme = params.genomes[params.genome]?.readme
    if (anno_readme && file(anno_readme).exists()) {
        file("${params.outdir}/genome/").mkdirs()
        file(anno_readme).copyTo("${params.outdir}/genome/")
    }

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK(
        ch_input,
        params.seq_center
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // SUBWORKFLOW: Read QC and trim adapters
    //
    FASTQ_FASTQC_UMITOOLS_UMITRANSFER_TRIMGALORE(
        INPUT_CHECK.out.reads,
        params.skip_fastqc || params.skip_qc,
        params.with_umi,
        params.skip_umi_extract,
        params.skip_trimming,
        params.umi_discard_read,
        params.min_trimmed_reads,
        params.hardtrim5_length,
        params.hardtrim3_length
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_UMITRANSFER_TRIMGALORE.out.fastqc_zip.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_UMITRANSFER_TRIMGALORE.out.trim_zip.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_UMITRANSFER_TRIMGALORE.out.trim_log.collect { it[1] })
    ch_versions = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_UMITRANSFER_TRIMGALORE.out.versions)

    //
    // SUBWORKFLOW: Alignment with BWA & BAM QC
    //
    ch_genome_bam = Channel.empty()
    ch_genome_bam_index = Channel.empty()
    if (params.aligner == 'bwa') {
        FASTQ_ALIGN_BWA(
            FASTQ_FASTQC_UMITOOLS_UMITRANSFER_TRIMGALORE.out.reads,
            ch_bwa_index,
            params.sort_bam,
            ch_fasta
        )
        ch_genome_bam = FASTQ_ALIGN_BWA.out.bam
        ch_genome_bam_index = FASTQ_ALIGN_BWA.out.bai
        ch_samtools_stats_summary = ch_samtools_stats_summary.mix(FASTQ_ALIGN_BWA.out.stats)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_BWA.out.stats.collect { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_BWA.out.flagstat.collect { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_BWA.out.idxstats.collect { it[1] })
        ch_versions = ch_versions.mix(FASTQ_ALIGN_BWA.out.versions.first())
    }

    //
    // SUBWORKFLOW: Alignment with Bowtie2 & BAM QC
    //
    if (params.aligner == 'bowtie2') {
        FASTQ_ALIGN_BOWTIE2(
            FASTQ_FASTQC_UMITOOLS_UMITRANSFER_TRIMGALORE.out.reads,
            ch_bowtie2_index,
            params.save_unaligned,
            params.sort_bam,
            ch_fasta
        )
        ch_genome_bam = FASTQ_ALIGN_BOWTIE2.out.bam
        ch_genome_bam_index = FASTQ_ALIGN_BOWTIE2.out.bai
        ch_samtools_stats_summary = ch_samtools_stats_summary.mix(FASTQ_ALIGN_BOWTIE2.out.stats)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_BOWTIE2.out.stats.collect { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_BOWTIE2.out.flagstat.collect { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_BOWTIE2.out.idxstats.collect { it[1] })
        ch_versions = ch_versions.mix(FASTQ_ALIGN_BOWTIE2.out.versions.first())
    }

    //
    // SUBWORKFLOW: Alignment with Chromap & BAM QC
    //
    if (params.aligner == 'chromap') {
        FASTQ_ALIGN_CHROMAP(
            FASTQ_FASTQC_UMITOOLS_UMITRANSFER_TRIMGALORE.out.reads,
            ch_chromap_index,
            ch_fasta,
            [],
            [],
            [],
            []
        )

        ch_genome_bam = FASTQ_ALIGN_CHROMAP.out.bam
        ch_genome_bam_index = FASTQ_ALIGN_CHROMAP.out.bai
        ch_samtools_stats_summary = ch_samtools_stats_summary.mix(FASTQ_ALIGN_CHROMAP.out.stats)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_CHROMAP.out.stats.collect { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_CHROMAP.out.flagstat.collect { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_CHROMAP.out.idxstats.collect { it[1] })
        ch_versions = ch_versions.mix(FASTQ_ALIGN_CHROMAP.out.versions.first())
    }

    //
    // SUBWORKFLOW: Alignment with STAR & BAM QC
    //
    if (params.aligner == 'star') {
        FASTQ_ALIGN_STAR(
            FASTQ_FASTQC_UMITOOLS_UMITRANSFER_TRIMGALORE.out.reads,
            ch_star_index,
            ch_gtf,
            true,
            params.seq_platform ?: '',
            params.seq_center ?: '',
            ch_fasta,
            Channel.value([[:], []])
        )
        ch_genome_bam = FASTQ_ALIGN_STAR.out.bam
        ch_genome_bam_index = FASTQ_ALIGN_STAR.out.bai
        ch_transcriptome_bam = FASTQ_ALIGN_STAR.out.bam_transcript
        ch_samtools_stats_summary = ch_samtools_stats_summary.mix(FASTQ_ALIGN_STAR.out.stats)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_STAR.out.stats.collect { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_STAR.out.flagstat.collect { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_STAR.out.idxstats.collect { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_STAR.out.log_final.collect { it[1] })
        ch_versions = ch_versions.mix(FASTQ_ALIGN_STAR.out.versions)
    }

    if (params.aligner == 'hisat2') {
        FASTQ_ALIGN_HISAT2(
            FASTQ_FASTQC_UMITOOLS_UMITRANSFER_TRIMGALORE.out.reads,
            ch_hisat2_index,
            ch_splicesites,
            ch_fasta
        )
        ch_genome_bam = FASTQ_ALIGN_HISAT2.out.bam
        ch_genome_bam_index = FASTQ_ALIGN_HISAT2.out.bai
        ch_samtools_stats_summary = ch_samtools_stats_summary.mix(FASTQ_ALIGN_HISAT2.out.stats)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_HISAT2.out.stats.collect { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_HISAT2.out.flagstat.collect { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_HISAT2.out.idxstats.collect { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_HISAT2.out.summary.collect { it[1] })
        ch_versions = ch_versions.mix(FASTQ_ALIGN_HISAT2.out.versions)
    }

    //
    // MODULE: Merge resequenced BAM files
    //
    ch_genome_bam
        .map { meta, bam ->
            def meta_clone = meta.clone()
            meta_clone.remove('read_group')
            meta_clone.id = meta_clone.id.split('_')[0..-2].join('_')
            [meta_clone, bam]
        }
        .groupTuple(by: 0)
        .map { it ->
            [it[0], it[1].flatten()]
        }
        .set { ch_sort_bam }

    PICARD_MERGESAMFILES(
        ch_sort_bam
    )
    ch_merged_bam = PICARD_MERGESAMFILES.out.bam
    ch_versions = ch_versions.mix(PICARD_MERGESAMFILES.out.versions.first())

    SAMTOOLS_INDEX(
        ch_merged_bam
    )
    ch_merged_bam_bai = ch_merged_bam.join(SAMTOOLS_INDEX.out.bai, by: [0])
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    BAM_STATS_SAMTOOLS(
        ch_merged_bam_bai,
        ch_fasta
    )
    ch_samtools_stats_summary = ch_samtools_stats_summary.mix(BAM_STATS_SAMTOOLS.out.stats)
    ch_multiqc_files = ch_multiqc_files.mix(BAM_STATS_SAMTOOLS.out.stats.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(BAM_STATS_SAMTOOLS.out.flagstat.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(BAM_STATS_SAMTOOLS.out.idxstats.collect { it[1] })
    ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)


    // TODO: change this so UMI dedup is evaluated per sample and not for the whole pipeline run
    if (params.with_umi) {

        //
        // MODULE: Preseq coverage analysis
        //
        // TODO: this is done on the bams with spike-in included
        if (!params.skip_preseq) {
            PRESEQ_LCEXTRAP(
                ch_merged_bam
            )
            ch_multiqc_files = ch_multiqc_files.mix(PRESEQ_LCEXTRAP.out.lc_extrap.collect { it[1] })
            ch_versions = ch_versions.mix(PRESEQ_LCEXTRAP.out.versions.first())
        }

        //
        // SUBWORKFLOW: Deduplicate BAM files
        //
        ch_transcriptome_bam = Channel.empty()
        ch_transcriptome_fasta = Channel.empty()
        BAM_DEDUP_UMI(
            ch_merged_bam_bai,
            [],
            params.umi_dedup_tool,
            params.get_dedup_stats,
            false,
            ch_transcriptome_bam,
            ch_transcriptome_fasta
        )
        ch_dedup_bam = BAM_DEDUP_UMI.out.bam
        ch_dedup_index = BAM_DEDUP_UMI.out.bai
        ch_dedup_umi_stats = BAM_DEDUP_UMI.out.stats
        ch_dedup_umi_flagstat = BAM_DEDUP_UMI.out.flagstat
        ch_dedup_umi_idxstats = BAM_DEDUP_UMI.out.idxstats
        ch_multiqc_files = ch_multiqc_files.mix(ch_dedup_umi_stats.collect { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(ch_dedup_umi_flagstat.collect { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(ch_dedup_umi_idxstats.collect { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(BAM_DEDUP_UMI.out.dedup_log.collect { it[1] })
        ch_versions = ch_versions.mix(BAM_DEDUP_UMI.out.versions)
    }
    else {
        //
        // SUBWORKFLOW: Mark duplicates & filter BAM files
        //
        BAM_MARKDUPLICATES_PICARD(
            ch_merged_bam,
            ch_fasta,
            ch_fai
        )
        ch_dedup_bam = BAM_MARKDUPLICATES_PICARD.out.bam
        ch_dedup_index = BAM_MARKDUPLICATES_PICARD.out.bai
        ch_samtools_stats_summary = ch_samtools_stats_summary.mix(BAM_MARKDUPLICATES_PICARD.out.stats)
        ch_multiqc_files = ch_multiqc_files.mix(BAM_MARKDUPLICATES_PICARD.out.stats.collect { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(BAM_MARKDUPLICATES_PICARD.out.flagstat.collect { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(BAM_MARKDUPLICATES_PICARD.out.idxstats.collect { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(BAM_MARKDUPLICATES_PICARD.out.metrics.collect { it[1] })
        ch_versions = ch_versions.mix(BAM_MARKDUPLICATES_PICARD.out.versions)

        //
        // MODULE: Preseq coverage analysis
        //
        // TODO: this is done on the bams with spike-in included
        if (!params.skip_preseq) {
            PRESEQ_LCEXTRAP(
                ch_dedup_bam
            )
            ch_multiqc_files = ch_multiqc_files.mix(PRESEQ_LCEXTRAP.out.lc_extrap.collect { it[1] })
            ch_versions = ch_versions.mix(PRESEQ_LCEXTRAP.out.versions.first())
        }
    }


    //
    // SUBWORKFLOW: Filter BAM file with SAMBAMBA
    //
    BAM_FILTER_SAMBAMBA_FLT1(
        ch_dedup_bam.join(ch_dedup_index, by: 0),
        Channel.value([[:], []]),
        ch_fasta
    )
    ch_filtered_bam = BAM_FILTER_SAMBAMBA_FLT1.out.bam
    ch_filtered_index = BAM_FILTER_SAMBAMBA_FLT1.out.bai
    ch_samtools_stats_summary = ch_samtools_stats_summary.mix(BAM_FILTER_SAMBAMBA_FLT1.out.stats)
    ch_multiqc_files = ch_multiqc_files.mix(BAM_FILTER_SAMBAMBA_FLT1.out.stats.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(BAM_FILTER_SAMBAMBA_FLT1.out.flagstat.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(BAM_FILTER_SAMBAMBA_FLT1.out.idxstats.collect { it[1] })
    ch_versions = ch_versions.mix(BAM_FILTER_SAMBAMBA_FLT1.out.versions)

    //
    // MODULE: Extract total mapped reads from flagstats
    //
    BAM_FLAGSTAT_MAPPED_FLT1(
        BAM_FILTER_SAMBAMBA_FLT1.out.flagstat
    )
    ch_versions = ch_versions.mix(BAM_FLAGSTAT_MAPPED_FLT1.out.versions)

    // Extract the total mapped reads from the text file
    BAM_FLAGSTAT_MAPPED_FLT1.out.txt
        .map { meta, total ->
            [meta, total.splitCsv(header: false)[0][0]]
        }
        .set { ch_flT1_total }

    // Add the total_mapped_reads to the bams' and bais' metas
    ch_filtered_bam
        .join(ch_filtered_index, by: 0)
        .combine(ch_flT1_total, by: 0)
        .map { meta, bam, bai, total ->
            def meta_clone = meta.clone()
            meta_clone.flT1_total_mapped_reads = total.toDouble()
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

    //
    // SUBWORKFLOW: Spike-in splitting
    //
    // TODO: if fasta and gtf are specified but not genome, val keep_genome_string in
    // BAM_SPLIT_BY_GENOME will fail
    ch_filtered_exo_bam = Channel.empty()
    ch_filtered_exo_index = Channel.empty()
    if (params.spikein_genome) {
        BAM_SPIKEIN_SPLIT(
            ch_filtered_bam,
            ch_fasta,
            Channel.value([[:], []]),
            params.genome,
            params.spikein_genome
        )
        ch_filtered_bam = BAM_SPIKEIN_SPLIT.out.bam
        ch_filtered_exo_bam = BAM_SPIKEIN_SPLIT.out.exo_bam
        ch_filtered_index = BAM_SPIKEIN_SPLIT.out.bai
        ch_filtered_exo_index = BAM_SPIKEIN_SPLIT.out.exo_bai
        ch_samtools_stats_summary = ch_samtools_stats_summary.mix(BAM_SPIKEIN_SPLIT.out.stats)
        ch_samtools_stats_summary = ch_samtools_stats_summary.mix(BAM_SPIKEIN_SPLIT.out.exo_stats)
        ch_multiqc_files = ch_multiqc_files.mix(BAM_SPIKEIN_SPLIT.out.stats.collect { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(BAM_SPIKEIN_SPLIT.out.exo_stats.collect { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(BAM_SPIKEIN_SPLIT.out.flagstat.collect { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(BAM_SPIKEIN_SPLIT.out.exo_flagstat.collect { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(BAM_SPIKEIN_SPLIT.out.idxstats.collect { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(BAM_SPIKEIN_SPLIT.out.exo_idxstats.collect { it[1] })
        ch_versions = ch_versions.mix(BAM_SPIKEIN_SPLIT.out.versions.first())

        //
        // MODULE: Extract total mapped reads from flagstats
        //
        BAM_FLAGSTAT_MAPPED_FLT2(
            BAM_SPIKEIN_SPLIT.out.flagstat.mix(BAM_SPIKEIN_SPLIT.out.exo_flagstat)
        )
        ch_versions = ch_versions.mix(BAM_FLAGSTAT_MAPPED_FLT2.out.versions)

        // Extract the total mapped reads from the text file
        BAM_FLAGSTAT_MAPPED_FLT2.out.txt
            .map { meta, total ->
                [meta, total.splitCsv(header: false)[0][0]]
            }
            .set { ch_flT2_total }

        // Add the total_mapped_reads both endo and exo bams' and bais' metas
        ch_filtered_bam
            .mix(ch_filtered_exo_bam)
            .join(ch_filtered_index.mix(ch_filtered_exo_index), by: 0)
            .combine(ch_flT2_total, by: 0)
            .map { meta, bam, bai, total ->
                def meta_clone = meta.clone()
                meta_clone.flT2_total_mapped_reads = total.toDouble()
                [meta_clone, bam, bai]
            }
            .set { ch_filtered2_endo_exo_bam_bai }

        // Create a new channel with just the BAMs    
        ch_filtered2_endo_exo_bam_bai
            .map { meta, bam, bai ->
                [meta, bam]
            }
            .branch { meta, bam ->
                endo: meta.genome == params.genome
                exo: meta.genome == params.spikein_genome
            }
            .set { ch_filtered2_bam }

        // Create a new channel with just the indexes
        ch_filtered2_endo_exo_bam_bai
            .map { meta, bam, bai ->
                [meta, bai]
            }
            .branch { meta, bai ->
                endo: meta.genome == params.genome
                exo: meta.genome == params.spikein_genome
            }
            .set { ch_filtered2_bai }

        ch_filtered_bam = ch_filtered2_bam.endo
        ch_filtered_exo_bam = ch_filtered2_bam.exo
        ch_filtered_index = ch_filtered2_bai.endo
        ch_filtered_exo_index = ch_filtered2_bai.exo
    }
    else {
        // If no spike-in genome add genome to metas
        ch_filtered_bam
            .join(ch_filtered_index, by: 0)
            .map { meta, bam, bai ->
                def meta_clone = meta.clone()
                meta_clone.genome = params.genome
                [meta_clone, bam, bai]
            }
            .set { ch_filtered_bam_bai }

        ch_filtered_bam = ch_filtered_bam_bai.map { meta, bam, bai -> [meta, bam] }
        ch_filtered_index = ch_filtered_bam_bai.map { meta, bam, bai -> [meta, bai] }
    }


    //
    // SUBWORKFLOW: Allocation of multimappers
    //
    if (params.multimap_allocation_method && params.multimap_allocation_method != 'chromap') {
        BAM_ALLOCATE_MULTIMAPPERS_ENDO(
            ch_filtered_bam,
            ch_fasta,
            params.multimap_allocation_method
        )
        ch_filtered_bam = BAM_ALLOCATE_MULTIMAPPERS_ENDO.out.bam
        ch_filtered_index = BAM_ALLOCATE_MULTIMAPPERS_ENDO.out.bai
        ch_samtools_stats_summary = ch_samtools_stats_summary.mix(BAM_ALLOCATE_MULTIMAPPERS_ENDO.out.stats)
        ch_multiqc_files = ch_multiqc_files.mix(BAM_ALLOCATE_MULTIMAPPERS_ENDO.out.stats.collect { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(BAM_ALLOCATE_MULTIMAPPERS_ENDO.out.flagstat.collect { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(BAM_ALLOCATE_MULTIMAPPERS_ENDO.out.idxstats.collect { it[1] })
        ch_versions = ch_versions.mix(BAM_ALLOCATE_MULTIMAPPERS_ENDO.out.versions)

        if (params.allocate_exogenous) {
            BAM_ALLOCATE_MULTIMAPPERS_EXO(
                ch_filtered_exo_bam,
                ch_fasta,
                params.multimap_allocation_method
            )
            ch_filtered_exo_bam = BAM_ALLOCATE_MULTIMAPPERS_EXO.out.bam
            ch_filtered_exo_index = BAM_ALLOCATE_MULTIMAPPERS_EXO.out.bai
            ch_samtools_stats_summary = ch_samtools_stats_summary.mix(BAM_ALLOCATE_MULTIMAPPERS_EXO.out.stats)
            ch_multiqc_files = ch_multiqc_files.mix(BAM_ALLOCATE_MULTIMAPPERS_EXO.out.stats.collect { it[1] })
            ch_multiqc_files = ch_multiqc_files.mix(BAM_ALLOCATE_MULTIMAPPERS_EXO.out.flagstat.collect { it[1] })
            ch_multiqc_files = ch_multiqc_files.mix(BAM_ALLOCATE_MULTIMAPPERS_EXO.out.idxstats.collect { it[1] })
            ch_versions = ch_versions.mix(BAM_ALLOCATE_MULTIMAPPERS_EXO.out.versions)
        }
    }

    // Mix the exogenous and endogenous BAM and index files
    ch_filtered_bam
        .mix(ch_filtered_exo_bam)
        .set { ch_filtered_bam }

    ch_filtered_index
        .mix(ch_filtered_exo_index)
        .set { ch_filtered_index }

    // Split the BAM and indexes into atacseq and other (for shifting)
    ch_filtered_bam
        .branch { meta, bam ->
            atacseq: meta.exp_type == 'atacseq'
            other: meta.exp_type != 'atacseq'
        }
        .set { ch_filtered_bam }

    ch_filtered_index
        .branch { meta, index ->
            atacseq: meta.exp_type == 'atacseq'
            other: meta.exp_type != 'atacseq'
        }
        .set { ch_filtered_index }

    //
    // SUBWORKFLOW: Shift ATAC-seq reads
    //
    BAM_SHIFT_READS(
        ch_filtered_bam.atacseq.join(ch_filtered_index.atacseq, by: 0),
        ch_fasta
    )
    ch_filtered_bam = ch_filtered_bam.other.mix(BAM_SHIFT_READS.out.bam)
    ch_filtered_index = ch_filtered_index.other.mix(BAM_SHIFT_READS.out.bai)
    ch_filtered_bam_bai = ch_filtered_bam.join(ch_filtered_index, by: 0)
    ch_versions = ch_versions.mix(BAM_SHIFT_READS.out.versions)

    if (!params.skip_flT3) {
        //
        // MODULE: Final filtering of BAM file with SAMBAMBA (quality filtering)
        //
        // TODO: fix that the same blacklist is used for both the endogenous and exogenous BAM files
        BAM_FILTER_SAMBAMBA_FLT3(
            ch_filtered_bam.join(ch_filtered_index, by: 0),
            Channel.value([[:], []]),
            ch_fasta
        )
        ch_filtered_bam = BAM_FILTER_SAMBAMBA_FLT3.out.bam
        ch_filtered_index = BAM_FILTER_SAMBAMBA_FLT3.out.bai
        ch_samtools_stats_summary = ch_samtools_stats_summary.mix(BAM_FILTER_SAMBAMBA_FLT3.out.stats)
        ch_multiqc_files = ch_multiqc_files.mix(BAM_FILTER_SAMBAMBA_FLT3.out.stats.collect { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(BAM_FILTER_SAMBAMBA_FLT3.out.flagstat.collect { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(BAM_FILTER_SAMBAMBA_FLT3.out.idxstats.collect { it[1] })
        ch_versions = ch_versions.mix(BAM_FILTER_SAMBAMBA_FLT3.out.versions)

        //
        // MODULE: Extract total mapped reads from flagstats
        //
        BAM_FLAGSTAT_MAPPED_FLT3(
            BAM_FILTER_SAMBAMBA_FLT3.out.flagstat
        )
        ch_versions = ch_versions.mix(BAM_FLAGSTAT_MAPPED_FLT3.out.versions)

        // Extract the total mapped reads from the text file
        BAM_FLAGSTAT_MAPPED_FLT3.out.txt
            .map { meta, total ->
                [meta, total.splitCsv(header: false)[0][0]]
            }
            .set { ch_flT3_total }

        // Add the total_mapped_reads to the bams' and bais' metas
        ch_filtered_bam
            .combine(ch_filtered_index, by: 0)
            .map { meta, bam, bai ->
                [meta, bam, bai]
            }
            .combine(ch_flT3_total, by: 0)
            .map { meta, bam, bai, total ->
                def meta_clone = meta.clone()
                meta_clone.flT3_total_mapped_reads = total.toDouble()
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
    }

    ch_pre_flTbl_bam = Channel.empty()
    ch_pre_flTbl_index = Channel.empty()
    if (!params.skip_flTbl) {

        // These are to later run TE counting on both pre- and post-blacklist-filtering BAM files
        ch_pre_flTbl_bam = ch_filtered_bam
        ch_pre_flTbl_index = ch_filtered_index

        // Separating endogenous and exogenous samples
        // TODO: could add a param "exo_blacklist" to use a different blacklist
        // for exogenous samples instead of skipping them
        ch_filtered_bam_bai
            .branch { meta, bam, bai ->
                endo: meta.genome == params.genome
                exo: meta.genome == params.spikein_genome
            }
            .set { ch_flt_bam_bai_by_genome }

        ch_flt_bam_bai_by_genome
            .exo
            .map { meta, bam, bai ->
                def meta_clone = meta.clone()
                if (meta.flT3_total_mapped_reads) {
                    meta_clone.flTbl_total_mapped_reads = meta.flT3_total_mapped_reads
                } else {
                    meta_clone.flTbl_total_mapped_reads = meta.flT2_total_mapped_reads
                }
                [meta_clone, bam, bai]
            }
            .set { ch_flt_bam_bai_by_genome_exo }
   
        //
        // SUBWORKFLOW: Filter BAM file with SAMBAMBA using blacklist (whitelist)
        //
        BAM_FILTER_BLACKLIST(
            ch_flt_bam_bai_by_genome.endo,
            ch_whitelist,
            ch_fasta
        )
        ch_filtered_bam = BAM_FILTER_BLACKLIST.out.bam.mix(ch_flt_bam_bai_by_genome_exo.map { meta, bam, bai -> [meta, bam] })
        ch_filtered_index = BAM_FILTER_BLACKLIST.out.bai.mix(ch_flt_bam_bai_by_genome_exo.map { meta, bam, bai -> [meta, bai] })
        ch_filtered_bam_bai = ch_filtered_bam.join(ch_filtered_index, by: 0)
        ch_samtools_stats_summary = ch_samtools_stats_summary.mix(BAM_FILTER_BLACKLIST.out.stats)
        ch_multiqc_files = ch_multiqc_files.mix(BAM_FILTER_BLACKLIST.out.multiqc_files)
        ch_versions = ch_versions.mix(BAM_FILTER_BLACKLIST.out.versions)
    }

    //
    // MODULE: Picard post alignment QC
    //
    if (!params.skip_picard_metrics) {
        PICARD_COLLECTMULTIPLEMETRICS(
            ch_filtered_bam_bai,
            ch_fasta,
            ch_fai
        )
        ch_multiqc_files = ch_multiqc_files.mix(PICARD_COLLECTMULTIPLEMETRICS.out.metrics.collect { it[1] })
        ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions.first())
    }

    //
    // MODULE: Phantompeaktools strand cross-correlation and QC metrics
    //
    if (!params.skip_spp) {
        PHANTOMPEAKQUALTOOLS(
            ch_filtered_bam
        )
        ch_multiqc_files = ch_multiqc_files.mix(PHANTOMPEAKQUALTOOLS.out.spp.collect { it[1] })
        ch_versions = ch_versions.mix(PHANTOMPEAKQUALTOOLS.out.versions.first())

        //
        // MODULE: MultiQC custom content for Phantompeaktools
        //
        MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS(
            PHANTOMPEAKQUALTOOLS.out.spp.join(PHANTOMPEAKQUALTOOLS.out.rdata, by: [0]),
            ch_spp_nsc_header,
            ch_spp_rsc_header,
            ch_spp_correlation_header
        )
        ch_multiqc_files = ch_multiqc_files.mix(MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS.out.nsc.collect { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS.out.rsc.collect { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS.out.correlation.collect { it[1] })
        ch_versions = ch_versions.mix(MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS.out.versions.first())
    }

    if (params.bam_downsampling_method) {
        //
        // SUBWORKFLOW: Downsample IP and input control BAM files
        //
        BAM_DOWNSAMPLE(
            ch_filtered_bam_bai,
            ch_fasta,
            ch_fai,
            params.genome,
            params.spikein_genome,
            params.bam_downsampling_method,
            params.downsampling_endo_threshold,
            params.downsampling_exo_threshold,
            params.dSp_use_flT2_total
        )
        ch_filtered_bam = BAM_DOWNSAMPLE.out.bam
        ch_filtered_index = BAM_DOWNSAMPLE.out.bai
        ch_filtered_bam_bai = ch_filtered_bam.join(ch_filtered_index, by: 0)
        ch_samtools_stats_summary = ch_samtools_stats_summary.mix(BAM_DOWNSAMPLE.out.stats)
        ch_multiqc_files = ch_multiqc_files.mix(BAM_DOWNSAMPLE.out.stats.collect { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(BAM_DOWNSAMPLE.out.flagstat.collect { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(BAM_DOWNSAMPLE.out.idxstats.collect { it[1] })
        ch_versions = ch_versions.mix(BAM_DOWNSAMPLE.out.versions.first())
    }

    //
    // SUBWORKFLOW: Normalized bigWig coverage tracks
    //
    BAM_NORMALIZE_BIGWIG_DEEPTOOLS(
        ch_filtered_bam_bai,
        ch_chrom_sizes_endo,
        ch_chrom_sizes_exo,
        params.coverage_bin_size,
        params.genome,
        params.spikein_genome,
        params.skip_srpm,
        params.skip_cisrpm,
        params.skip_cisrpmsoi,
        params.skip_plot_profile,
        params.rpm_use_flT2_total,
        params.srpm_use_flT2_total,
        params.cisrpm_use_flT2_total
    )
    ch_versions = ch_versions.mix(BAM_NORMALIZE_BIGWIG_DEEPTOOLS.out.versions)


    if (!params.skip_plot_profile) {
        //
        // MODULE: deepTools matrix generation for plotting
        //
        DEEPTOOLS_COMPUTEMATRIX(
            BAM_NORMALIZE_BIGWIG_DEEPTOOLS.out.bigwig_binsize1,
            ch_gene_bed.map { it[1] }
        )
        ch_versions = ch_versions.mix(DEEPTOOLS_COMPUTEMATRIX.out.versions.first())

        //
        // MODULE: deepTools profile plots
        //
        DEEPTOOLS_PLOTPROFILE(
            DEEPTOOLS_COMPUTEMATRIX.out.matrix
        )
        ch_multiqc_files = ch_multiqc_files.mix(DEEPTOOLS_PLOTPROFILE.out.table.collect { it[1] })
        ch_versions = ch_versions.mix(DEEPTOOLS_PLOTPROFILE.out.versions.first())

        //
        // MODULE: deepTools heatmaps
        //
        DEEPTOOLS_PLOTHEATMAP(
            DEEPTOOLS_COMPUTEMATRIX.out.matrix
        )
        ch_versions = ch_versions.mix(DEEPTOOLS_PLOTHEATMAP.out.versions.first())
    }

    // Removing the exogenous samples from the filtered_bam_bai channel
    ch_filtered_bam = ch_filtered_bam.filter { it[0].genome == params.genome }
    ch_filtered_index = ch_filtered_index.filter { it[0].genome == params.genome }
    ch_filtered_bam_bai = ch_filtered_bam_bai.filter { it[0].genome == params.genome }

    //
    // SUBWORKFLOW: Counting reads in transposable elements
    //
    if (!params.skip_te_counting) {
        TE_COUNTING(
            // Here we run TE counting on both pre- and post-blacklist-filtering BAM files
            ch_filtered_bam.mix(ch_pre_flTbl_bam.filter { it[0].genome == params.genome }),
            ch_fasta,
            false,
            ch_tecount_gene_index,
            ch_tecount_te_index,
            ch_telocal_gene_index,
            ch_telocal_te_index,
            params.skip_telocal
        )
        ch_versions = ch_versions.mix(TE_COUNTING.out.versions.first())
    }

    //
    // SUBWORKFLOW: Call consensus regions with Consenrich and ROCCO
    //
    if (!params.skip_consenrich) {
        BAM_PEAKS_CALL_QC_ANNOTATE_CONSENRICH_HOMER(
            ch_filtered_bam_bai,
            ch_chrom_sizes_endo,
            ch_blacklist,
            ch_sparsebed.ifEmpty([[:], []]),
            ch_active_regions.ifEmpty([[:], []]),
            ch_rocco_params,
            ch_effective_gsize
        )
        ch_versions = ch_versions.mix(BAM_PEAKS_CALL_QC_ANNOTATE_CONSENRICH_HOMER.out.versions.first())
    }

    //
    // Create channel for downstream processes: [ meta, [ ip_bam, control_bam ] [ ip_bai, control_bai ] ]
    //

    // Branch channels based on if input control is present
    ch_filtered_bam_bai
        .branch { meta, bam, bai ->
            ips_with_control: meta.control
                return [meta.control, meta.antibody, meta, bam, bai]
            ips_wo_control: !meta.control && !meta.is_control
                return [meta.id, meta.antibody, meta, bam, bai]
            controls: !meta.control && meta.is_control
                return [meta.id, meta, bam, bai]
        }
        .set { ch_bam_bai_by_type }

    // For non-downsampled files, duplicate input controls for each antibody 
    ch_bam_bai_by_type
        .controls
        .branch { id, meta, bam, bai ->
            dsp: meta.control_of_antibody && meta.dSp_total_mapped_reads
                return [id, meta.control_of_antibody, meta, bam, bai]
            not_dsp: !meta.control_of_antibody && !meta.dSp_total_mapped_reads
                return [id, meta, bam, bai]
        }
        .set { ch_bam_controls }
    
    ch_bam_controls
        .not_dsp
        .combine(ch_bam_bai_by_type.ips_with_control, by: 0) // combine by control id only
        .map { control_id, control_meta, control_bam, control_bai, ip_antibody, ip_meta, ip_bam, ip_bai ->
            def meta_clone = control_meta.clone()
            meta_clone.control_of_antibody = ip_antibody
            [ control_id, meta_clone.control_of_antibody, meta_clone, control_bam, control_bai ]
        }
        .unique()
        .set { ch_bam_controls_not_dsp }

    ch_bam_bai_by_type
        .ips_with_control
        .combine(ch_bam_controls.dsp.mix(ch_bam_controls_not_dsp), by: [0,1])
        .map { control_id, antibody, ip_meta, ip_bam, ip_bai, control_meta, control_bam, control_bai ->
            [ control_id, antibody, ip_meta, ip_bam, ip_bai, control_bam, control_bai ]
        }
        .mix(ch_bam_bai_by_type.ips_wo_control)
        // ips_wo_control do not have control_bam (it[5]) and control_bai (it[6])
        .map { it -> [ it[2], [ it[3] ] + [ (it[5] ?: []) ], [ it[4] ] + [ (it[6] ?: []) ] ] }
        .set { ch_ip_control_bam_bai }

    // TODO: Print to file for debuggin
    ch_ip_control_bam_bai
        .map { meta, bams, bais ->
            "${meta}\t${bams}\t${bais}"
        }
        .collectFile(name: 'ch_ip_control_bam_bai.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug")

    //
    // MODULE: deepTools plotFingerprint joint QC for IP and control
    //
    if (!params.skip_plot_fingerprint) {
        DEEPTOOLS_PLOTFINGERPRINT(
            ch_ip_control_bam_bai
        )
        ch_multiqc_files = ch_multiqc_files.mix(DEEPTOOLS_PLOTFINGERPRINT.out.matrix.collect { it[1] })
        ch_versions = ch_versions.mix(DEEPTOOLS_PLOTFINGERPRINT.out.versions.first())
    }

    // Create channels: [ meta, ip_bam, control_bam ]
    ch_ip_control_bam_bai
        .map { meta, bams, bais ->
            [meta, bams[0], (bams[1] ?: [])]
        }
        .set { ch_ip_control_bam }


    // separate samples based on meta.exp_type
    ch_ip_control_bam_cs = Channel.empty()
    ch_ip_control_bam_cs = ch_ip_control_bam.filter { !(it[0].exp_type in ['scarseq', 'ChIP-exo', 'OK-seq']) }
    //ch_ip_control_bam_cs = ch_ip_control_bam.filter { it[0].exp_type != 'scarseq' && it[0].exp_type != 'ChIP-exo' }

    // TODO: Print to file for debuggin
    ch_ip_control_bam_cs
        .map { meta, ip_bam, control_bam ->
            "${meta.id}\t${ip_bam}\t${control_bam}"
        }
        .collectFile(name: 'ch_ip_control_bam_cs.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug")

    //
    // SUBWORKFLOW: Call peaks with epic2, annotate with HOMER and perform downstream QC
    //
    ch_epic2_peaks = Channel.empty()
    ch_epic2_frip_multiqc = Channel.empty()
    ch_epic2_peak_count_multiqc = Channel.empty()
    ch_epic2_plot_homer_annotatepeaks_tsv = Channel.empty()
    if (!params.skip_epic2) {
        BAM_PEAKS_CALL_QC_ANNOTATE_EPIC2_HOMER(
            ch_filtered_bam.filter { !(it[0].exp_type in ['scarseq', 'ChIP-exo', 'OK-seq']) },
            ch_fasta,
            ch_gtf,
            ch_chrom_sizes_endo,
            ch_effective_gfraction,
            ".annotatePeaks.txt",
            ch_epic2_peak_count_header,
            ch_epic2_frip_score_header,
            ch_epic2_peak_annotation_header,
            params.skip_peak_annotation,
            params.skip_peak_qc
        )
        ch_epic2_peaks = BAM_PEAKS_CALL_QC_ANNOTATE_EPIC2_HOMER.out.peaks
        ch_epic2_frip_multiqc = BAM_PEAKS_CALL_QC_ANNOTATE_EPIC2_HOMER.out.frip_multiqc
        ch_epic2_peak_count_multiqc = BAM_PEAKS_CALL_QC_ANNOTATE_EPIC2_HOMER.out.peak_count_multiqc
        ch_epic2_plot_homer_annotatepeaks_tsv = BAM_PEAKS_CALL_QC_ANNOTATE_EPIC2_HOMER.out.plot_homer_annotatepeaks_tsv
        ch_versions = ch_versions.mix(BAM_PEAKS_CALL_QC_ANNOTATE_EPIC2_HOMER.out.versions)
    }


    if (!params.skip_edd) {
        //
        // MODULE: Call peaks with EDD
        //
        EDD (
            ch_ip_control_bam_cs,
            ch_chrom_sizes_endo,
            ch_blacklist
        )
        ch_versions = ch_versions.mix(EDD.out.versions.first())
    }

    //
    // SUBWORKFLOW: Call peaks with MACS3, annotate with HOMER and perform downstream QC
    //
    BAM_PEAKS_CALL_QC_ANNOTATE_MACS3_HOMER(
        ch_ip_control_bam_cs,
        ch_fasta,
        ch_gtf,
        ch_chrom_sizes_endo,
        ch_effective_gsize,
        "_peaks.annotatePeaks.txt",
        ch_peak_count_header,
        ch_frip_score_header,
        ch_peak_annotation_header,
        params.narrow_peak,
        params.skip_peak_annotation,
        params.skip_peak_qc,
        params.skip_bdgcmp
    )
    ch_multiqc_files = ch_multiqc_files.mix(BAM_PEAKS_CALL_QC_ANNOTATE_MACS3_HOMER.out.frip_multiqc.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(BAM_PEAKS_CALL_QC_ANNOTATE_MACS3_HOMER.out.peak_count_multiqc.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(BAM_PEAKS_CALL_QC_ANNOTATE_MACS3_HOMER.out.plot_homer_annotatepeaks_tsv.collect { it[1] })
    ch_versions = ch_versions.mix(BAM_PEAKS_CALL_QC_ANNOTATE_MACS3_HOMER.out.versions)

    //
    //  Consensus peaks analysis
    //
    ch_macs3_consensus_bed_lib = Channel.empty()
    ch_macs3_consensus_txt_lib = Channel.empty()
    if (!params.skip_consensus_peaks) {
        // Create channels: [ antibody, [ ip_bams ] ]
        ch_ip_control_bam_cs
            .map { meta, ip_bam, control_bam ->
                [meta.antibody, ip_bam]
            }
            .groupTuple()
            .set { ch_antibody_bams }

        BED_CONSENSUS_QUANTIFY_QC_BEDTOOLS_FEATURECOUNTS_DESEQ2(
            BAM_PEAKS_CALL_QC_ANNOTATE_MACS3_HOMER.out.peaks,
            ch_antibody_bams,
            ch_fasta.map { it[1] },
            ch_gtf.map { it[1] },
            ch_deseq2_pca_header,
            ch_deseq2_clustering_header,
            params.narrow_peak,
            params.skip_peak_annotation,
            params.skip_deseq2_qc
        )
        ch_macs3_consensus_bed_lib = BED_CONSENSUS_QUANTIFY_QC_BEDTOOLS_FEATURECOUNTS_DESEQ2.out.consensus_bed
        ch_macs3_consensus_txt_lib = BED_CONSENSUS_QUANTIFY_QC_BEDTOOLS_FEATURECOUNTS_DESEQ2.out.consensus_txt
        ch_multiqc_files = ch_multiqc_files.mix(BED_CONSENSUS_QUANTIFY_QC_BEDTOOLS_FEATURECOUNTS_DESEQ2.out.featurecounts_summary.collect { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(BED_CONSENSUS_QUANTIFY_QC_BEDTOOLS_FEATURECOUNTS_DESEQ2.out.deseq2_qc_pca_multiqc.collect { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(BED_CONSENSUS_QUANTIFY_QC_BEDTOOLS_FEATURECOUNTS_DESEQ2.out.deseq2_qc_dists_multiqc.collect { it[1] })
        ch_versions = ch_versions.mix(BED_CONSENSUS_QUANTIFY_QC_BEDTOOLS_FEATURECOUNTS_DESEQ2.out.versions)
    }

    //
    // SUBWORKFLOW: Call peaks with Genrich, annotate with HOMER and perform downstream QC
    //
    ch_genrich_peaks = Channel.empty()
    if (!params.skip_genrich) {
        BAM_PEAKS_CALL_QC_ANNOTATE_GENRICH_HOMER(
            ch_filtered_bam.filter { !(it[0].exp_type in ['scarseq', 'ChIP-exo', 'OK-seq']) },
            ch_fasta,
            ch_gtf,
            ch_blacklist,
            ".annotatePeaks.txt",
            ch_gr_peak_count_header,
            ch_gr_frip_score_header,
            ch_gr_peak_annotation_header,
            params.narrow_peak,
            params.skip_peak_annotation,
            params.skip_peak_qc
        )
        ch_genrich_peaks = BAM_PEAKS_CALL_QC_ANNOTATE_GENRICH_HOMER.out.peaks
        ch_multiqc_files = ch_multiqc_files.mix(BAM_PEAKS_CALL_QC_ANNOTATE_GENRICH_HOMER.out.frip_multiqc.collect { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(BAM_PEAKS_CALL_QC_ANNOTATE_GENRICH_HOMER.out.peak_count_multiqc.collect { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(BAM_PEAKS_CALL_QC_ANNOTATE_GENRICH_HOMER.out.plot_homer_annotatepeaks_tsv.collect { it[1] })
        ch_versions = ch_versions.mix(BAM_PEAKS_CALL_QC_ANNOTATE_GENRICH_HOMER.out.versions)
    }

    //
    // SUBWORKFLOW: Call peaks with MACE (for ChIP-exo samples)
    //
    if (!params.skip_mace) {
        BAM_PEAKS_CALL_QC_ANNOTATE_MACE_HOMER(
            ch_filtered_bam_bai.filter { it[0].exp_type == 'ChIP-exo' },
            ch_fasta,
            ch_gtf,
            ch_chrom_sizes_endo,
            ".annotatePeaks.txt",
            ch_mace_peak_count_header,
            ch_mace_frip_score_header,
            ch_mace_peak_annotation_header,
            params.skip_peak_annotation,
            params.skip_peak_qc
        )
        ch_multiqc_files = ch_multiqc_files.mix(BAM_PEAKS_CALL_QC_ANNOTATE_MACE_HOMER.out.frip_multiqc.collect { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(BAM_PEAKS_CALL_QC_ANNOTATE_MACE_HOMER.out.peak_count_multiqc.collect { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(BAM_PEAKS_CALL_QC_ANNOTATE_MACE_HOMER.out.plot_homer_annotatepeaks_tsv.collect { it[1] })
        ch_versions = ch_versions.mix(BAM_PEAKS_CALL_QC_ANNOTATE_MACE_HOMER.out.versions)
    }


    ch_filtered_bam_ss = Channel.empty()
    ch_filtered_bam_ss = ch_filtered_bam.filter { it[0].exp_type in ['scarseq', 'OK-seq'] }

    // TODO: remove when optional inputs to subworkflows are implemented
    // Make ch_chrom_sizes_endo empty if there are no scarseq samples
    // This is to avoid unnecessarily running modules in the BAM_CREATE_PARTITIONS
    ch_chrom_sizes_endo
        .combine(ch_filtered_bam_ss)
        .first()
        .map { sizes_meta, sizes, ss_meta, ss_bam ->
            [sizes_meta, sizes]
        }
        .set { ch_chrom_sizes_endo_ss }

    //
    // SUBWORKFLOW: SCAR-seq and OK-seq analysis: partitioning of reads
    //
    ch_partition_smooth = Channel.empty()
    BAM_CREATE_PARTITIONS(
        ch_filtered_bam_ss,
        ch_chrom_sizes_endo_ss,
        ch_blacklist,
        ch_initiation_zones,
        params.rpm_use_flT2_total,
        params.smooth_radius,
        params.derivative_radius,
        params.zero_crossing_radius
    )
    ch_partition_smooth = BAM_CREATE_PARTITIONS.out.tab
    ch_versions = ch_versions.mix(BAM_CREATE_PARTITIONS.out.versions)

    //
    // SUBWORKFLOW: Create SAMtools summary table
    //
    SAMTOOLS_STATS_SUMMARY(
        ch_samtools_stats_summary,
        params.genome,
        params.spikein_genome ?: Channel.value([])
    )
    ch_versions = ch_versions.mix(SAMTOOLS_STATS_SUMMARY.out.versions)

    //
    // MODULE: Create IGV session
    //
    if (!params.skip_igv) {
        IGV(
            params.aligner,
            params.multimap_allocation_method ? params.multimap_allocation_method == 'chromap' ? 'chromap_allocation' : params.multimap_allocation_method + '/' : '',
            params.narrow_peak ? 'narrow_peak' : 'broad_peak',
            ch_fasta.map { it[1] },
            BAM_NORMALIZE_BIGWIG_DEEPTOOLS.out.bigwig_endo.collect { it[1] }.ifEmpty([]),
            ch_epic2_peaks.collect { it[1] }.ifEmpty([]),
            ch_genrich_peaks.collect { it[1] }.ifEmpty([]),
            BAM_PEAKS_CALL_QC_ANNOTATE_MACS3_HOMER.out.peaks.collect { it[1] }.ifEmpty([]),
            ch_macs3_consensus_bed_lib.collect { it[1] }.ifEmpty([]),
            ch_macs3_consensus_txt_lib.collect { it[1] }.ifEmpty([])
        )
        ch_versions = ch_versions.mix(IGV.out.versions)
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'glseq_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    if (!params.skip_multiqc) {

        // Load MultiQC configuration files
        ch_multiqc_config = Channel.fromPath("${projectDir}/assets/multiqc_config.yml", checkIfExists: true)
        ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()
        ch_multiqc_logo = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo) : Channel.empty()

        // Prepare the workflow summary
        ch_workflow_summary = Channel.value(
                paramsSummaryMultiqc(
                    paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
                )
            )
            .collectFile(name: 'workflow_summary_mqc.yaml')

        // Prepare the methods section
        // ch_methods_description = Channel.value(
        //     methodsDescriptionText(
        //         params.multiqc_methods_description
        //             ? file(params.multiqc_methods_description)
        //             : file("$projectDir/workflows/assets/multiqc/methods_description_template.yml", checkIfExists: true)
        //     )
        // ).collectFile(name: 'methods_description_mqc.yaml')

        // Add summary, versions, and methods to the MultiQC input file list
        ch_multiqc_files = ch_multiqc_files
            .mix(ch_workflow_summary)
            .mix(ch_collated_versions)
        // .mix(ch_methods_description)


        // Provide MultiQC with rename patterns to ensure it uses sample names
        // for single-techrep samples not processed by CAT_FASTQ, and trims out
        // _raw or _trimmed

        // ch_name_replacements = ch_fastq
        //     .map{ meta, reads ->
        //         def name1 = file(reads[0][0]).simpleName + "\t" + meta.id + '_1'
        //         def fastqcnames = meta.id + "_raw\t" + meta.id + "\n" + meta.id + "_trimmed\t" + meta.id
        //         if (reads[0][1] ){
        //             def name2 = file(reads[0][1]).simpleName + "\t" + meta.id + '_2'
        //             def fastqcnames1 = meta.id + "_raw_1\t" + meta.id + "_1\n" + meta.id + "_trimmed_1\t" + meta.id + "_1"
        //             def fastqcnames2 = meta.id + "_raw_2\t" + meta.id + "_2\n" + meta.id + "_trimmed_2\t" + meta.id + "_2"
        //             return [ name1, name2, fastqcnames1, fastqcnames2 ]
        //         } else{
        //             return [ name1, fastqcnames ]
        //         }
        //     }
        //     .flatten()
        //     .collectFile(name: 'name_replacement.txt', newLine: true)

        MULTIQC(
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList(),
            [],
            []
        )
        ch_multiqc_report = MULTIQC.out.report
    }

    emit:
    multiqc_report = ch_multiqc_report // channel: /path/to/multiqc_report.html
    versions       = ch_versions // channel: [ path(versions.yml) ]
}