/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { IGV                                 } from '../modules/local/igv/main'
include { MULTIQC                             } from '../modules/local/multiqc/main'
include { MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS } from '../modules/local/multiqc_custom_phantompeakqualtools/main'
include { 
    BAM_FLAGSTAT_MAPPED as BAM_FLAGSTAT_MAPPED_FLT1
    BAM_FLAGSTAT_MAPPED as BAM_FLAGSTAT_MAPPED_FLT2
    BAM_FLAGSTAT_MAPPED as BAM_FLAGSTAT_MAPPED_FLT3
                                              } from '../modules/local/bam_flagstat_mapped/main'
include { CONSENRICH             } from '../modules/local/consenrich/main'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_grothlab_glseq_pipeline'
include { INPUT_CHECK            } from '../subworkflows/local/input_check/main'
include {
    BAM_FILTER_SAMBAMBA as BAM_FILTER_SAMBAMBA_FLT1
    BAM_FILTER_SAMBAMBA as BAM_FILTER_SAMBAMBA_FLT3
                                } from '../subworkflows/local/bam_filter_sambamba/main'
include { BAM_SPIKEIN_SPLIT   } from '../subworkflows/local/bam_spikein_split/main'
include { FASTQ_FASTQC_UMITOOLS_UMITRANSFER_TRIMGALORE      } from '../subworkflows/local/fastq_fastqc_umitools_umitransfer_trimgalore/main'
// include { BAM_BEDGRAPH_BIGWIG_BEDTOOLS_UCSC                       } from '../subworkflows/local/bam_bedgraph_bigwig_bedtools_ucsc/main'
include { BAM_PEAKS_CALL_QC_ANNOTATE_MACS3_HOMER                  } from '../subworkflows/local/bam_peaks_call_qc_annotate_macs3_homer/main'
include { BAM_PEAKS_CALL_QC_ANNOTATE_GENRICH_HOMER                  } from '../subworkflows/local/bam_peaks_call_qc_annotate_genrich_homer/main'
include { BED_CONSENSUS_QUANTIFY_QC_BEDTOOLS_FEATURECOUNTS_DESEQ2 } from '../subworkflows/local/bed_consensus_quantify_qc_bedtools_featurecounts_deseq2/main'
include { BAM_CREATE_SCAR_PARTITIONS } from '../subworkflows/local/bam_create_scar_partitions/main'
include { BAM_ALLOCATE_MULTIMAPPERS as BAM_ALLOCATE_MULTIMAPPERS_ENDO } from '../subworkflows/local/bam_allocate_multimappers/main'
include { BAM_ALLOCATE_MULTIMAPPERS as BAM_ALLOCATE_MULTIMAPPERS_EXO } from '../subworkflows/local/bam_allocate_multimappers/main'
include { BAM_SHIFT_READS            } from '../subworkflows/local/bam_shift_reads/main'
include { SAMTOOLS_STATS_SUMMARY                    } from '../subworkflows/local/samtools_stats_summary/main'
// include { COUNT_READS_IN_BINS                        } from '../subworkflows/local/count_reads_in_bins/main'
// include { BAM_CHORSEQ_RRPM                           } from '../subworkflows/local/bam_chorseq_rrpm/main'
include { BAM_NORMALIZE_BIGWIG_DEEPTOOLS           } from '../subworkflows/local/bam_normalize_bigwig_deeptools/main'
include { BAM_DOWNSAMPLE                            } from '../subworkflows/local/bam_downsample/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { SAMTOOLS_INDEX                } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_COLLATE              } from '../modules/nf-core/samtools/collate/main'
include { PICARD_MERGESAMFILES          } from '../modules/nf-core/picard/mergesamfiles/main'
include { PICARD_COLLECTMULTIPLEMETRICS } from '../modules/nf-core/picard/collectmultiplemetrics/main'
include { PRESEQ_LCEXTRAP               } from '../modules/nf-core/preseq/lcextrap/main'
include { PHANTOMPEAKQUALTOOLS          } from '../modules/nf-core/phantompeakqualtools/main'
include { DEEPTOOLS_COMPUTEMATRIX       } from '../modules/nf-core/deeptools/computematrix/main'
include { DEEPTOOLS_PLOTPROFILE         } from '../modules/nf-core/deeptools/plotprofile/main'
include { DEEPTOOLS_PLOTHEATMAP         } from '../modules/nf-core/deeptools/plotheatmap/main'
include { DEEPTOOLS_PLOTFINGERPRINT     } from '../modules/nf-core/deeptools/plotfingerprint/main'
include { KHMER_UNIQUEKMERS             } from '../modules/nf-core/khmer/uniquekmers/main'
//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//

// include { FASTQ_FASTQC_UMITOOLS_TRIMGALORE      } from '../subworkflows/nf-core/fastq_fastqc_umitools_trimgalore'
include { FASTQ_ALIGN_BWA                   } from '../subworkflows/nf-core/fastq_align_bwa'
include { FASTQ_ALIGN_BOWTIE2               } from '../subworkflows/nf-core/fastq_align_bowtie2'
include { FASTQ_ALIGN_CHROMAP               } from '../subworkflows/nf-core/fastq_align_chromap'
include { FASTQ_ALIGN_STAR                  } from '../subworkflows/nf-core/fastq_align_star'
include { FASTQ_ALIGN_HISAT2                } from '../subworkflows/nf-core/fastq_align_hisat2'                                                                                                                                                                            
include { BAM_MARKDUPLICATES_PICARD         } from '../subworkflows/nf-core/bam_markduplicates_picard'
include { BAM_DEDUP_UMI                     } from '../subworkflows/nf-core/bam_dedup_umi'
include { BAM_STATS_SAMTOOLS                } from '../subworkflows/nf-core/bam_stats_samtools'
include { BAM_SORT_STATS_SAMTOOLS           } from '../subworkflows/nf-core/bam_sort_stats_samtools'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Header files for MultiQC
ch_spp_nsc_header           = file("$projectDir/assets/multiqc/spp_nsc_header.txt", checkIfExists: true)
ch_spp_rsc_header           = file("$projectDir/assets/multiqc/spp_rsc_header.txt", checkIfExists: true)
ch_spp_correlation_header   = file("$projectDir/assets/multiqc/spp_correlation_header.txt", checkIfExists: true)
ch_peak_count_header        = file("$projectDir/assets/multiqc/peak_count_header.txt", checkIfExists: true)
ch_gr_peak_count_header     = file("$projectDir/assets/multiqc/gr_peak_count_header.txt", checkIfExists: true)
ch_frip_score_header        = file("$projectDir/assets/multiqc/frip_score_header.txt", checkIfExists: true)
ch_gr_frip_score_header     = file("$projectDir/assets/multiqc/gr_frip_score_header.txt", checkIfExists: true)
ch_peak_annotation_header   = file("$projectDir/assets/multiqc/peak_annotation_header.txt", checkIfExists: true)
ch_gr_peak_annotation_header = file("$projectDir/assets/multiqc/gr_peak_annotation_header.txt", checkIfExists: true)
ch_deseq2_pca_header        = Channel.value(file("$projectDir/assets/multiqc/deseq2_pca_header.txt", checkIfExists: true))
ch_deseq2_clustering_header = Channel.value(file("$projectDir/assets/multiqc/deseq2_clustering_header.txt", checkIfExists: true))

// Save AWS IGenomes file containing annotation version
def anno_readme = params.genomes[ params.genome ]?.readme
if (anno_readme && file(anno_readme).exists()) {
    file("${params.outdir}/genome/").mkdirs()
    file(anno_readme).copyTo("${params.outdir}/genome/")
}

// // Info required for completion email and summary
// def multiqc_report = []

workflow GLSEQ {

    take:
    ch_input         // channel: path(sample_sheet.csv)
    ch_versions      // channel: [ path(versions.yml) ]
    ch_fasta         // channel: path(genome.fa)
    ch_fai           // channel: path(genome.fai)
    ch_gtf           // channel: path(genome.gtf)
    ch_gene_bed      // channel: path(gene.beds)
    ch_chrom_sizes   // channel: path(chrom.sizes)
    ch_chrom_sizes_endo // path(chrom.sizes.endo)
    ch_chrom_sizes_exo
    ch_scaffolds     // channel: val(scaffolds)
    ch_filtered_bed  // channel: path(filtered.bed)
    ch_blacklist     // channel: path(blacklist.bed)
    ch_sparsebed     // channel: path(sparse.bed)
    ch_initiation_zones // channel: path(initiation_zones)
    ch_bwa_index     // channel: path(bwa/index/)
    ch_bowtie2_index // channel: path(bowtie2/index)
    ch_chromap_index // channel: path(chromap.index)
    ch_star_index    // channel: path(star/index/)
    ch_hisat2_index  // channel: path(hisat2/index)
    ch_splicesites   // channel: path(splicesites)

    main:
    ch_multiqc_files = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input,
        params.seq_center
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    // TODO: OPTIONAL, you can use nf-validation plugin to create an input channel from the samplesheet with Channel.fromSamplesheet("input")
    // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
    // ! There is currently no tooling to help you write a sample sheet schema

    //
    // SUBWORKFLOW: Read QC and trim adapters
    //
    FASTQ_FASTQC_UMITOOLS_UMITRANSFER_TRIMGALORE (
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
    ch_versions = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_UMITRANSFER_TRIMGALORE.out.versions)

    //
    // SUBWORKFLOW: Alignment with BWA & BAM QC
    //
    ch_genome_bam        = Channel.empty()
    ch_genome_bam_index  = Channel.empty()
    ch_samtools_stats    = Channel.empty()
    ch_samtools_flagstat = Channel.empty()
    ch_samtools_idxstats = Channel.empty()
    if (params.aligner == 'bwa') {
        FASTQ_ALIGN_BWA (
            FASTQ_FASTQC_UMITOOLS_UMITRANSFER_TRIMGALORE.out.reads,
            ch_bwa_index,
            params.sort_bam,
            ch_fasta.map{ it[1] }.first()

        )
        ch_genome_bam        = FASTQ_ALIGN_BWA.out.bam
        ch_genome_bam_index  = FASTQ_ALIGN_BWA.out.bai
        ch_samtools_stats    = FASTQ_ALIGN_BWA.out.stats
        ch_samtools_flagstat = FASTQ_ALIGN_BWA.out.flagstat
        ch_samtools_idxstats = FASTQ_ALIGN_BWA.out.idxstats
        ch_versions = ch_versions.mix(FASTQ_ALIGN_BWA.out.versions.first())
    }

    //
    // SUBWORKFLOW: Alignment with Bowtie2 & BAM QC
    //
    // TODO: using first() to convert the tuple to a value channel and make it consumable
    if (params.aligner == 'bowtie2') {
        FASTQ_ALIGN_BOWTIE2 (
            FASTQ_FASTQC_UMITOOLS_UMITRANSFER_TRIMGALORE.out.reads,
            ch_bowtie2_index.first(),
            params.save_unaligned,
            params.sort_bam,
            ch_fasta.first()
        )
        ch_genome_bam        = FASTQ_ALIGN_BOWTIE2.out.bam
        ch_genome_bam_index  = FASTQ_ALIGN_BOWTIE2.out.bai
        ch_samtools_stats    = FASTQ_ALIGN_BOWTIE2.out.stats
        ch_samtools_flagstat = FASTQ_ALIGN_BOWTIE2.out.flagstat
        ch_samtools_idxstats = FASTQ_ALIGN_BOWTIE2.out.idxstats
        ch_versions = ch_versions.mix(FASTQ_ALIGN_BOWTIE2.out.versions.first())
    }

    //
    // SUBWORKFLOW: Alignment with Chromap & BAM QC
    //
    if (params.aligner == 'chromap') {
        FASTQ_ALIGN_CHROMAP (
            FASTQ_FASTQC_UMITOOLS_UMITRANSFER_TRIMGALORE.out.reads,
            ch_chromap_index.first(),
            ch_fasta.first(),
            [],
            [],
            [],
            []
        )

        ch_genome_bam        = FASTQ_ALIGN_CHROMAP.out.bam
        ch_genome_bam_index  = FASTQ_ALIGN_CHROMAP.out.bai
        ch_samtools_stats    = FASTQ_ALIGN_CHROMAP.out.stats
        ch_samtools_flagstat = FASTQ_ALIGN_CHROMAP.out.flagstat
        ch_samtools_idxstats = FASTQ_ALIGN_CHROMAP.out.idxstats
        ch_versions = ch_versions.mix(FASTQ_ALIGN_CHROMAP.out.versions.first())
    }

    //
    // SUBWORKFLOW: Alignment with STAR & BAM QC
    //
    if (params.aligner == 'star') {
        FASTQ_ALIGN_STAR (
            FASTQ_FASTQC_UMITOOLS_UMITRANSFER_TRIMGALORE.out.reads,
            ch_star_index.first(),
            ch_gtf.first(),
            true,
            params.seq_platform ?: '',
            params.seq_center ?: '',
            ch_fasta.first(),
            Channel.of([[:], []])

        )
        ch_genome_bam        = FASTQ_ALIGN_STAR.out.bam
        ch_genome_bam_index  = FASTQ_ALIGN_STAR.out.bai
        ch_transcriptome_bam = FASTQ_ALIGN_STAR.out.bam_transcript
        ch_samtools_stats    = FASTQ_ALIGN_STAR.out.stats
        ch_samtools_flagstat = FASTQ_ALIGN_STAR.out.flagstat
        ch_samtools_idxstats = FASTQ_ALIGN_STAR.out.idxstats
        ch_star_multiqc      = FASTQ_ALIGN_STAR.out.log_final

        ch_versions = ch_versions.mix(FASTQ_ALIGN_STAR.out.versions)
    }

    if (params.aligner == 'hisat2') {
        FASTQ_ALIGN_HISAT2 (
            FASTQ_FASTQC_UMITOOLS_UMITRANSFER_TRIMGALORE.out.reads,
            ch_hisat2_index.first(),
            ch_splicesites.first(),
            ch_fasta.first()
        )
        ch_genome_bam        = FASTQ_ALIGN_HISAT2.out.bam
        ch_genome_bam_index  = FASTQ_ALIGN_HISAT2.out.bai
        ch_samtools_stats    = FASTQ_ALIGN_HISAT2.out.stats
        ch_samtools_flagstat = FASTQ_ALIGN_HISAT2.out.flagstat
        ch_samtools_idxstats = FASTQ_ALIGN_HISAT2.out.idxstats
        ch_versions = ch_versions.mix(FASTQ_ALIGN_HISAT2.out.versions)
    }

    //
    // MODULE: Merge resequenced BAM files
    //
    ch_genome_bam
        .map {
            meta, bam ->
                def meta_clone = meta.clone()
                meta_clone.remove('read_group')
                meta_clone.id = meta_clone.id.split('_')[0..-2].join('_')
                [ meta_clone, bam ]
        }
        .groupTuple(by: 0)
        .map {
            it ->
                [ it[0], it[1].flatten() ]
        }
        .set { ch_sort_bam }

    PICARD_MERGESAMFILES (
        ch_sort_bam,
    )
    ch_merged_bam = PICARD_MERGESAMFILES.out.bam
    ch_versions = ch_versions.mix(PICARD_MERGESAMFILES.out.versions.first().ifEmpty(null))

    SAMTOOLS_INDEX (
        ch_merged_bam
    )
    ch_merged_bam_bai = ch_merged_bam.join(SAMTOOLS_INDEX.out.bai, by: [0])
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    BAM_STATS_SAMTOOLS (
        ch_merged_bam_bai,
        ch_fasta.first()
    )
    ch_merged_bam_stats = BAM_STATS_SAMTOOLS.out.stats
    ch_merged_bam_flagstat = BAM_STATS_SAMTOOLS.out.flagstat
    ch_merged_bam_idxstats = BAM_STATS_SAMTOOLS.out.idxstats
    ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)

    //
    // MODULE: Preseq coverage analysis
    //
    // TODO: this is done on the bams with spike-in included
    ch_preseq_multiqc = Channel.empty()
    if (!params.skip_preseq) {
        PRESEQ_LCEXTRAP (
            ch_merged_bam
        )
        ch_preseq_multiqc = PRESEQ_LCEXTRAP.out.lc_extrap
        ch_versions = ch_versions.mix(PRESEQ_LCEXTRAP.out.versions.first())
    }
    
    ch_dedup_umi_stats = Channel.empty()
    ch_dedup_umi_flagstat = Channel.empty()
    ch_dedup_umi_idxstats = Channel.empty()
    ch_dedup_umi_deduplog = Channel.empty()
    ch_mk_stats = Channel.empty()
    ch_mk_flagstat = Channel.empty()
    ch_mk_idxstats = Channel.empty()
    ch_mk_metrics = Channel.empty()
    // TODO: change this so UMI dedup is evaluated per sample and not for the whole pipeline run
    if (params.with_umi) {

        //
        // SUBWORKFLOW: Deduplicate BAM files
        //
        ch_transcriptome_bam = Channel.empty()
        ch_transcriptome_fasta = Channel.empty()
        BAM_DEDUP_UMI (
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
        ch_dedup_umi_deduplog = BAM_DEDUP_UMI.out.dedup_log

        // TODO: print ch_dedup_umi_deduplog for debugging
        ch_dedup_umi_deduplog
            .map {
                meta, deduplog ->
                   "${meta}\t${deduplog}"
            }
            .collectFile( name: 'ch_dedup_umi_deduplog.txt', newLine: true, sort: false, storeDir: "${params.outdir}") 
        
        ch_versions = ch_versions.mix(BAM_DEDUP_UMI.out.versions.first())

    } else {
        //
        // SUBWORKFLOW: Mark duplicates & filter BAM files
        //
        // TODO: using first() to convert the tuple to a value channel and make it consumable
        BAM_MARKDUPLICATES_PICARD (
            ch_merged_bam,
            ch_fasta.first(),
            ch_fai.first()
        )
        ch_dedup_bam = BAM_MARKDUPLICATES_PICARD.out.bam
        ch_dedup_index = BAM_MARKDUPLICATES_PICARD.out.bai
        ch_mk_stats = BAM_MARKDUPLICATES_PICARD.out.stats
        ch_mk_flagstat = BAM_MARKDUPLICATES_PICARD.out.flagstat
        ch_mk_idxstats = BAM_MARKDUPLICATES_PICARD.out.idxstats
        ch_mk_metrics = BAM_MARKDUPLICATES_PICARD.out.metrics

        ch_versions = ch_versions.mix(BAM_MARKDUPLICATES_PICARD.out.versions)

        //
        // MODULE: Preseq coverage analysis
        //
        // TODO: this is done on the bams with spike-in included
        ch_preseq_multiqc = Channel.empty()
        if (!params.skip_preseq) {
            PRESEQ_LCEXTRAP (
                ch_dedup_bam
            )
            ch_preseq_multiqc = PRESEQ_LCEXTRAP.out.lc_extrap
            ch_versions = ch_versions.mix(PRESEQ_LCEXTRAP.out.versions.first())
        }
    }


    //
    // SUBWORKFLOW: Filter BAM file with SAMBAMBA
    //
    BAM_FILTER_SAMBAMBA_FLT1 (
        ch_dedup_bam.join(ch_dedup_index, by: 0),
        ch_filtered_bed.first(),
        ch_fasta.first()
    )
    ch_filtered_bam = BAM_FILTER_SAMBAMBA_FLT1.out.bam
    ch_filtered_index = BAM_FILTER_SAMBAMBA_FLT1.out.bai
    ch_filtered_stats = BAM_FILTER_SAMBAMBA_FLT1.out.stats
    ch_filtered_flagstat = BAM_FILTER_SAMBAMBA_FLT1.out.flagstat
    ch_filtered_idxstats = BAM_FILTER_SAMBAMBA_FLT1.out.idxstats
    ch_versions = ch_versions.mix(BAM_FILTER_SAMBAMBA_FLT1.out.versions)

    //
    // MODULE: Extract total mapped reads from flagstats
    //
    BAM_FLAGSTAT_MAPPED_FLT1 (
        ch_filtered_flagstat
    )
    ch_versions = ch_versions.mix(BAM_FLAGSTAT_MAPPED_FLT1.out.versions)

    // Extract the total mapped reads from the text file
    BAM_FLAGSTAT_MAPPED_FLT1.out.txt
        .map {
            meta, total ->
                [ meta, total.splitCsv(header:false)[0][0] ]
        }
        .set { ch_flT1_total }

    // Add the total_mapped_reads to the bams' metas
    ch_filtered_bam
        .map {
            meta, bam ->
                [ meta, bam ]
        }
        .combine(ch_flT1_total, by: 0)
        .map {
            meta, bam, total ->
                meta_clone = meta.clone()
                meta_clone.flT1_total_mapped_reads = total.toDouble()
                [ meta_clone, bam ]
        }
        .set { ch_filtered_bam }


    //
    // SUBWORKFLOW: Spike-in splitting
    //
    // TODO: if fasta and gtf are specified but not genome, val keep_genome_string in
    // BAM_SPLIT_BY_GENOME will fail
    ch_filtered2_flagstat = Channel.empty()
    ch_filtered2_stats = Channel.empty()
    ch_filtered2_idxstats = Channel.empty()
    ch_filtered_exo_bam = Channel.empty()
    ch_filtered_exo_index = Channel.empty()
    ch_filtered2_exo_flagstat = Channel.empty()
    ch_filtered2_exo_stats = Channel.empty()
    ch_filtered2_exo_idxstats = Channel.empty()
    ch_filtered2_endo_exo_bam = Channel.empty()
    if (params.spikein_genome) {
        BAM_SPIKEIN_SPLIT (
            ch_filtered_bam,
            ch_fasta.first(),
            ch_filtered_bed.first(),
            params.genome,
            params.spikein_genome
        )

        ch_filtered_bam = BAM_SPIKEIN_SPLIT.out.bam
        ch_filtered_index = BAM_SPIKEIN_SPLIT.out.bai
        ch_filtered2_stats = BAM_SPIKEIN_SPLIT.out.stats
        ch_filtered2_flagstat = BAM_SPIKEIN_SPLIT.out.flagstat
        ch_filtered2_idxstats = BAM_SPIKEIN_SPLIT.out.idxstats

        ch_filtered_exo_bam = BAM_SPIKEIN_SPLIT.out.exo_bam
        ch_filtered_exo_index = BAM_SPIKEIN_SPLIT.out.exo_bai
        ch_filtered2_exo_stats = BAM_SPIKEIN_SPLIT.out.exo_stats
        ch_filtered2_exo_flagstat = BAM_SPIKEIN_SPLIT.out.exo_flagstat
        ch_filtered2_exo_idxstats = BAM_SPIKEIN_SPLIT.out.exo_idxstats
        ch_versions = ch_versions.mix(BAM_SPIKEIN_SPLIT.out.versions.first())
    
        //
        // MODULE: Extract total mapped reads from flagstats
        //
        BAM_FLAGSTAT_MAPPED_FLT2 (
            ch_filtered2_flagstat.mix(ch_filtered2_exo_flagstat)
        )
        ch_versions = ch_versions.mix(BAM_FLAGSTAT_MAPPED_FLT2.out.versions)

        // Extract the total mapped reads from the text file
        BAM_FLAGSTAT_MAPPED_FLT2.out.txt
            .map {
                meta, total ->
                    [ meta, total.splitCsv(header:false)[0][0] ]
            }
            .set { ch_flT2_total }

        // Add the total_mapped_reads both endo and exo bams' and bais' metas
        ch_filtered_bam
            .mix(ch_filtered_exo_bam)
            .join(ch_filtered_index.mix(ch_filtered_exo_index), by: 0)
            .combine(ch_flT2_total, by: 0)
            .map {
                meta, bam, bai, total ->
                    meta_clone = meta.clone()
                    meta_clone.flT2_total_mapped_reads = total.toDouble()
                    [ meta_clone, bam, bai ]
            }
            .set { ch_filtered2_endo_exo_bam_bai }

        // Create a new channel with just the BAMs    
        ch_filtered2_endo_exo_bam_bai
            .map { meta, bam, bai ->
                [ meta, bam ]
            }
            .branch { meta, bam ->
                endo: meta.genome == params.genome
                exo: meta.genome == params.spikein_genome
            }
            .set { ch_filtered2_bam }
        
        // Create a new channel with just the indexes
        ch_filtered2_endo_exo_bam_bai
            .map { meta, bam, bai ->
                [ meta, bai ]
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

    //
    // SUBWORKFLOW: Allocation of multimappers
    //
    ch_allocated_flagstat = Channel.empty()
    ch_allocated_stats = Channel.empty()
    ch_allocated_idxstats = Channel.empty()
    if (params.allocate_n_multimappers && params.allocation_method != 'chromap') {

        BAM_ALLOCATE_MULTIMAPPERS_ENDO (
            ch_filtered_bam,
            ch_fasta,
            params.allocation_method
        )
        ch_filtered_bam = BAM_ALLOCATE_MULTIMAPPERS_ENDO.out.bam
        ch_filtered_index = BAM_ALLOCATE_MULTIMAPPERS_ENDO.out.bai
        ch_allocated_flagstat = BAM_ALLOCATE_MULTIMAPPERS_ENDO.out.flagstat
        ch_allocated_stats = BAM_ALLOCATE_MULTIMAPPERS_ENDO.out.stats
        ch_allocated_idxstats = BAM_ALLOCATE_MULTIMAPPERS_ENDO.out.idxstats
        ch_versions = ch_versions.mix(BAM_ALLOCATE_MULTIMAPPERS_ENDO.out.versions)
    
        ch_exo_allocated_flagstat = Channel.empty()
        ch_exo_allocated_stats = Channel.empty()
        ch_exo_allocated_idxstats = Channel.empty()
        if (params.allocate_exogenous) {
            BAM_ALLOCATE_MULTIMAPPERS_EXO (
                ch_filtered_exo_bam,
                ch_fasta,
                params.allocation_method
            )
            ch_filtered_exo_bam = BAM_ALLOCATE_MULTIMAPPERS_EXO.out.bam
            ch_filtered_exo_index = BAM_ALLOCATE_MULTIMAPPERS_EXO.out.bai
            ch_exo_allocated_flagstat = BAM_ALLOCATE_MULTIMAPPERS_EXO.out.flagstat
            ch_exo_allocated_stats = BAM_ALLOCATE_MULTIMAPPERS_EXO.out.stats
            ch_exo_allocated_idxstats = BAM_ALLOCATE_MULTIMAPPERS_EXO.out.idxstats
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
    BAM_SHIFT_READS (
        ch_filtered_bam.atacseq.join(ch_filtered_index.atacseq, by: 0),
        ch_fasta.first()
    )
    ch_filtered_bam = ch_filtered_bam.other.mix(BAM_SHIFT_READS.out.bam)
    ch_filtered_index = ch_filtered_index.other.mix(BAM_SHIFT_READS.out.bai)
    ch_versions = ch_versions.mix(BAM_SHIFT_READS.out.versions)

    //
    // MODULE: Final filtering of BAM file with SAMBAMBA (quality filtering)
    //
    // TODO: the same blacklist is used for both the endogenous and exogenous BAM files
    BAM_FILTER_SAMBAMBA_FLT3 (
        ch_filtered_bam.join(ch_filtered_index, by: 0),
        ch_filtered_bed.first(),
        ch_fasta.first()
    )
    ch_filtered_bam = BAM_FILTER_SAMBAMBA_FLT3.out.bam
    ch_filtered_index = BAM_FILTER_SAMBAMBA_FLT3.out.bai
    ch_filtered3_stats = BAM_FILTER_SAMBAMBA_FLT3.out.stats
    ch_filtered3_flagstat = BAM_FILTER_SAMBAMBA_FLT3.out.flagstat
    ch_filtered3_idxstats = BAM_FILTER_SAMBAMBA_FLT3.out.idxstats
    ch_versions = ch_versions.mix(BAM_FILTER_SAMBAMBA_FLT3.out.versions)

     //
    // MODULE: Extract total mapped reads from flagstats
    //
    BAM_FLAGSTAT_MAPPED_FLT3 (
        ch_filtered3_flagstat
    )
    ch_versions = ch_versions.mix(BAM_FLAGSTAT_MAPPED_FLT3.out.versions)

    // Extract the total mapped reads from the text file
    BAM_FLAGSTAT_MAPPED_FLT3.out.txt
        .map {
            meta, total ->
                [ meta, total.splitCsv(header:false)[0][0] ]
        }
        .set { ch_flT3_total }

    // Add the total_mapped_reads to the bams' metas
    ch_filtered_bam
        .combine(ch_filtered_index, by: 0)
        .map {
            meta, bam, bai ->
                [ meta, bam, bai ]
        }
        .combine(ch_flT3_total, by: 0)
        .map {
            meta, bam, bai, total ->
                meta_clone = meta.clone()
                meta_clone.flT3_total_mapped_reads = total.toDouble()
                [ meta_clone, bam, bai ]
        }
        .set { ch_filtered_bam_bai }

    //
    // MODULE: Picard post alignment QC
    //
    // TODO: using first() to convert the tuple to a value channel and make it consumable
    ch_picardcollectmultiplemetrics_multiqc = Channel.empty()
    if (!params.skip_picard_metrics) {
        PICARD_COLLECTMULTIPLEMETRICS (
            ch_filtered_bam.join(ch_filtered_index, by: [0]),
            ch_fasta.first(),
            ch_fai.first()
        )
        ch_picardcollectmultiplemetrics_multiqc = PICARD_COLLECTMULTIPLEMETRICS.out.metrics
        ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions.first())
    }

    //
    // MODULE: Phantompeaktools strand cross-correlation and QC metrics
    //
    ch_phantompeakqualtools_spp_multiqc                 = Channel.empty()
    ch_multiqc_phantompeakqualtools_nsc_multiqc         = Channel.empty()
    ch_multiqc_phantompeakqualtools_rsc_multiqc         = Channel.empty()
    ch_multiqc_phantompeakqualtools_correlation_multiqc = Channel.empty()
    if (!params.skip_spp) {
        PHANTOMPEAKQUALTOOLS (
            ch_filtered_bam
        )
        ch_phantompeakqualtools_spp_multiqc           = PHANTOMPEAKQUALTOOLS.out.spp
        ch_versions = ch_versions.mix(PHANTOMPEAKQUALTOOLS.out.versions.first())

        //
        // MODULE: MultiQC custom content for Phantompeaktools
        //
        MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS (
            PHANTOMPEAKQUALTOOLS.out.spp.join(PHANTOMPEAKQUALTOOLS.out.rdata, by: [0]),
            ch_spp_nsc_header,
            ch_spp_rsc_header,
            ch_spp_correlation_header
        )
        ch_multiqc_phantompeakqualtools_nsc_multiqc         = MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS.out.nsc
        ch_multiqc_phantompeakqualtools_rsc_multiqc         = MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS.out.rsc
        ch_multiqc_phantompeakqualtools_correlation_multiqc = MULTIQC_CUSTOM_PHANTOMPEAKQUALTOOLS.out.correlation
    }

    //
    // SUBWORKFLOW: Normalized bigWig coverage tracks
    //
    BAM_NORMALIZE_BIGWIG_DEEPTOOLS (
        ch_filtered_bam_bai,
        ch_chrom_sizes_endo.first(),
        ch_chrom_sizes_exo.first(),
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


    ch_deeptoolsplotprofile_multiqc = Channel.empty()
    if (!params.skip_plot_profile) {
        //
        // MODULE: deepTools matrix generation for plotting
        //
        DEEPTOOLS_COMPUTEMATRIX (
            BAM_NORMALIZE_BIGWIG_DEEPTOOLS.out.bigwig_binsize1,
            ch_gene_bed.map{ it[1] }.first()
        )
        ch_versions = ch_versions.mix(DEEPTOOLS_COMPUTEMATRIX.out.versions.first())

        //
        // MODULE: deepTools profile plots
        //
        DEEPTOOLS_PLOTPROFILE (
            DEEPTOOLS_COMPUTEMATRIX.out.matrix
        )
        ch_deeptoolsplotprofile_multiqc = DEEPTOOLS_PLOTPROFILE.out.table
        ch_versions = ch_versions.mix(DEEPTOOLS_PLOTPROFILE.out.versions.first())

        //
        // MODULE: deepTools heatmaps
        //
        DEEPTOOLS_PLOTHEATMAP (
            DEEPTOOLS_COMPUTEMATRIX.out.matrix
        )
        ch_versions = ch_versions.mix(DEEPTOOLS_PLOTHEATMAP.out.versions.first())
    }

    // Here we remove the exogenous samples from the filtered_bam_bai channel
    ch_filtered_bam = ch_filtered_bam.filter { it[0].genome == params.genome }
    ch_filtered_index = ch_filtered_index.filter { it[0].genome == params.genome }
    ch_filtered_bam_bai = ch_filtered_bam_bai.filter { it[0].genome == params.genome }

    ch_ds_stats = Channel.empty()
    ch_ds_flagstat = Channel.empty()
    ch_ds_idxstats = Channel.empty()
    if (params.bam_downsampling_method) {
        //
        // SUBWORKFLOW: Downsample IP and control BAM files
        //
        BAM_DOWNSAMPLE (
            ch_filtered_bam_bai,
            ch_fasta.first(),
            ch_fai.first(),
            params.bam_downsampling_method
        )
        ch_filtered_bam = BAM_DOWNSAMPLE.out.bam
        ch_filtered_index = BAM_DOWNSAMPLE.out.bai
        ch_ds_stats = BAM_DOWNSAMPLE.out.stats
        ch_ds_flagstat = BAM_DOWNSAMPLE.out.flagstat
        ch_ds_idxstats = BAM_DOWNSAMPLE.out.idxstats
        ch_versions = ch_versions.mix(BAM_DOWNSAMPLE.out.versions.first())
    }
    
    // Branch channels based on if input control is present
    ch_filtered_bam
        .join(ch_filtered_index, by: 0)
        .branch { meta, bam, bai ->
            ips_with_control: meta.control
                return [ meta.control, meta, [ bam ], [ bai ] ]
            ips_wo_control: !meta.control && !meta.is_control
                return [ meta.id, meta, [ bam ], [ bai ] ]
            controls: !meta.control && meta.is_control
                return [ meta.id, [ bam ], [ bai ] ]
        }
        .set { ch_bam_by_type }

    // Create channel for Consenrich: [ meta, [ip_bams_merged_reps], [ip_bais_merged_reps], [control_bams_merged_reps], [control_bais_merged_reps] ]
    ch_bam_by_type.ips_with_control
        .combine(ch_bam_by_type.controls, by: 0)
        .mix(ch_bam_by_type.ips_wo_control)
        .map { control_id, ip_meta, ip_bam, ip_bai, control_bam, control_bai ->
            def meta_clone = ip_meta.clone()
            meta_clone.id = meta_clone.id - ~/_REP\d+$/
            meta_clone.control = meta_clone.control - ~/_REP\d+$/
            [ meta_clone.id, meta_clone, ip_bam, ip_bai, control_bam ?: [], control_bai ?: [] ]
        }
        .groupTuple()
        .map {
            id, metas, ip_bams, ip_bais, control_bams, control_bais ->
                [ metas[0], ip_bams.flatten(), ip_bais.flatten(), control_bams.flatten(), control_bais.flatten() ]
        }
        .set { ch_ip_control_bam_bai_merged_reps }

    // TODO: Print to file for debuggin
    ch_ip_control_bam_bai_merged_reps
        .map {
            meta, bams, bais, control_bams, control_bais ->
                "${meta}\t${bams}\t${bais}\t${control_bams}\t${control_bais}"
        }
        .collectFile( name: 'ch_ip_control_bam_bai_merged_reps.txt', newLine: true, sort: false, storeDir: "${params.outdir}" )
        
    //
    // MODULE: Call consensus regions with Consenrich
    //
    // CONSENRICH (
    //     ch_ip_control_bam_bai_merged_reps,
    //     ch_chrom_sizes.map{ it[1] }.first(),
    //     ch_blacklist.map{ it[1] }.first(),
    //     ch_sparsebed.map{ it[1] }.first(),
    //     []
    // )
    // ch_versions = ch_versions.mix(CONSENRICH.out.versions.first())
    
    //
    // Create channel for deepTools plotFingerprint: [ meta, [ ip_bam, control_bam ] [ ip_bai, control_bai ] ]
    //
    ch_bam_by_type.ips_with_control
        .combine(ch_bam_by_type.controls, by: 0)
        .mix(ch_bam_by_type.ips_wo_control)
        .map { it -> [ it[1], it[2] + (it[4] ?: []), it[3] + (it[5] ?: []) ] }
        .set { ch_ip_control_bam_bai }
    
    // TODO: Print to file for debuggin
    ch_ip_control_bam_bai
        .map {
            meta, bams, bais ->
                "${meta}\t${bams}\t${bais}"
        }
        .collectFile( name: 'ch_ip_control_bam_bai.txt', newLine: true, sort: false, storeDir: "${params.outdir}" )

    //
    // MODULE: deepTools plotFingerprint joint QC for IP and control
    //
    ch_deeptoolsplotfingerprint_multiqc = Channel.empty()
    if (!params.skip_plot_fingerprint) {
        DEEPTOOLS_PLOTFINGERPRINT (
            ch_ip_control_bam_bai
        )
        ch_deeptoolsplotfingerprint_multiqc = DEEPTOOLS_PLOTFINGERPRINT.out.matrix
        ch_versions = ch_versions.mix(DEEPTOOLS_PLOTFINGERPRINT.out.versions.first())
    }

    // Create channels: [ meta, ip_bam, control_bam ]
    ch_ip_control_bam_bai
        .map {
            meta, bams, bais ->
                [ meta , bams[0], (bams[1] ?: []) ]
        }
        .set { ch_ip_control_bam }

    
    // separate samples based on meta.exp_type
    ch_ip_control_bam_cs = Channel.empty()
    ch_ip_control_bam_cs = ch_ip_control_bam.filter { it[0].exp_type != 'scarseq' }

    // TODO: Print to file for debuggin
    ch_ip_control_bam_cs
        .map {
            meta, ip_bam, control_bam ->
                "${meta.id}\t${ip_bam}\t${control_bam}"
        }
        .collectFile( name: 'ch_ip_control_bam_cs.txt', newLine: true, sort: false, storeDir: "${params.outdir}" )


    //
    // MODULE: Calculate genome size with khmer
    //

    // TODO: genome size is calculated with khmer even when not needed (no chipseq samples)
    // this is an ugly workaround (https://github.com/nextflow-io/nextflow/discussions/5102#discussioncomment-9939140)
    ch_macs_gsize                     = Channel.empty()
    ch_custompeaks_frip_multiqc       = Channel.empty()
    ch_custompeaks_count_multiqc      = Channel.empty()
    ch_plothomerannotatepeaks_multiqc = Channel.empty()
    ch_subreadfeaturecounts_multiqc   = Channel.empty()
    if (!params.macs_gsize) { // && need_macs_gsize) {
        KHMER_UNIQUEKMERS (
            ch_fasta,
            params.read_length
        )
        ch_macs_gsize = KHMER_UNIQUEKMERS.out.kmers.map { it[1].text.trim() }
    }

    //
    // SUBWORKFLOW: Call peaks with MACS3, annotate with HOMER and perform downstream QC
    //
    BAM_PEAKS_CALL_QC_ANNOTATE_MACS3_HOMER (
        ch_ip_control_bam_cs,
        ch_fasta.map{ it[1] }.first(),
        ch_gtf.map{ it[1] }.first(),
        ch_chrom_sizes_endo.first(),
        ch_blacklist.first(),
        ch_macs_gsize.first(),
        "_peaks.annotatePeaks.txt", // TODO: check if this is correct
        ch_peak_count_header,
        ch_frip_score_header,
        ch_peak_annotation_header,
        params.narrow_peak,
        params.skip_peak_annotation,
        params.skip_peak_qc,
        params.skip_edd
    )
    ch_versions = ch_versions.mix(BAM_PEAKS_CALL_QC_ANNOTATE_MACS3_HOMER.out.versions)

    //
    //  Consensus peaks analysis
    //
    ch_macs3_consensus_bed_lib   = Channel.empty()
    ch_macs3_consensus_txt_lib   = Channel.empty()
    ch_deseq2_pca_multiqc        = Channel.empty()
    ch_deseq2_clustering_multiqc = Channel.empty()
    if (!params.skip_consensus_peaks) {
        // Create channels: [ antibody, [ ip_bams ] ]
        ch_ip_control_bam_cs
            .map {
                meta, ip_bam, control_bam ->
                    [ meta.antibody, ip_bam ]
            }
            .groupTuple()
            .set { ch_antibody_bams }

        BED_CONSENSUS_QUANTIFY_QC_BEDTOOLS_FEATURECOUNTS_DESEQ2 (
            BAM_PEAKS_CALL_QC_ANNOTATE_MACS3_HOMER.out.peaks,
            ch_antibody_bams,
            ch_fasta.map{ it[1] },
            ch_gtf.map{ it[1] },
            ch_deseq2_pca_header,
            ch_deseq2_clustering_header,
            params.narrow_peak,
            params.skip_peak_annotation,
            params.skip_deseq2_qc
        )
        ch_macs3_consensus_bed_lib       = BED_CONSENSUS_QUANTIFY_QC_BEDTOOLS_FEATURECOUNTS_DESEQ2.out.consensus_bed
        ch_macs3_consensus_txt_lib       = BED_CONSENSUS_QUANTIFY_QC_BEDTOOLS_FEATURECOUNTS_DESEQ2.out.consensus_txt
        ch_subreadfeaturecounts_multiqc  = BED_CONSENSUS_QUANTIFY_QC_BEDTOOLS_FEATURECOUNTS_DESEQ2.out.featurecounts_summary
        ch_deseq2_pca_multiqc            = BED_CONSENSUS_QUANTIFY_QC_BEDTOOLS_FEATURECOUNTS_DESEQ2.out.deseq2_qc_pca_multiqc
        ch_deseq2_clustering_multiqc     = BED_CONSENSUS_QUANTIFY_QC_BEDTOOLS_FEATURECOUNTS_DESEQ2.out.deseq2_qc_dists_multiqc
        ch_versions = ch_versions.mix(BED_CONSENSUS_QUANTIFY_QC_BEDTOOLS_FEATURECOUNTS_DESEQ2.out.versions)
    }

    //
    // SUBWORKFLOW: Call peaks with Genrich, annotate with HOMER and perform downstream QC
    //
    ch_genrich_frip_multiqc = Channel.empty()
    ch_genrich_peak_count_multiqc = Channel.empty()
    ch_genrich_plot_homer_annotatepeaks_tsv = Channel.empty()
    if (!params.skip_genrich) {
        BAM_PEAKS_CALL_QC_ANNOTATE_GENRICH_HOMER (
            ch_filtered_bam.filter { it[0].exp_type != 'scarseq' },
            ch_fasta.first(),
            ch_gtf.map{ it[1] }.first(),
            ch_blacklist.map{ it[1] }.first(),
            ".annotatePeaks.txt",
            ch_gr_peak_count_header,
            ch_gr_frip_score_header,
            ch_gr_peak_annotation_header,
            params.narrow_peak,
            params.skip_peak_annotation,
            params.skip_peak_qc
        )
        ch_genrich_frip_multiqc = BAM_PEAKS_CALL_QC_ANNOTATE_GENRICH_HOMER.out.frip_multiqc
        ch_genrich_peak_count_multiqc = BAM_PEAKS_CALL_QC_ANNOTATE_GENRICH_HOMER.out.peak_count_multiqc
        ch_genrich_plot_homer_annotatepeaks_tsv = BAM_PEAKS_CALL_QC_ANNOTATE_GENRICH_HOMER.out.plot_homer_annotatepeaks_tsv
        ch_versions = ch_versions.mix(BAM_PEAKS_CALL_QC_ANNOTATE_GENRICH_HOMER.out.versions)
    }


    ch_filtered_bam_ss = Channel.empty()
    ch_filtered_bam_ss = ch_filtered_bam.filter { it[0].exp_type == 'scarseq' }

    // Make ch_chrom_sizes_endo empty if there are no scarseq samples
    // This is to avoid unnecessarily running modules in the next subworkflow
    ch_chrom_sizes_endo
        .combine(ch_filtered_bam_ss)
        .first()
        .map {
            sizes_meta, sizes, ss_meta, ss_bam ->
            [sizes_meta, sizes]
        }
        .set { ch_chrom_sizes_endo_ss }

    //
    // SUBWORKFLOW: SCAR-seq analysis: partitioning of reads
    //
    //https://github.com/nextflow-io/nextflow/issues/1052
    ch_scar_smooth = Channel.empty()
    BAM_CREATE_SCAR_PARTITIONS (
        ch_filtered_bam_ss,
        ch_chrom_sizes_endo_ss.first(),
        ch_blacklist.first(),
        ch_initiation_zones.first()
        //ch_scaffolds
    )
    ch_scar_smooth = BAM_CREATE_SCAR_PARTITIONS.out.tab
    ch_versions = ch_versions.mix(BAM_CREATE_SCAR_PARTITIONS.out.versions)

    // Create channel containing all samtools_stats files
    ch_samtools_stats
        .mix(ch_merged_bam_stats)
        .mix(ch_dedup_umi_stats)
        .mix(ch_mk_stats)
        .mix(ch_filtered_stats)
        .mix(ch_filtered2_stats)
        .mix(ch_filtered2_exo_stats)
        .mix(ch_allocated_stats)
        .mix(ch_filtered3_stats)
        .mix(ch_ds_stats)
        .set { ch_samtools_stats_summary }

    //
    // SUBWORKFLOW: Create SAMtools summary table
    //
    SAMTOOLS_STATS_SUMMARY (
        ch_samtools_stats_summary,
        params.genome,
        params.spikein_genome ?: Channel.of([])
    )
    ch_versions = ch_versions.mix(SAMTOOLS_STATS_SUMMARY.out.versions)

    //
    // MODULE: Create IGV session
    //
    if (!params.skip_igv) {
        IGV (
            params.aligner,
            params.allocate_n_multimappers ? params.allocation_method == 'chromap' ? 'chromap_allocation' : params.allocation_method + '/' : '',
            params.narrow_peak ? 'narrow_peak' : 'broad_peak',
            ch_fasta.map{ it[1] },
            BAM_NORMALIZE_BIGWIG_DEEPTOOLS.out.bigwig_endo.collect{it[1]}.ifEmpty([]),
            BAM_PEAKS_CALL_QC_ANNOTATE_MACS3_HOMER.out.peaks.collect{it[1]}.ifEmpty([]),
            ch_macs3_consensus_bed_lib.collect{it[1]}.ifEmpty([]),
            ch_macs3_consensus_txt_lib.collect{it[1]}.ifEmpty([])
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
        ch_multiqc_config        = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
        ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath( params.multiqc_config ): Channel.empty()
        ch_multiqc_logo          = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo )  : Channel.empty()
        summary_params           = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
        ch_workflow_summary      = Channel.value(paramsSummaryMultiqc(summary_params))
        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)

        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList(),

            FASTQ_FASTQC_UMITOOLS_UMITRANSFER_TRIMGALORE.out.fastqc_zip.collect{it[1]}.ifEmpty([]),
            FASTQ_FASTQC_UMITOOLS_UMITRANSFER_TRIMGALORE.out.trim_zip.collect{it[1]}.ifEmpty([]),
            FASTQ_FASTQC_UMITOOLS_UMITRANSFER_TRIMGALORE.out.trim_log.collect{it[1]}.ifEmpty([]),

            ch_samtools_stats.collect{it[1]}.ifEmpty([]),
            ch_samtools_flagstat.collect{it[1]}.ifEmpty([]),
            ch_samtools_idxstats.collect{it[1]}.ifEmpty([]),

            ch_merged_bam_stats.collect{it[1]}.ifEmpty([]),
            ch_merged_bam_flagstat.collect{it[1]}.ifEmpty([]),
            ch_merged_bam_idxstats.collect{it[1]}.ifEmpty([]),

            ch_preseq_multiqc.collect{it[1]}.ifEmpty([]),

            ch_dedup_umi_stats.collect{it[1]}.ifEmpty([]),
            ch_dedup_umi_flagstat.collect{it[1]}.ifEmpty([]),
            ch_dedup_umi_idxstats.collect{it[1]}.ifEmpty([]),
            ch_dedup_umi_deduplog.collect{it[1]}.ifEmpty([]),

            ch_mk_stats.collect{it[1]}.ifEmpty([]),
            ch_mk_flagstat.collect{it[1]}.ifEmpty([]),
            ch_mk_idxstats.collect{it[1]}.ifEmpty([]),
            ch_mk_metrics.collect{it[1]}.ifEmpty([]),

            ch_filtered_stats.collect{it[1]}.ifEmpty([]),
            ch_filtered_flagstat.collect{it[1]}.ifEmpty([]),
            ch_filtered_idxstats.collect{it[1]}.ifEmpty([]),

            ch_filtered2_stats.collect{it[1]}.ifEmpty([]),
            ch_filtered2_flagstat.collect{it[1]}.ifEmpty([]),
            ch_filtered2_idxstats.collect{it[1]}.ifEmpty([]),

            ch_allocated_flagstat.collect{it[1]}.ifEmpty([]),
            ch_allocated_stats.collect{it[1]}.ifEmpty([]),
            ch_allocated_idxstats.collect{it[1]}.ifEmpty([]),

            ch_picardcollectmultiplemetrics_multiqc.collect{it[1]}.ifEmpty([]),


            ch_deeptoolsplotprofile_multiqc.collect{it[1]}.ifEmpty([]),
            ch_deeptoolsplotfingerprint_multiqc.collect{it[1]}.ifEmpty([]),

            ch_phantompeakqualtools_spp_multiqc.collect{it[1]}.ifEmpty([]),
            ch_multiqc_phantompeakqualtools_nsc_multiqc.collect{it[1]}.ifEmpty([]),
            ch_multiqc_phantompeakqualtools_rsc_multiqc.collect{it[1]}.ifEmpty([]),
            ch_multiqc_phantompeakqualtools_correlation_multiqc.collect{it[1]}.ifEmpty([]),

            BAM_PEAKS_CALL_QC_ANNOTATE_MACS3_HOMER.out.frip_multiqc.collect{it[1]}.ifEmpty([]),
            BAM_PEAKS_CALL_QC_ANNOTATE_MACS3_HOMER.out.peak_count_multiqc.collect{it[1]}.ifEmpty([]),
            BAM_PEAKS_CALL_QC_ANNOTATE_MACS3_HOMER.out.plot_homer_annotatepeaks_tsv.collect().ifEmpty([]),
            ch_subreadfeaturecounts_multiqc.collect{it[1]}.ifEmpty([]),

            ch_deseq2_pca_multiqc.collect().ifEmpty([]),
            ch_deseq2_clustering_multiqc.collect().ifEmpty([]),

            ch_genrich_frip_multiqc.collect{it[1]}.ifEmpty([]),
            ch_genrich_peak_count_multiqc.collect{it[1]}.ifEmpty([]),
            ch_genrich_plot_homer_annotatepeaks_tsv.collect().ifEmpty([])
            
        )
        ch_multiqc_report = MULTIQC.out.report
    }

    emit:
    multiqc_report = ch_multiqc_report  // channel: /path/to/multiqc_report.html
    versions       = ch_versions       // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/