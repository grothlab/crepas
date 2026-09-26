//
// FASTQ_ALIGN: Align fastq files to a reference genome
//

include { FASTQ_ALIGN_BWA as FASTQ_ALIGN_BWAMEM1     } from '../../../subworkflows/nf-core/fastq_align_bwa'
include { FASTQ_ALIGN_BWAMEM2                        } from '../../../subworkflows/local/fastq_align_bwamem2'
include { FASTQ_ALIGN_MINIBWA                        } from '../../../subworkflows/local/fastq_align_minibwa'
include { FASTQ_ALIGN_BOWTIE2                        } from '../../../subworkflows/nf-core/fastq_align_bowtie2'
include { FASTQ_ALIGN_BOWTIE                         } from '../../../subworkflows/local/fastq_align_bowtie'
include { FASTQ_ALIGN_STROBEALIGN                    } from '../../../subworkflows/local/fastq_align_strobealign'
include { FASTQ_ALIGN_CHROMAP                        } from '../../../subworkflows/nf-core/fastq_align_chromap'
include { FASTQ_ALIGN_STAR                           } from '../../../subworkflows/nf-core/fastq_align_star'
include { FASTQ_ALIGN_HISAT2                         } from '../../../subworkflows/nf-core/fastq_align_hisat2'
include { FASTQ_ALIGN_MINIMAP2                       } from '../../../subworkflows/local/fastq_align_minimap2'

workflow FASTQ_ALIGN {
    take:
    ch_reads            // channel: [mandatory] meta, reads
    ch_fasta_fai
    aligner
    ch_bwa_index
    ch_bwamem2_index
    ch_minibwa_index
    ch_bowtie_index
    ch_bowtie2_index
    ch_chromap_index
    ch_star_index
    ch_hisat2_index
    ch_minimap2_index
    ch_strobealign_index
    ch_gtf
    ch_splicesites
    save_unaligned
    val_sort_bam        // boolean: always false -> aligners emit unsorted BAM, sorted downstream by BAM_SORT_STATS_SAMTOOLS

    main:

    ch_multiqc_files = channel.empty()
    ch_samtools_stats_summary = channel.empty()


    if (aligner == 'bwa') {
        FASTQ_ALIGN_BWAMEM1 (
            ch_reads,
            ch_bwa_index,
            false,
            ch_fasta_fai
        )
        ch_genome_bam = FASTQ_ALIGN_BWAMEM1.out.bam
        ch_genome_bam_index = FASTQ_ALIGN_BWAMEM1.out.index
        ch_samtools_stats_summary = ch_samtools_stats_summary.mix(FASTQ_ALIGN_BWAMEM1.out.stats)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_BWAMEM1.out.stats.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_BWAMEM1.out.flagstat.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_BWAMEM1.out.idxstats.collect { it -> it[1] })
    }

    if (aligner == 'bwamem2') {
        FASTQ_ALIGN_BWAMEM2 (
            ch_reads,
            ch_bwamem2_index,
            val_sort_bam,
            ch_fasta_fai
        )
        ch_genome_bam = FASTQ_ALIGN_BWAMEM2.out.bam
        ch_genome_bam_index = FASTQ_ALIGN_BWAMEM2.out.index
        ch_samtools_stats_summary = ch_samtools_stats_summary.mix(FASTQ_ALIGN_BWAMEM2.out.stats)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_BWAMEM2.out.stats.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_BWAMEM2.out.flagstat.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_BWAMEM2.out.idxstats.collect { it -> it[1] })
    }

    if (aligner == 'minibwa') {
        FASTQ_ALIGN_MINIBWA (
            ch_reads,
            ch_minibwa_index,
            val_sort_bam,
            ch_fasta_fai
        )
        ch_genome_bam = FASTQ_ALIGN_MINIBWA.out.bam
        ch_genome_bam_index = FASTQ_ALIGN_MINIBWA.out.index
        ch_samtools_stats_summary = ch_samtools_stats_summary.mix(FASTQ_ALIGN_MINIBWA.out.stats)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_MINIBWA.out.stats.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_MINIBWA.out.flagstat.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_MINIBWA.out.idxstats.collect { it -> it[1] })
    }


    //
    // SUBWORKFLOW: Alignment with Bowtie & BAM QC
    //
    if (aligner == 'bowtie') {
        FASTQ_ALIGN_BOWTIE (
            ch_reads,
            ch_bowtie_index,
            save_unaligned,
            ch_fasta_fai
        )
        ch_genome_bam = FASTQ_ALIGN_BOWTIE.out.bam
        ch_genome_bam_index = FASTQ_ALIGN_BOWTIE.out.index
        ch_samtools_stats_summary = ch_samtools_stats_summary.mix(FASTQ_ALIGN_BOWTIE.out.stats)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_BOWTIE.out.stats.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_BOWTIE.out.flagstat.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_BOWTIE.out.idxstats.collect { it -> it[1] })
    }

    //
    // SUBWORKFLOW: Alignment with Bowtie2 & BAM QC
    //
    if (aligner == 'bowtie2') {
        FASTQ_ALIGN_BOWTIE2 (
            ch_reads,
            ch_bowtie2_index,
            save_unaligned,
            val_sort_bam,
            ch_fasta_fai
        )
        ch_genome_bam = FASTQ_ALIGN_BOWTIE2.out.bam
        ch_genome_bam_index = FASTQ_ALIGN_BOWTIE2.out.index
        ch_samtools_stats_summary = ch_samtools_stats_summary.mix(FASTQ_ALIGN_BOWTIE2.out.stats)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_BOWTIE2.out.stats.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_BOWTIE2.out.flagstat.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_BOWTIE2.out.idxstats.collect { it -> it[1] })
    }

    //
    // SUBWORKFLOW: Alignment with strobealign
    //
    if (aligner == 'strobealign') {
        FASTQ_ALIGN_STROBEALIGN (
            ch_reads,
            ch_strobealign_index,
            ch_fasta_fai,
            val_sort_bam
        )
        ch_genome_bam = FASTQ_ALIGN_STROBEALIGN.out.bam
        ch_genome_bam_index = FASTQ_ALIGN_STROBEALIGN.out.index
        ch_samtools_stats_summary = ch_samtools_stats_summary.mix(FASTQ_ALIGN_STROBEALIGN.out.stats)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_STROBEALIGN.out.stats.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_STROBEALIGN.out.flagstat.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_STROBEALIGN.out.idxstats.collect { it -> it[1] })
    }

    //
    // SUBWORKFLOW: Alignment with Chromap & BAM QC
    //
    if (aligner == 'chromap') {
        FASTQ_ALIGN_CHROMAP (
            ch_reads,
            ch_chromap_index,
            ch_fasta_fai,
            [],
            [],
            [],
            [],
            true
        )
        ch_genome_bam = FASTQ_ALIGN_CHROMAP.out.bam
        ch_genome_bam_index = FASTQ_ALIGN_CHROMAP.out.index
        ch_samtools_stats_summary = ch_samtools_stats_summary.mix(FASTQ_ALIGN_CHROMAP.out.stats)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_CHROMAP.out.stats.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_CHROMAP.out.flagstat.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_CHROMAP.out.idxstats.collect { it -> it[1] })
    }

    //
    // SUBWORKFLOW: Alignment with STAR & BAM QC
    //
    if (aligner == 'star') {
        FASTQ_ALIGN_STAR (
            ch_reads,
            ch_star_index,
            ch_gtf,
            true,
            ch_fasta_fai,
            channel.value([[:], []])
        )
        ch_genome_bam = FASTQ_ALIGN_STAR.out.bam
        ch_genome_bam_index = FASTQ_ALIGN_STAR.out.index
        ch_samtools_stats_summary = ch_samtools_stats_summary.mix(FASTQ_ALIGN_STAR.out.stats)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_STAR.out.stats.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_STAR.out.flagstat.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_STAR.out.idxstats.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_STAR.out.log_final.collect { it -> it[1] })
    }

    if (aligner == 'hisat2') {
        FASTQ_ALIGN_HISAT2 (
            ch_reads,
            ch_hisat2_index,
            ch_splicesites,
            ch_fasta_fai,
            save_unaligned
        )
        ch_genome_bam = FASTQ_ALIGN_HISAT2.out.bam
        ch_genome_bam_index = FASTQ_ALIGN_HISAT2.out.index
        ch_samtools_stats_summary = ch_samtools_stats_summary.mix(FASTQ_ALIGN_HISAT2.out.stats)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_HISAT2.out.stats.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_HISAT2.out.flagstat.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_HISAT2.out.idxstats.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_HISAT2.out.summary.collect { it -> it[1] })
    }

    if (aligner == 'minimap2') {
        FASTQ_ALIGN_MINIMAP2 (
            ch_reads,
            ch_minimap2_index,
            ch_fasta_fai,
            true,
            'bai',
            false,
            true
        )
        ch_genome_bam = FASTQ_ALIGN_MINIMAP2.out.bam
        ch_genome_bam_index = FASTQ_ALIGN_MINIMAP2.out.index
        ch_samtools_stats_summary = ch_samtools_stats_summary.mix(FASTQ_ALIGN_MINIMAP2.out.stats)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_MINIMAP2.out.stats.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_MINIMAP2.out.flagstat.collect { it -> it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_ALIGN_MINIMAP2.out.idxstats.collect { it -> it[1] })
    }

    emit:
    bam                     = ch_genome_bam                 // channel: [ val(meta), path(bam) ]
    bai                     = ch_genome_bam_index           // channel: [ val(meta), path(bai) ]
    samtools_stats_summary  = ch_samtools_stats_summary     // channel: [ val(meta), path(summary) ]
    multiqc_files           = ch_multiqc_files              // channel: [ path(multiqc files ]
}
