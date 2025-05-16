include { PICARD_DOWNSAMPLESAM } from '../../../modules/local/picard/downsamplesam/main'
include { SAMTOOLS_INDEX } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS } from '../../../subworkflows/nf-core/bam_stats_samtools/main'

workflow BAM_DOWNSAMPLE {

    take:
    ch_bam_bai        // channel: [ val(meta), [ bam ] , [ bai ] ]
    ch_fasta          // channel: [ val(meta), [ fasta ] ]
    ch_fai            // channel: [ val(meta), [ fai ] ]
    downsampling_method // string: e.g. 'min_by_type'

    main:

    ch_versions = Channel.empty()

    if (downsampling_method == 'min_by_type') {

        // Split ChIPs and inputs
        ch_bam_bai
            .branch { meta, bam, bai ->
                ips: !meta.is_control
                controls: meta.is_control
            }
            .set { ch_bam_bai_b }

        // Get the minimum number of reads among the ChIPs
        ch_bam_bai_b
            .ips
            .min { it[0].total_mapped_reads }
            .map { meta, bam, bai -> meta.total_mapped_reads }
            .set { ch_min_ips }

        // Get the downsample proportion for each ChIP based on the minimum
        ch_bam_bai_b
            .ips
            .combine(ch_min_ips)
            .map { meta, bam, bai, min ->
                def meta_clone = meta.clone()
                meta_clone.downsampleA_prob = min / meta_clone.total_mapped_reads
                [ meta_clone, bam, bai ]
            }
            .set { ch_bam_bai_ips_ds }

        // Get the minimum number of reads among the inputs
        ch_bam_bai_b
            .controls
            .min { it[0].total_mapped_reads }
            .map { meta, bam, bai -> meta.total_mapped_reads }
            .set { ch_min_controls }

        // Get the downsample proportion for each input based on the minimum
        ch_bam_bai_b
            .controls
            .combine(ch_min_controls)
            .map { meta, bam, bai, min ->
                def meta_clone = meta.clone()
                meta_clone.downsampleA_prob = min / meta_clone.total_mapped_reads
                [ meta_clone, bam, bai ]
            }
            .set { ch_bam_bai_controls_ds }

        // Merge back the downsampled ChIPs and inputs
        ch_bam_bai_ips_ds
            .mix(ch_bam_bai_controls_ds)
            .map { 
                meta, bam, bai ->
                    [ meta, bam ]
            }
            .set { ch_bam_ds }

        // TODO: Remove this: Print channel for debugging
        ch_bam_ds
            .map {
                meta, bam ->
                    "${meta}\t${bam}"
            }
            .collectFile( name: 'ch_bam_ds.txt', newLine: true, sort: false, storeDir: "${params.outdir}" )

    }

    //
    // MODULE: Downsample BAMs
    //
    PICARD_DOWNSAMPLESAM (
        ch_bam_ds,
        ch_fasta,
        ch_fai
    )
    ch_ds_bam = PICARD_DOWNSAMPLESAM.out.bam
    ch_versions = ch_versions.mix(PICARD_DOWNSAMPLESAM.out.versions)

    //
    // MODULE: Index BAMs
    //
    SAMTOOLS_INDEX (
        ch_ds_bam
    )
    ch_ds_index = SAMTOOLS_INDEX.out.bai
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    //
    // SUBWORKFLOW: Run SAMtools stats, flagstat and idxstats
    //
    BAM_STATS_SAMTOOLS (
        ch_ds_bam.join(ch_ds_index, by: [0]),
        ch_fasta
    )
    ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)


    emit:
    bam      = ch_ds_bam   // channel: [ val(meta), [ bam ] ]
    bai      = ch_ds_index           // channel: [ val(meta), [ index ] ]
    stats    = BAM_STATS_SAMTOOLS.out.stats // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]

    versions = ch_versions                    // channel: [ versions.yml ]
}

