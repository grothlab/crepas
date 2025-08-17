include { FILE_SORT as SIZES_SORT } from '../../../modules/local/file_sort/main'
include { CONSENRICH           } from '../../../modules/local/consenrich/main'
include { ROCCO               } from '../../../modules/local/rocco/main'


workflow BAM_PEAKS_CALL_QC_ANNOTATE_CONSENRICH_HOMER {
    take:
    ch_bam_bai            // channel: [ val(meta), [ bam ], [ bai ] ]
    ch_chrom_sizes        // channel: [ val(meta), [ chrom_sizes ] ]
    ch_blacklist          // channel: [ val(meta), [ blacklist ] ]
    ch_sparsebed          // channel: [ val(meta), [ sparsebed ] ]
    ch_active_regions     // channel: [ val(meta), [ active_regions ] ]
    ch_rocco_params      // channel: [ val(meta), path(rocco_params) ]
    ch_effective_gsize

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: Sort chromosome sizes (avoid Consenrich errors)
    //
    SIZES_SORT (
        ch_chrom_sizes,
        'sizes'
    )
    ch_chrom_sizes = SIZES_SORT.out.sorted
    ch_versions = ch_versions.mix(SIZES_SORT.out.versions.first())


    // Branch channels based on if input control is present
    ch_bam_bai
        .branch { meta, bam, bai ->
            ips_with_control: meta.control
                return [meta.control, meta.antibody, meta, bam, bai]
            ips_wo_control: !meta.control && !meta.is_control
                return [meta.id, meta.antibody, meta, bam, bai]
            controls: !meta.control && meta.is_control
                return [meta.id, meta, bam, bai]
        }
        .set { ch_bam_by_type }

    // For non-downsampled files, duplicate input controls for each antibody 
    ch_bam_by_type
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
        .combine(ch_bam_by_type.ips_with_control, by: 0) // combine by control id only
        .map { control_id, control_meta, control_bam, control_bai, ip_antibody, ip_meta, ip_bam, ip_bai ->
            def meta_clone = control_meta.clone()
            meta_clone.control_of_antibody = ip_antibody
            [ control_id, meta_clone.control_of_antibody, meta_clone, control_bam, control_bai ]
        }
        .unique()
        .set { ch_bam_controls_not_dsp }
    
    // Create channel for Consenrich: [ meta, [ip_bams_merged_reps], [ip_bais_merged_reps], [control_bams_merged_reps], [control_bais_merged_reps] ]
    ch_bam_by_type
        .ips_with_control
        .combine(ch_bam_controls.dsp.mix(ch_bam_controls_not_dsp), by: [0, 1])
        .map { control_id, antibody, ip_meta, ip_bam, ip_bai, control_meta, control_bam, control_bai ->
            [ control_id, antibody, ip_meta, ip_bam, ip_bai, control_bam, control_bai ]
        }
        .mix(ch_bam_by_type.ips_wo_control)
        .map { it ->
            def meta_clone = it[2].clone()
            meta_clone.id = meta_clone.id - ~/_REP\d+$/
            meta_clone.control = meta_clone.control - ~/_REP\d+$/
            // ips_wo_control do not have control_bam (it[5]) and control_bai (it[6])
            [meta_clone.id, it[1], meta_clone, it[3], it[4], it[5] ?: [], it[6] ?: []]
        }
        .groupTuple(by: [0, 1])
        .map { id, antibody, metas, ip_bams, ip_bais, control_bams, control_bais ->
            [metas[0], ip_bams.flatten(), ip_bais.flatten(), control_bams.flatten(), control_bais.flatten()]
        }
        .set { ch_ip_control_bam_bai_merged_reps }

    // TODO: Print to file for debuggin
    ch_ip_control_bam_bai_merged_reps
        .map { meta, bams, bais, control_bams, control_bais ->
            "${meta}\t${bams}\t${bais}\t${control_bams}\t${control_bais}"
        }
        .collectFile(name: 'ch_ip_control_bam_bai_merged_reps.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_PEAKS_CALL_QC_ANNOTATE_CONSENRICH_HOMER" )


    //
    // MODULE: Integrate ChIPs and their input controls
    //
    CONSENRICH (
        ch_ip_control_bam_bai_merged_reps,
        ch_chrom_sizes,
        ch_blacklist,
        ch_sparsebed,
        ch_active_regions
    )
    ch_versions = ch_versions.mix(CONSENRICH.out.versions.first())

    
    // Create channel: [ val(meta), ch_csr_signal, bamlist_txt ]
    CONSENRICH
        .out
        .signal_track
        .combine(ch_ip_control_bam_bai_merged_reps, by: 0)
        .map { meta, ch_csr_signal, ip_bams, ip_bais, control_bams, control_bais ->
            // TODO: check if ip_bams and control_bams should be interleaved or
            //       in the same order as in consenrich. Here we just mix them:
            [ meta, ch_csr_signal, ip_bams + control_bams, ip_bais + control_bais ]
        }
        .set { ch_csr_signal_bamlist }


    // TODO: Print for debugging
    ch_csr_signal_bamlist
        .map { meta, ch_csr_signal, bamlist, bailist ->
            "${meta}\t${ch_csr_signal}\t${bamlist}\t${bailist}"
        }
        .collectFile( name: 'ch_csr_signal_bamlist.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_PEAKS_CALL_QC_ANNOTATE_CONSENRICH_HOMER" )

    //
    // MODULE: Call peaks
    //
    ROCCO (
        ch_csr_signal_bamlist,
        ch_chrom_sizes,
        ch_rocco_params,
        ch_effective_gsize
    )
    ch_versions = ch_versions.mix(ROCCO.out.versions.first())

    emit:

    rocco_bed   = ROCCO.out.bed            // channel: [ val(meta), path(rocco_bed) ]

    versions    = ch_versions              // channel: [ versions.yml ]
}
