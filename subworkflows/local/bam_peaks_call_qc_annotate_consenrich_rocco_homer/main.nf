include { FILE_SORT as SIZES_SORT } from '../../../modules/local/file_sort/main'
include { CONSENRICH           } from '../../../modules/local/consenrich/main'
include { ROCCO               } from '../../../modules/local/rocco/main'


workflow BAM_PEAKS_CALL_QC_ANNOTATE_CONSENRICH_ROCCO_HOMER {
    take:
    ch_bam_bai            // channel: [ val(meta), [ bam ], [ bai ] ]
    ch_chrom_sizes        // channel: [ val(meta), [ chrom_sizes ] ]
    ch_blacklist          // channel: [ val(meta), [ blacklist ] ]
    ch_sparsebed          // channel: [ val(meta), [ sparsebed ] ]
    ch_active_regions     // channel: [ val(meta), [ active_regions ] ]
    ch_rocco_params      // channel: [ val(meta), path(rocco_params) ]
    ch_effective_gsize

    main:

    //
    // MODULE: Sort chromosome sizes (avoid Consenrich errors)
    //
    SIZES_SORT (
        ch_chrom_sizes,
        'sizes'
    )
    ch_chrom_sizes = SIZES_SORT.out.sorted


    // Branch channels based on if input control is present
    ch_bam_bai
        .branch { meta, bam, bai ->
            ips_with_ipcontrol: meta.input_control
                return [meta.input_control, meta.antibody, meta, bam, bai]
            ips_wo_ipcontrol: !meta.input_control && !meta.is_input_control
                return [meta.id, meta.antibody, meta, bam, bai]
            ipcontrols: !meta.input_control && meta.is_input_control
                return [meta.id, meta.input_control_of_antibody, meta, bam, bai]
        }
        .set { ch_bam_by_type }

    // Create channel for Consenrich: [ meta, [ip_bams_merged_reps], [ip_bais_merged_reps], [ipcontrol_bams_merged_reps], [ipcontrol_bais_merged_reps] ]
    ch_bam_by_type
        .ips_with_ipcontrol
        .combine(ch_bam_by_type.ipcontrols, by: [0, 1])
        .map { ipcontrol_id, antibody, ip_meta, ip_bam, ip_bai, ipcontrol_meta, ipcontrol_bam, ipcontrol_bai ->
            [ ipcontrol_id, antibody, ip_meta, ip_bam, ip_bai, ipcontrol_bam, ipcontrol_bai ]
        }
        .mix(ch_bam_by_type.ips_wo_ipcontrol)
        .map { it ->
            def meta_clone = it[2].clone()
            meta_clone.id = meta_clone.id - ~/_bRep_.*$/
            meta_clone.input_control = meta_clone.input_control - ~/_bRep_.*$/
            // ips_wo_ipcontrol do not have ipcontrol_bam (it[5]) and ipcontrol_bai (it[6])
            [meta_clone.id, it[1], meta_clone, it[3], it[4], it[5] ?: [], it[6] ?: []]
        }
        .groupTuple(by: [0, 1])
        .map { id, antibody, metas, ip_bams, ip_bais, ipcontrol_bams, ipcontrol_bais ->
            [metas[0], ip_bams.flatten(), ip_bais.flatten(), ipcontrol_bams.flatten(), ipcontrol_bais.flatten()]
        }
        .set { ch_ip_ipcontrol_bam_bai_merged_reps }

    // TODO: Print to file for debuggin
    ch_ip_ipcontrol_bam_bai_merged_reps
        .map { meta, bams, bais, ipcontrol_bams, ipcontrol_bais ->
            "${meta}\t${bams}\t${bais}\t${ipcontrol_bams}\t${ipcontrol_bais}"
        }
        .collectFile(name: 'ch_ip_ipcontrol_bam_bai_merged_reps.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_PEAKS_CALL_QC_ANNOTATE_CONSENRICH_HOMER" )


    //
    // MODULE: Integrate ChIPs and their input ipcontrols
    //
    CONSENRICH (
        ch_ip_ipcontrol_bam_bai_merged_reps,
        ch_chrom_sizes,
        ch_blacklist,
        ch_sparsebed,
        ch_active_regions
    )


    // Create channel: [ val(meta), ch_csr_signal, bamlist_txt ]
    CONSENRICH
        .out
        .signal_track
        .combine(ch_ip_ipcontrol_bam_bai_merged_reps, by: 0)
        .map { meta, ch_csr_signal, ip_bams, ip_bais, ipcontrol_bams, ipcontrol_bais ->
            // TODO: check if ip_bams and ipcontrol_bams should be interleaved or
            //       in the same order as in consenrich. Here we just mix them:
            [ meta, ch_csr_signal, ip_bams + ipcontrol_bams, ip_bais + ipcontrol_bais ]
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

    emit:

    consenrich_signal       = CONSENRICH.out.signal_track   // channel: [ val(meta), path(consenrich_signal_track.bw) ]
    consenrich_residuals    = CONSENRICH.out.residuals_track // channel
    consenrich_eratio       = CONSENRICH.out.eratio_track   // channel: [ val(meta), path(consenrich_eratio_track.bw) ]
    rocco_peaks             = ROCCO.out.bed            // channel: [ val(meta), path(rocco_bed) ]
}
