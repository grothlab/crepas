
//
// Call peaks with epic2, annotate with HOMER and perform downstream QC
//

include { BAM_REMOVE_SCAFFOLDS      } from '../../../modules/local/bam_remove_scaffolds/main'
include { DANPOS2_DPEAK             } from '../../../modules/local/danpos2/dpeak/main'
include { DANPOS2_DPOS              } from '../../../modules/local/danpos2/dpos/main'
include { DANPOS2_DREGION           } from '../../../modules/local/danpos2/dregion/main'
include { DANPOS2_DTRIPLE           } from '../../../modules/local/danpos2/dtriple/main'

workflow BAM_PEAKS_CALL_QC_ANNOTATE_DANPOS2_HOMER {
    take:
    ch_bam                            // channel: [ val(meta), [ ip_bam ], [ control_bam ] ]

    main:

    //
    // MODULE: Remove scaffolds from BAM files
    //
    BAM_REMOVE_SCAFFOLDS (
        ch_bam
    )

    // Branch channels based on if input control is present
    BAM_REMOVE_SCAFFOLDS.out.bam
        .branch { meta, bam ->
            ips_with_ipcontrol: meta.input_control
                return [meta.input_control, meta.antibody, meta, bam]
            ips_wo_ipcontrol: !meta.input_control && !meta.is_input_control
                return [meta.id, meta.antibody, meta, bam]
            ipcontrols: !meta.input_control && meta.is_input_control
                return [meta.id, meta.input_control_of_antibody, meta, bam]
        }
        .set { ch_bam_by_type }


    // Create channel: [ meta, [ip_bams_merged_reps], [ipcontrol_bams_merged_reps] ]
    ch_bam_by_type
        .ips_with_ipcontrol
        .combine(ch_bam_by_type.ipcontrols, by: [0, 1])
        .map { ipcontrol_id, antibody, ip_meta, ip_bam, ipcontrol_meta, ipcontrol_bam ->
            [ ipcontrol_id, antibody, ip_meta, ip_bam, ipcontrol_bam ]
        }
        .mix(ch_bam_by_type.ips_wo_ipcontrol)
        // ips_wo_ipcontrol do not have control_bam (it[4])
        .map { it ->
            def meta_clone = it[2].clone()
            meta_clone.id = meta_clone.id - ~/_bRep_.*$/
            meta_clone.input_control = meta_clone.input_control - ~/_bRep_.*$/
            [ meta_clone.id, it[1], meta_clone, it[3], it[4] ?: [] ]
        }
        .groupTuple(by: [0, 1])
        .map {
            id, antibody, metas, ip_bams, ipcontrol_bams ->
                // , [], [], [], [] ] is because DANPOS2 modules expect 7 inputs
                [ metas[0], ip_bams.flatten(), ipcontrol_bams.flatten().unique(), [], [], [], [] ]
        }
        .set { ch_ip_control_bam_merged_reps }

    // TODO: Print to file for debuggin
    ch_ip_control_bam_merged_reps
        .map {
            meta, ip_bams, ipcontrol_bams, treatment_count, control, control_input, control_count ->
                "${meta}\t${ip_bams}\t${ipcontrol_bams}\t${treatment_count}\t${control}\t${control_input}\t${control_count}"
        }
        .collectFile( name: 'ch_ip_control_bam_merged_reps.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_PEAKS_CALL_QC_ANNOTATE_DANPOS2_HOMER" )


    ch_dpeak_pooled_xls = channel.empty()
    //
    // MODULE: Call peaks with DANPOS2 dpeak
    //
    DANPOS2_DPEAK (
        ch_ip_control_bam_merged_reps.filter { it -> it[0].peak_caller == 'dpeak' }
    )
    ch_dpeak_pooled_xls = DANPOS2_DPEAK.out.pooled_xls

    //
    // MODULE: call peaks with DANPOS2 dpos
    //
    DANPOS2_DPOS (
        ch_ip_control_bam_merged_reps.filter { it -> it[0].peak_caller == 'dpos' }
    )

    //
    // MODULE: call peaks with DANPOS2 dregion
    //
    DANPOS2_DREGION (
        ch_ip_control_bam_merged_reps.filter { it -> it[0].peak_caller == 'dregion' }
    )

    //
    // MODULE: call peaks with DANPOS2 dtriple
    //
    DANPOS2_DTRIPLE (
        ch_ip_control_bam_merged_reps.filter { it -> it[0].peak_caller == 'dtriple' }
    )


    emit:

    dpeak_pooled_xls    = ch_dpeak_pooled_xls

}
