//
// Call peaks with epic2
//

include { EPIC2         } from '../../../modules/local/epic2/epic2/main'
include { EPIC2_TO_PEAK } from '../../../modules/local/epic2_to_peak/main'

workflow BAM_PEAKS_CALL_EPIC2 {
    take:
    ch_bam_bai                        // channel: [ val(meta), bam, bai ]
    ch_chrom_sizes                    // channel: [ val(meta), path(chrom_sizes) ]
    efective_gfraction                // float: effective genome fraction

    main:

    // Branch channels based on if input control is present
    ch_bam_bai
        .branch { meta, bam, bai ->
            ips_with_ipcontrol: meta.input_control
                return [ meta.input_control, meta.antibody, meta, bam, bai ]
            ips_wo_ipcontrol: !meta.input_control && !meta.is_input_control
                return [ meta.id, meta.antibody, meta, bam, bai ]
            ipcontrols: !meta.input_control && meta.is_input_control
                return [ meta.id, meta.input_control_of_antibody, meta, bam, bai ]
        }
        .set { ch_bam_by_type }

    // Create channel: [ meta, [ip_bams_merged_reps], [ip_bais_merged_reps], [ipcontrol_bams_merged_reps], [ipcontrol_bais_merged_reps] ]
    ch_bam_by_type
        .ips_with_ipcontrol
        .combine(ch_bam_by_type.ipcontrols, by: [0, 1])
        .map { ipcontrol_id, antibody, ip_meta, ip_bam, ip_bai, ipcontrol_meta, ipcontrol_bam, ipcontrol_bai ->
            [ ipcontrol_id, antibody, ip_meta, ip_bam, ip_bai, ipcontrol_bam, ipcontrol_bai ]
        }
        .mix(ch_bam_by_type.ips_wo_ipcontrol)
        // ips_wo_ipcontrol do not have ipcontrol_bam (it[5]) and ipcontrol_bai (it[6])
        .map { it ->
            def meta_clone = it[2].clone()
            meta_clone.id = meta_clone.id - ~/_bRep_.*$/
            meta_clone.input_control = meta_clone.input_control - ~/_bRep_.*$/
            [ meta_clone.id, it[1], meta_clone, it[3], it[4], it[5] ?: [], it[6] ?: [] ]
        }
        .groupTuple(by: [0, 1])
        .map {
            id, antibody, metas, ip_bams, ip_bais, ipcontrol_bams, ipcontrol_bais ->
                [ metas[0], ip_bams.flatten(), ip_bais.flatten(), ipcontrol_bams.flatten().unique(), ipcontrol_bais.flatten().unique() ]
        }
        .set { ch_ip_control_bam_bai_merged_reps }

    // TODO: Print to file for debuggin
    ch_ip_control_bam_bai_merged_reps
        .map {
            meta, ip_bams, ip_bais, ipcontrol_bams, ipcontrol_bais ->
                "${meta}\t${ip_bams}\t${ip_bais}\t${ipcontrol_bams}\t${ipcontrol_bais}"
        }
        .collectFile( name: 'ch_ip_control_bam_merged_reps.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_PEAKS_CALL_EPIC2" )


    //
    // Call peaks with epic2
    //
    EPIC2 (
        ch_ip_control_bam_bai_merged_reps,
        ch_chrom_sizes,
        efective_gfraction

    )

    //
    // Convert epic2 output to diffusePeak and broadPeak files
    //
    EPIC2_TO_PEAK (
        EPIC2.out.bed
    )

    //
    // Filter out samples with 0 epic2 peaks called
    //
    EPIC2_TO_PEAK
        .out
        .diffusepeak
        .filter {
            meta, peaks ->
                peaks.size() > 0
        }
        .set { ch_epic2_peaks }

    // Create channels: [ meta, ip_bam, peaks ]
    // Peaks are called on the pooled replicates, so join with the channel EPIC2 was given
    ch_ip_control_bam_bai_merged_reps
        .join(ch_epic2_peaks, by: 0)
        .map {
            meta, ip_bams, _ip_bais, _ipcontrol_bams, _ipcontrol_bais, peaks ->
                [ meta, ip_bams, peaks ]
        }
        // Split channel by ip_bam
        .transpose()
        .map {
            meta, ip_bam, peaks ->
                def meta_clone = meta.clone()
                meta_clone.id = ip_bam.getSimpleName()
                [ meta_clone, ip_bam, peaks ]
        }
        .set { ch_bam_peak }

    // Create channel: [ meta, broadpeak ], for the peak QC plots
    EPIC2_TO_PEAK
        .out
        .broadpeak
        .filter {
            meta, peaks ->
                peaks.size() > 0
        }
        .set { ch_epic2_broadpeaks }

    emit:
    peaks     = ch_epic2_peaks          // channel: [ val(meta), path(diffusepeak) ]
    bam_peaks = ch_bam_peak             // channel: [ val(meta), path(ip_bam), path(diffusepeak) ]
    qc_peaks  = ch_epic2_broadpeaks     // channel: [ val(meta), path(broadpeak) ]
    bed       = EPIC2.out.bed           // channel: [ val(meta), path(bed) ]
}
