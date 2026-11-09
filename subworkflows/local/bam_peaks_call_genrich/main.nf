//
// Call peaks with Genrich
//

include { SAMTOOLS_SORT } from '../../../modules/nf-core/samtools/sort/main'
include { GENRICH       } from '../../../modules/nf-core/genrich/main'

workflow BAM_PEAKS_CALL_GENRICH {
    take:
    ch_bam                            // channel: [ val(meta), [ ip_bam ], [ ipcontrol_bam ] ]
    ch_fasta_fai                      // channel: [ val(meta), path(fasta), path(fai) ]
    ch_blacklist                      // channel: [ val(meta), path(bed) ]

    main:

    SAMTOOLS_SORT ( ch_bam, ch_fasta_fai, '' )

    SAMTOOLS_SORT
        .out
        .bam
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
        .map { it ->
            def meta_clone = it[2].clone()
            meta_clone.id = meta_clone.id - ~/_bRep_.*$/
            meta_clone.input_control = meta_clone.input_control - ~/_bRep_.*$/
            [meta_clone.id, it[1], meta_clone, it[3], it[4] ?: []]
        }
        .groupTuple(by: [0, 1])
        .map { id, antibody, metas, ip_bams, ipcontrol_bams ->
            // Genrich matches controls to replicates by position, so a replicate without a control would shift the rest
            if (ipcontrol_bams.any { ctrl -> ctrl } && !ipcontrol_bams.every { ctrl -> ctrl }) {
                error("ERROR: Sample '${id}' has an input control for some of its replicates only. Genrich matches input controls to replicates by position, so either all replicates of a sample need an input control or none. Add the missing input controls, remove them all, or leave `genrich` out of `--peak_callers`.")
            }
            [metas[0], ip_bams.flatten(), ipcontrol_bams.flatten()]
        }
        .set { ch_ip_ipcontrol_bam_merged_reps }

    // TODO: Print to file for debuggin
    ch_ip_ipcontrol_bam_merged_reps
        .map { meta, ip_bams, ipcontrol_bams ->
            "${meta}\t${ip_bams}\t${ipcontrol_bams}"
        }
        .collectFile(name: 'ch_ip_ipcontrol_bam_merged_reps.txt', newLine: true, sort: false, storeDir: "${params.outdir}")


    //
    // Call peaks with Genrich
    //
    GENRICH (
        ch_ip_ipcontrol_bam_merged_reps,
        ch_blacklist.map { it -> it[1] }
    )

    //
    // Filter out samples with 0 Genrich peaks called
    //
    GENRICH.out.peak.filter { meta, peaks ->
        peaks.size() > 0
    }.set { ch_gr_peaks }

    // Create channels: [ meta, ip_bam, peaks ]
    // Peaks are called on the pooled replicates, so join with the channel GENRICH was given
    ch_ip_ipcontrol_bam_merged_reps
        .join(ch_gr_peaks, by: 0)
        .map { meta, ip_bams, _ipcontrol_bams, peaks ->
            [meta, ip_bams, peaks]
        }
        .transpose()
        .map { meta, ip_bam, peaks ->
            def meta_clone = meta.clone()
            meta_clone.id = ip_bam.getSimpleName()
            [meta_clone, ip_bam, peaks]
        }
        .set { ch_bam_peak }

    emit:
    peaks            = ch_gr_peaks                    // channel: [ val(meta), path(peaks) ]
    bam_peaks        = ch_bam_peak                    // channel: [ val(meta), path(ip_bam), path(peaks) ]
    qc_peaks         = ch_gr_peaks                    // channel: [ val(meta), path(peaks) ]
    bed_intervals    = GENRICH.out.bed_intervals      // channel: [ val(meta), path(bed) ]
    bedgraph_pileup  = GENRICH.out.bedgraph_pileup    // channel: [ val(meta), path(bedgraph) ]
    bedgraph_pvalues = GENRICH.out.bedgraph_pvalues   // channel: [ val(meta), path(bedgraph) ]
    duplicates       = GENRICH.out.duplicates         // channel: [ val(meta), path(bed) ]
}
