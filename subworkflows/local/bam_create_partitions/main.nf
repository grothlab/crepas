include { BAM_SPLIT_BY_STRAND                                       } from '../../../modules/local/bam_split_by_strand/main'
include { SAMTOOLS_INDEX                                            } from '../../../modules/nf-core/samtools/index/main'
include { BEDTOOLS_GENOMECOV                                        } from '../../../modules/nf-core/bedtools/genomecov/main'
include { FILE_SORT as BEDGRAPH_SORT                                    } from '../../../modules/local/file_sort/main'
include { BEDTOOLS_MAKEWINDOWS                                      } from '../../../modules/nf-core/bedtools/makewindows/main'
include { BIGTOOLS_BIGWIGAVERAGEOVERBED                             } from '../../../modules/local/bigtools/bigwigaverageoverbed/main'
include { BEDGRAPH_NORMALIZE                                            } from '../../../modules/local/bedgraph_normalize/main'
include { BEDGRAPH_SIGNAL_MINUS_INPUT                                   } from '../../../modules/local/bedgraph_signal_minus_input/main'
include { PARTITION_OR_RFD_SMOOTH                                          } from '../../../modules/local/partition_or_rfd_smooth/main'
include { COLLECT_PARTITIONS                                          } from '../../../modules/local/collect_partitions/main'
include { UCSC_BEDGRAPHTOBIGWIG as UCSC_BEDGRAPHTOBIGWIG_WINDOWS } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'
include { UCSC_BEDGRAPHTOBIGWIG as UCSC_BEDGRAPHTOBIGWIG_PARTITIONS } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'
include { PARTITION_PLOT                                      } from '../../../modules/local/partition_plot/main'


workflow BAM_CREATE_PARTITIONS {

    take:
    ch_bam                  // channel: [ val(meta), [ bam ] ]
    ch_chrom_sizes          // channel: [ bed ]
    ch_blacklist            // channel: [ val(meta), [ bed ] ]
    ch_okseq_rfd_file       // channel: [ val(meta), [ bed ] ]
    ch_initiation_zones     // channel: [ val(meta), [ bed ] ]
    rpm_use_flT2_total      // string: comma-separated list of antibodies for which to use flT2_total_mapped_reads instead of flT3_total_mapped_reads for RPM normalization
    smooth_radius
    derivative_radius
    zero_crossing_radius

    main:

    ch_versions = Channel.empty()

    // TODO: print for debugging
    ch_bam
        .map { meta, bam ->
            "${meta}\t${bam}"
        }
        .collectFile( name: '1_scar_ch_bam.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_CREATE_PARTITIONS")


    //
    // MODULE: Split BAMs by strand (forward and reverse)
    //
    BAM_SPLIT_BY_STRAND ( ch_bam )
    ch_versions = ch_versions.mix(BAM_SPLIT_BY_STRAND.out.versions.first())

    // Add strand to the meta information
    BAM_SPLIT_BY_STRAND
        .out
        .f_bam
        .map {
            meta, f_bam ->
                def meta_clone = meta.clone()
                meta_clone.strand = 'forward'
                [ meta_clone, f_bam ]
        }
        .set { ch_f_bam }

    // TODO: print for debugging
    ch_f_bam
        .map { meta, f_bam ->
            "${meta}\t${f_bam}"
        }
        .collectFile( name: '2_scar_ch_f_bam.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_CREATE_PARTITIONS")

    BAM_SPLIT_BY_STRAND
        .out
        .r_bam
        .map {
            meta, r_bam ->
                def meta_clone = meta.clone()
                meta_clone.strand = 'reverse'
                [ meta_clone, r_bam ]
        }
        .set { ch_r_bam }
    
    // TODO: print for debugging
    ch_r_bam
        .map { meta, r_bam ->
            "${meta}\t${r_bam}"
        }
        .collectFile( name: '3_scar_ch_r_bam.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_CREATE_PARTITIONS")

    ch_bam = ch_f_bam.mix(ch_r_bam)

    // TODO: print for debugging
    ch_bam
        .map { meta, bam ->
            "${meta}\t${bam}"
        }
        .collectFile( name: '4_scar_ch_bam.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_CREATE_PARTITIONS")

    //
    // MODULE: Index BAM files per strand
    //
    SAMTOOLS_INDEX (
        ch_bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    // Creating channel: [ val(meta), [ bam ], [ scale ] ] 
    ch_bam
        .map {
            meta, bam ->
                [ meta, bam, 1 ]
        }
        .set { ch_bam_scale }

    // TODO: print for debugging
    ch_bam_scale
        .map { meta, bam, scale ->
            "${meta}\t${bam}\t${scale}"
        }
        .collectFile( name: '5_scar_ch_bam_scale.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_CREATE_PARTITIONS")

    //
    // MODULE: Calculate genome coverage
    //
    BEDTOOLS_GENOMECOV (
        ch_bam_scale,
        ch_chrom_sizes.map { it[1] },
        'bdg',
        true
    )
    ch_versions  = ch_versions.mix(BEDTOOLS_GENOMECOV.out.versions.first())

    //
    // MODULE: Sort the bedgraph so that it works with bedgraphtobigwig
    //
    BEDGRAPH_SORT (
        BEDTOOLS_GENOMECOV.out.genomecov,
        'bedgraph'
    )
    ch_versions = ch_versions.mix(BEDGRAPH_SORT.out.versions.first())

    //
    // MODULE: Convert bedgraph to bigwig
    //
    UCSC_BEDGRAPHTOBIGWIG_WINDOWS (
        BEDGRAPH_SORT.out.sorted,
        ch_chrom_sizes.map { it[1] }
    )
    ch_bigwig = UCSC_BEDGRAPHTOBIGWIG_WINDOWS.out.bigwig
    ch_versions = ch_versions.mix(UCSC_BEDGRAPHTOBIGWIG_WINDOWS.out.versions.first())

    // TODO: print for debugging
    ch_bigwig
        .map { meta, bigwig ->
            "${meta}\t${bigwig}"
        }
        .collectFile( name: '6_scar_ch_bigwig.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_CREATE_PARTITIONS")

    //
    // MODULE: Create genomic windows
    //
    ch_windows = Channel.empty()
    BEDTOOLS_MAKEWINDOWS (
        ch_chrom_sizes
    )
    ch_windows = BEDTOOLS_MAKEWINDOWS.out.bed
    // count number of lines in the windows file
    ch_num_windows = ch_windows.map { meta, windows -> windows.countLines() }
    ch_versions = ch_versions.mix(BEDTOOLS_MAKEWINDOWS.out.versions)

    // TODO: print for debugging
    ch_num_windows
        .map { num_windows ->
            "${num_windows}"
        }
        .collectFile( name: '7_scar_ch_num_windows.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_CREATE_PARTITIONS")

    //
    // MODULE: Calculate average coverage over windows
    //
    BIGTOOLS_BIGWIGAVERAGEOVERBED (
        ch_bigwig,
        ch_windows
    )
    ch_bwaob = BIGTOOLS_BIGWIGAVERAGEOVERBED.out.bed
    ch_versions = ch_versions.mix(BIGTOOLS_BIGWIGAVERAGEOVERBED.out.versions.first())

    // TODO: print for debugging
    ch_bwaob
        .map { meta, bwaob ->
            "${meta}\t${bwaob}"
        }
        .collectFile( name: '8_scar_ch_bwaob.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_CREATE_PARTITIONS")


    // RPM normalization factors
    // num_windows is used to add a pseudocount to the RPM normalization factor (prevent division by zero in partition_or_rfd_smooth)
    ch_bwaob
        .combine(ch_num_windows)
        .map { meta, bwaob, num_windows ->
            def meta_clone = meta.clone()
            // samples have meta.antibody, while input controls have meta.input_control_of_antibody
            def antibody_to_use = meta.antibody ?: meta.input_control_of_antibody
            // If the sample was downsampled before, we want to use dSp_total_mapped_reads
            if (meta_clone.dSp_total_mapped_reads) {
                meta_clone.norm_factor_val = 1e6 / meta_clone.dSp_total_mapped_reads
                meta_clone.norm_factor_val_used = 'dSp_total_mapped_reads'
            // if antibody_to_use is in the list of antibodies or there is no flT3, use flT2 or flT1, otherwise use flT3
            } else if (rpm_use_flT2_total && antibody_to_use in rpm_use_flT2_total.split(',').collect { it.trim() } || !meta_clone.flT3_total_mapped_reads) {
                if (meta_clone.flT2_total_mapped_reads) {
                    meta_clone.norm_factor_val = 1e6 / (meta_clone.flT2_total_mapped_reads + num_windows)
                    meta_clone.norm_factor_val_used = 'flT2_total_mapped_reads'
                } else {
                    // Samples without spike-in wouldn't have flT2_total_mapped_reads, so we use flT1_total_mapped_reads instead
                    meta_clone.norm_factor_val = 1e6 / (meta_clone.flT1_total_mapped_reads + num_windows)
                    meta_clone.norm_factor_val_used = 'flT1_total_mapped_reads'
                }
            } else {
                meta_clone.norm_factor_val = 1e6 / (meta_clone.flT3_total_mapped_reads + num_windows)
                meta_clone.norm_factor_val_used = 'flT3_total_mapped_reads'
            }
            meta_clone.norm_factor_type = 'rpm'
            [ meta_clone, bwaob ]
        }
        .set { ch_bwaob_rpm }

    // TODO: print for debugging
    ch_bwaob_rpm
        .map { meta, bwaob ->
            "${meta}\t${bwaob}"
        }
        .collectFile( name: '9_scar_ch_bwaob_rpm.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_CREATE_PARTITIONS")

    //
    // MODULE: Normalize strands
    //
    BEDGRAPH_NORMALIZE (
        ch_bwaob_rpm
    )
    ch_norm = BEDGRAPH_NORMALIZE.out.bedgraph
    ch_versions = ch_versions.mix(BEDGRAPH_NORMALIZE.out.versions.first())

    // TODO: print for debugging
    ch_norm
        .map { meta, bdg ->
            "${meta}\t${bdg}"
        }
        .collectFile( name: '10_scar_ch_norm.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_CREATE_PARTITIONS")

    // for each of the strands, subtract the input from the sample
    ch_norm
        .map { meta, bdg ->
            // samples can have meta.antibody, while controls can have meta.input_control_of_antibody (if downsampling was performed)
            def antibody_to_use = meta.antibody ?: meta.input_control_of_antibody
            [ meta, antibody_to_use, bdg ]
        }
        .branch { meta, antibody, bdg ->
            scar_with_ipcontrol: meta.input_control
                return [ meta.input_control, antibody, meta.strand, meta, bdg ]
            ipcontrol: !meta.input_control && meta.is_input_control
                return [ meta.id, antibody, meta.strand, bdg ]
        }
        .set { ch_norm_by_type }

    // create channel: [ val(meta), [ scar_bdg ], [ input_bdg ] ]
    ch_norm_by_type
        .scar_with_ipcontrol
        .combine(ch_norm_by_type.ipcontrol, by: [0, 1, 2]) // combine by id, antibody, and strand
        .map { ipcontrol_id, antibody, strand, scar_meta, scar_bdg, input_bdg ->
            def meta_clone = scar_meta.clone()
                meta_clone.signal_minus_input = true
                [ meta_clone, scar_bdg, input_bdg ]
        }
        .set { ch_norm_scar_input }


    // TODO: print for debugging
    ch_norm_scar_input
        .map { meta, scar_bdg, input_bdg ->
            "${meta}\t${scar_bdg}\t${input_bdg}"
        }
        .collectFile( name: '11_scar_ch_norm_scar_input.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_CREATE_PARTITIONS")


    //
    // MODULE: Substract input from sample
    //
    BEDGRAPH_SIGNAL_MINUS_INPUT (
        ch_norm_scar_input
    )
    ch_bdg_smi = BEDGRAPH_SIGNAL_MINUS_INPUT.out.bedgraph
    ch_versions = ch_versions.mix(BEDGRAPH_SIGNAL_MINUS_INPUT.out.versions.first())

    // TODO: print for debugging
    ch_bdg_smi
        .map { meta, bdg ->
            "${meta}\t${bdg}"
        }
        .collectFile( name: '12_scar_ch_bdg_smi.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_CREATE_PARTITIONS")

    // create channel: [ val(meta), [ bdg_fwd ], [ bdg_rev ] ]
    ch_norm
        .mix(ch_bdg_smi)
        // copy meta and remove meta.strand to then merge fwd and rev by meta
        .map {
            meta, bdg ->
                def meta_clone = meta.clone()
                meta_clone.remove('strand')
                [ meta_clone, meta, bdg ]
        }
        .branch { meta_clone, meta, bdg ->
            forward: meta.strand == 'forward'
            reverse: meta.strand == 'reverse'
        }
        .set { ch_norm_and_smi }

    ch_norm_and_smi.forward
        .combine(ch_norm_and_smi.reverse, by: 0)
        .map { meta_clone, meta_fwd, bdg_fwd, meta_rev, bdg_rev ->
            [ meta_clone, bdg_fwd, bdg_rev ]
        }
        .set { ch_norm_and_smi }

    // TODO: print for debugging
    ch_norm_and_smi
        .map { meta, bdg_fwd, bdg_rev ->
            "${meta}\t${bdg_fwd}\t${bdg_rev}"
        }
        .collectFile( name: '13_scar_ch_norm_and_smi.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_CREATE_PARTITIONS")


    // Create channel: [ val(meta), val(partition_or_rfd), [ f_tab ], [ r_tab ] ]
    ch_norm_and_smi
        .map { meta, bdg_fwd, bdg_rev ->
            // 'partition' is for SCAR-seq and 'RFD' for OK-seq
            def partition_or_rfd = meta.exp_type == 'SCAR-seq' ? 'partition' : meta.exp_type == 'OK-seq' ? 'RFD' : null
            [ meta, partition_or_rfd, bdg_fwd, bdg_rev ]
        }
        .set { ch_part_norm_and_smi }

    //
    // MODULE: Calculate partitions (RFD)
    //
    PARTITION_OR_RFD_SMOOTH (
        ch_part_norm_and_smi,
        smooth_radius,
        derivative_radius,
        zero_crossing_radius
    )
    ch_rfd = PARTITION_OR_RFD_SMOOTH.out.rfd
    ch_versions = ch_versions.mix(PARTITION_OR_RFD_SMOOTH.out.versions.first())

    // TODO: print for debugging
    ch_rfd
        .map { meta, rfd ->
            "${meta}\t${rfd}"
        }
        .collectFile( name: '14_scar_ch_rfd.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_CREATE_PARTITIONS")


    // Prepare bwaob channel for combine()
    ch_bwaob
        .map { meta, bwaob ->
            def meta_clone = meta.clone()
            meta_clone.remove('strand')
            [ meta.strand, meta_clone, bwaob ]
        }
        .branch { strand, meta, bwaob ->
            forward: strand == 'forward'
                return [ meta, bwaob ]
            reverse: strand == 'reverse'
                return [ meta, bwaob ]
        }
        .set { ch_bwaob_strands }

    // Prepare norm and smi channel for combine()
    ch_norm_and_smi
        .map { meta, norm_or_smi_fwd, norm_or_smi_rev ->
            def meta_clone = meta.clone()
            meta_clone.removeAll { it.key in ['norm_factor_val', 'norm_factor_val_used', 'norm_factor_type', 'signal_minus_input'] }
            [ meta_clone, meta, norm_or_smi_fwd, norm_or_smi_rev ]
        }
        .set { ch_norm_and_smi_to_combine }


    // TODO: print for debugging
    ch_norm_and_smi_to_combine
        .map { meta, meta_norm_or_smi, norm_or_smi_fwd, norm_or_smi_rev ->
            "${meta}\t${meta_norm_or_smi}\t${norm_or_smi_fwd}\t${norm_or_smi_rev}"
        }
        .collectFile( name: '15_scar_ch_norm_and_smi_to_combine.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_CREATE_PARTITIONS")

    // Create channel: [ meta, windows, bwaob_fwd, bwaob_rev, norm_or_smi_fwd, norm_or_smi_rev, rfd ]
    ch_bwaob_strands.forward
        .combine(ch_bwaob_strands.reverse, by: 0) // this creates channel: [ meta, bwaob_fwd, bwaob_rev ]
        .combine(ch_norm_and_smi_to_combine, by: 0) // this creates channel: [ meta, bwaob_fwd, bwaob_rev, meta_norm_or_smi, norm_or_smi_fwd, norm_or_smi_rev ]
        .map { meta, bwaob_fwd, bwaob_rev, meta_norm_or_smi, norm_or_smi_fwd, norm_or_smi_rev ->
            [ meta_norm_or_smi, bwaob_fwd, bwaob_rev, norm_or_smi_fwd, norm_or_smi_rev ]
        }
        .combine(ch_rfd, by: 0) // this creates channel: [ meta, bwaob_fwd, bwaob_rev, norm_or_smi_fwd, norm_or_smi_rev, rfd ]
        .combine(ch_windows) // this creates channel: [ meta, bwaob_fwd, bwaob_rev, norm_or_smi_fwd, norm_or_smi_rev, rfd, windows ]
        .map { meta, bwaob_fwd, bwaob_rev, norm_or_smi_fwd, norm_or_smi_rev, rfd, meta_windows, windows ->
            [ meta, windows, bwaob_fwd, bwaob_rev, norm_or_smi_fwd, norm_or_smi_rev, rfd ]
        }
        .set { ch_partitions }

    // TODO: print for debugging
    ch_partitions
        .map { meta, windows, bwaob_fwd, bwaob_rev, norm_or_smi_fwd, norm_or_smi_rev, rfd ->
            "${meta}\t${windows}\t${bwaob_fwd}\t${bwaob_rev}\t${norm_or_smi_fwd}\t${norm_or_smi_rev}\t${rfd}"
        }
        .collectFile( name: '16_scar_ch_partitions.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_CREATE_PARTITIONS")

    //
    // MODULE: Collect partitions
    //
    COLLECT_PARTITIONS (
        ch_partitions
    )
    ch_versions = ch_versions.mix(COLLECT_PARTITIONS.out.versions.first())

    COLLECT_PARTITIONS
        .out
        .tsv
        .branch { meta, tsv ->
            scar_with_ipcontrol: !meta.is_input_control && !meta.signal_minus_input
                return [ meta.input_control, meta, tsv ]
            ipcontrol: meta.is_input_control
                return [ meta.id, tsv ]
            minusipcontrol: meta.signal_minus_input
                return [ meta.id, tsv ]
        }
        .set { ch_partitions_by_type }

    ch_partitions_by_type
        .scar_with_ipcontrol
        .combine(ch_partitions_by_type.ipcontrol, by: 0)
        .map { ipcontrol_id, meta_scar, scar_tsv, input_tsv ->
            [ meta_scar.id, meta_scar, scar_tsv, input_tsv ]
        }
        .combine(ch_partitions_by_type.minusipcontrol, by: 0)
        .map { scar_id, meta_scar, scar_tsv, input_tsv, minusinput_tsv ->
            [ meta_scar, scar_tsv, input_tsv, minusinput_tsv ]
        }
        .set { ch_partitions_to_plot }


    // TODO: print for debugging
    ch_partitions_to_plot
        .map { meta_scar, scar_tsv, input_tsv, minusinput_tsv ->
            "${meta_scar}\t${scar_tsv}\t${input_tsv}\t${minusinput_tsv}"
        }
        .collectFile( name: '17_scar_ch_partitions_to_plot.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_CREATE_PARTITIONS")


    //
    // MODULE: Plot the final partition
    //
    PARTITION_PLOT (
        ch_partitions_to_plot,
        ch_blacklist,
        ch_okseq_rfd_file,
        ch_initiation_zones,
        ch_chrom_sizes
    )
    ch_versions = ch_versions.mix(PARTITION_PLOT.out.versions.first())

    //
    // MODULE: Convert the final partition bedgraph to bigwig
    //
    UCSC_BEDGRAPHTOBIGWIG_PARTITIONS (
        COLLECT_PARTITIONS.out.bdg,
        ch_chrom_sizes.map { it[1] }
    )
    ch_versions = ch_versions.mix(UCSC_BEDGRAPHTOBIGWIG_PARTITIONS.out.versions.first())

    emit:
    tab      = PARTITION_OR_RFD_SMOOTH.out.rfd       // channel: [ val(meta), [ tab ] ]
    versions = ch_versions                    // channel: [ versions.yml ]
}

