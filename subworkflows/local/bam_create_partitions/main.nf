include { BAM_SPLIT_BY_STRAND                                       } from '../../../modules/local/bam_split_by_strand/main'
include { SAMTOOLS_INDEX                                            } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS                                        } from '../../../subworkflows/nf-core/bam_stats_samtools/main'
include { BEDTOOLS_GENOMECOV                                        } from '../../../modules/nf-core/bedtools/genomecov/main'
include { FILE_SORT as BEDGRAPH_SORT                                } from '../../../modules/local/file_sort/main'
include { BEDTOOLS_MAKEWINDOWS                                      } from '../../../modules/nf-core/bedtools/makewindows/main'
include { FILE_SORT as WINDOWS_SORT                                 } from '../../../modules/local/file_sort/main'
include { UCSC_BIGWIGAVERAGEOVERBED                                 } from '../../../modules/nf-core/ucsc/bigwigaverageoverbed/main'
include { FILE_SORT as BWAOB_SORT                                   } from '../../../modules/local/file_sort/main'
include { BEDGRAPH_NORMALIZE as BWAOB_NORMALIZE                     } from '../../../modules/local/bedgraph_normalize/main'
include { BEDGRAPH_SIGNAL_MINUS_INPUT                               } from '../../../modules/local/bedgraph_signal_minus_input/main'
include { PARTITION_OR_RFD_SMOOTH                                   } from '../../../modules/local/partition_or_rfd_smooth/main'
include { COLLECT_PARTITIONS                                        } from '../../../modules/local/collect_partitions/main'
include { PARTITION_AVERAGE                                         } from '../../../modules/local/partition_average/main'
include { UCSC_BEDGRAPHTOBIGWIG as UCSC_BEDGRAPHTOBIGWIG_WINDOWS    } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'
include { UCSC_BEDGRAPHTOBIGWIG as UCSC_BEDGRAPHTOBIGWIG_PARTITIONS } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'
include { RFD_TO_IZ                                                 } from '../../../modules/local/rfd_to_iz/main'
include { PARTITION_OR_RFD_PLOT as PARTITION_PLOT_SAMPLE            } from '../../../modules/local/partition_or_rfd_plot/main'
include { PARTITION_OR_RFD_PLOT as PARTITION_PLOT_GROUP             } from '../../../modules/local/partition_or_rfd_plot/main'
include { PARTITION_OR_RFD_PLOT as RFD_PLOT                         } from '../../../modules/local/partition_or_rfd_plot/main'


workflow BAM_CREATE_PARTITIONS {

    take:
    ch_bam                  // channel: [ val(meta), [ bam ] ]
    ch_fasta_fai            // channel: [ val(meta), path(fasta), path(fai) ]
    ch_chrom_sizes          // channel: [ bed ]
    ch_blacklist            // channel: [ val(meta), [ bed ] ]
    ch_okseq_rfd_file       // channel: [ val(meta), [ bed ] ]
    ch_initiation_zones     // channel: [ val(meta), [ bed ] ]
    smooth_radius
    derivative_radius
    zero_crossing_radius
    skip_partition_group_plot // boolean: skip the grouped partition plots

    main:

    // TODO: print for debugging
    ch_bam
        .map { meta, bam ->
            "${meta}\t${bam}"
        }
        .collectFile( name: '1_scar_ch_bam.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_CREATE_PARTITIONS")


    //
    // MODULE: Split BAMs by strand (forward and reverse)
    //
    BAM_SPLIT_BY_STRAND (
        ch_bam.map { meta, bam -> [ meta, bam, meta.strandedness ] }
    )

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

    //
    // MODULE: Run samtools stats, flagstat and idxstats per strand
    //
    BAM_STATS_SAMTOOLS (
        ch_bam.join(SAMTOOLS_INDEX.out.index, by: 0),
        ch_fasta_fai
    )

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
        ch_chrom_sizes.map { it -> it[1] },
        'bdg',
        true
    )

    //
    // MODULE: Sort the bedgraph so that it works with bedgraphtobigwig
    //
    BEDGRAPH_SORT (
        BEDTOOLS_GENOMECOV.out.genomecov,
        'bedGraph'
    )

    //
    // MODULE: Convert bedgraph to bigwig
    //
    UCSC_BEDGRAPHTOBIGWIG_WINDOWS (
        BEDGRAPH_SORT.out.sorted,
        ch_chrom_sizes.map { it -> it[1] }
    )
    ch_bigwig = UCSC_BEDGRAPHTOBIGWIG_WINDOWS.out.bigwig

    // TODO: print for debugging
    ch_bigwig
        .map { meta, bigwig ->
            "${meta}\t${bigwig}"
        }
        .collectFile( name: '6_scar_ch_bigwig.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_CREATE_PARTITIONS")

    //
    // MODULE: Create genomic windows
    //
    ch_windows = channel.empty()
    BEDTOOLS_MAKEWINDOWS (
        ch_chrom_sizes
    )
    ch_windows = BEDTOOLS_MAKEWINDOWS.out.bed

    //
    // MODULE: Sort windows
    //
    WINDOWS_SORT (
        ch_windows,
        'bed'
    )
    ch_windows = WINDOWS_SORT.out.sorted

    // count number of lines in the windows file
    ch_num_windows = ch_windows.map { meta, windows -> windows.countLines() }

    // TODO: print for debugging
    ch_num_windows
        .map { num_windows ->
            "${num_windows}"
        }
        .collectFile( name: '7_scar_ch_num_windows.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_CREATE_PARTITIONS")

    ch_bigwig
        .combine(ch_windows)
        .map { meta, bigwig, meta_windows, windows ->
            [ meta, windows, bigwig ]
        }
        .set { ch_windows_bigwig }

    //
    // MODULE: Calculate average coverage over windows
    //
    UCSC_BIGWIGAVERAGEOVERBED (
        ch_windows_bigwig
    )
    ch_bwaob = UCSC_BIGWIGAVERAGEOVERBED.out.tab

    // TODO: print for debugging
    ch_bwaob
        .map { meta, bwaob ->
            "${meta}\t${bwaob}"
        }
        .collectFile( name: '8_scar_ch_bwaob.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_CREATE_PARTITIONS")

    //
    // MODULE: Sort BWAOB
    //
    BWAOB_SORT (
        ch_bwaob,
        'tab'
    )
    ch_bwaob = BWAOB_SORT.out.sorted

    // RPM normalization factors
    // num_windows is used to add a pseudocount to the RPM normalization factor (prevent division by zero in partition_or_rfd_smooth)
    ch_bwaob
        .combine(ch_num_windows)
        .map { meta, bwaob, num_windows ->
            def meta_clone = meta.clone()
            meta_clone.norm_factor_val = 1e6 / (meta[meta.ref_total_mapped_reads_for_rpm_key] + (num_windows * 2))
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
    // MODULE: Normalize the 4th column of the bwaob file (sum of values over all bases covered)
    //
    BWAOB_NORMALIZE (
        ch_bwaob_rpm,
        'tab'
    )
    ch_norm = BWAOB_NORMALIZE.out.normalized

    // TODO: print for debugging
    ch_norm
        .map { meta, bdg ->
            "${meta}\t${bdg}"
        }
        .collectFile( name: '10_scar_ch_norm.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_CREATE_PARTITIONS")

    // for each of the strands, subtract the input from the sample
    ch_norm
        .branch { meta, bdg ->
            scar_with_ipcontrol: meta.input_control
                return [ meta.input_control, meta.antibody, meta.strand, meta, bdg ]
            ipcontrol: !meta.input_control && meta.is_input_control
                return [ meta.id, meta.input_control_of_antibody, meta.strand, bdg ]
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
            // 'partition' is for SCAR-seq and eSPAN, and 'RFD' for OK-seq
            def partition_or_rfd = meta.exp_type in ['SCAR-seq', 'eSPAN'] ? 'partition' : meta.exp_type == 'OK-seq' ? 'RFD' : null
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
            meta_clone.removeAll { it -> it.key in ['norm_factor_val', 'norm_factor_type', 'signal_minus_input'] }
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
        .set { ch_to_collect }

    // TODO: print for debugging
    ch_to_collect
        .map { meta, windows, bwaob_fwd, bwaob_rev, norm_or_smi_fwd, norm_or_smi_rev, rfd ->
            "${meta}\t${windows}\t${bwaob_fwd}\t${bwaob_rev}\t${norm_or_smi_fwd}\t${norm_or_smi_rev}\t${rfd}"
        }
        .collectFile( name: '16_scar_ch_partitions.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_CREATE_PARTITIONS")

    //
    // MODULE: Collect partitions
    //
    COLLECT_PARTITIONS (
        ch_to_collect
    )

    //
    // Add meta information to the filtered partitions and bedgraph
    //
    COLLECT_PARTITIONS.out.filtered_tsv
        .map { meta, tsv ->
            def meta_clone = meta.clone()
            meta_clone.flT_by_counts = true
            [ meta_clone, tsv ]
        }
        .set { ch_partitions_filtered }

    COLLECT_PARTITIONS.out.filtered_bdg
        .map { meta, bdg ->
            def meta_clone = meta.clone()
            meta_clone.flT_by_counts = true
            [ meta_clone, bdg ]
        }
        .set { ch_partitions_filtered_bdg }


    // Create channel: [ val(meta), partitions_brep ]
    ch_partitions_brep = channel.empty()
    COLLECT_PARTITIONS.out.tsv
        .map { meta, partition ->
            def meta_clone = meta.clone()
            def antibody = meta.antibody ?: meta.input_control_of_antibody
            meta_clone.id = meta.id - ~/_bRep_.*$/
            meta_clone.input_control = meta.input_control - ~/_bRep_.*$/
            [ meta_clone.id, antibody, meta.signal_minus_input, meta_clone, partition ]
        }
        .groupTuple(by: [0, 1, 2])
        // remove elements where there is only one biological replicate
        .filter { id, antibody, smi, metas, partitions ->
            partitions.size() > 1
        }
        .map { id, antibody, smi, metas, partitions ->
            // Sort metas and partitions to ensure consistent order
            def sorted_metas = metas.sort { meta -> meta.brep }
            def sorted_partitions = partitions.sort { partition -> partition.name }
            def meta_clone = sorted_metas[0].clone()
            meta_clone.averaged_brep = true
            [meta_clone, sorted_partitions]
        }
        .set { ch_partitions_brep }

    // TODO: print for debugging
    ch_partitions_brep
        .map { meta, partitions ->
            "${meta}\t${partitions}"
        }
        .collectFile( name: '18_scar_ch_partitions_brep.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_CREATE_PARTITIONS")

    //
    // MODULE: Create average partition across biological replicates
    //
    PARTITION_AVERAGE (
        ch_partitions_brep
    )

    //
    // Add meta information to the filtered partitions and bedgraph
    //
    PARTITION_AVERAGE.out.filtered_tsv
        .map { meta, tsv ->
            def meta_clone = meta.clone()
            meta_clone.flT_by_counts = true
            [ meta_clone, tsv ]
        }
        .set { ch_part_avg_filtered }

    PARTITION_AVERAGE.out.filtered_bdg
        .map { meta, bdg ->
            def meta_clone = meta.clone()
            meta_clone.flT_by_counts = true
            [ meta_clone, bdg ]
        }
        .set { ch_part_avg_filtered_bdg }

    //
    // MODULE: Convert the final partition bedgraph to bigwig
    //
    UCSC_BEDGRAPHTOBIGWIG_PARTITIONS (
        ch_partitions_filtered_bdg.mix(ch_part_avg_filtered_bdg),
        ch_chrom_sizes.map { it -> it[1] }
    )

    // Mix individual and averaged RFD (OK-seq) files
    ch_partitions_filtered
        .mix(ch_part_avg_filtered)
        .filter { it -> it[0].exp_type == 'OK-seq' }
        .set { ch_okseq }

    //
    // MODULE: Process OK-seq RFD file to get initiation zones
    //
    // Note: this IZ file is not used for plotting in the PARTITION_PLOT module,
    // it is only created here when an OK-seq sample is being processed
    // The one used for plotting is generated in PREPARE_GENOME
    RFD_TO_IZ (
        ch_okseq,
        ch_blacklist,
        ch_chrom_sizes
    )

    ch_okseq
        .combine(RFD_TO_IZ.out.iz_rm_overlaps_bed, by: 0)
        .set { ch_okseq_and_iz }


    //
    // MODULE: Plot RFD profiles around initiation zones
    //
    RFD_PLOT (
        ch_okseq_and_iz.map { meta, okseq, iz -> [ meta, [], [], [], okseq, iz ] },
        ch_blacklist,
        ch_chrom_sizes
    )


    // Create channel: [ val(meta), [ scar_tsv ], [ input_tsv ], [ minusinput_tsv ] ]
    ch_partitions_filtered
        .mix(ch_part_avg_filtered)
        .filter { it -> !it[0].exp_type == 'OK-seq' }
        .branch { meta, tsv ->
            scar_with_ipcontrol: !meta.is_input_control && !meta.signal_minus_input && meta.input_control
                return [ meta.input_control, meta, tsv ]
            scar_wo_ipcontrol: !meta.is_input_control && !meta.signal_minus_input && !meta.input_control
                return [ meta, tsv, [], [] ]
            ipcontrol: meta.is_input_control
                return [ meta.id, tsv ]
            minusipcontrol: meta.signal_minus_input
                return [ meta.id, tsv ]
        }
        .set { ch_part_flt_by_type }

    ch_part_flt_by_type
        .scar_with_ipcontrol
        .combine(ch_part_flt_by_type.ipcontrol, by: 0)
        .map { ipcontrol_id, meta_scar, scar_tsv, input_tsv ->
            [ meta_scar.id, meta_scar, scar_tsv, input_tsv ]
        }
        .combine(ch_part_flt_by_type.minusipcontrol, by: 0)
        .map { scar_id, meta_scar, scar_tsv, input_tsv, minusinput_tsv ->
            [ meta_scar, scar_tsv, input_tsv, minusinput_tsv ]
        }
        .mix(ch_part_flt_by_type.scar_wo_ipcontrol)
        .set { ch_part_flt }

    ch_part_flt
        .combine(ch_okseq_rfd_file.map { it -> it[1] })
        .combine(ch_initiation_zones.map { it -> it[1] })
        .set { ch_part_flt_to_plot }


    // TODO: print for debugging
    ch_part_flt_to_plot
        .map { meta_scar, scar_tsv, input_tsv, minusinput_tsv, okseq_rfd_file, initiation_zones ->
            "${meta_scar}\t${scar_tsv}\t${input_tsv}\t${minusinput_tsv}\t${okseq_rfd_file}\t${initiation_zones}"
        }
        .collectFile( name: '17_scar_ch_partitions_to_plot.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_CREATE_PARTITIONS")


    //
    // MODULE: Plot the final partition
    //
    PARTITION_PLOT_SAMPLE (
        ch_part_flt_to_plot,
        ch_blacklist,
        ch_chrom_sizes
    )

    if (!skip_partition_group_plot) {

        ch_part_flt
            .map { meta, scar_tsv, input_tsv, minusinput_tsv ->
                [ meta.antibody, meta.averaged_brep, meta, scar_tsv, input_tsv, minusinput_tsv ]
            }
            .groupTuple(by: [0, 1])
            // remove elements where there is only one biological replicate
            .filter { antibody, averaged_brep, metas, scar_tsvs, input_tsvs, minusinput_tsvs ->
                scar_tsvs.size() > 1
            }
            .map { antibody, averaged_brep, metas, scar_tsvs, input_tsvs, minusinput_tsvs ->
                // Sort metas and tsvs to ensure consistent order
                def sorted_metas = metas.sort { meta -> meta.id }
                // .flatten() here prevents the *nested* empty lists caused by saples wo_ipcontrol, which is not a valid path input value
                def sorted_scar_tsvs = scar_tsvs.flatten().sort { scar_tsv -> scar_tsv.name }
                def sorted_input_tsvs = input_tsvs.flatten().sort { input_tsv -> input_tsv.name }
                def sorted_minusinput_tsvs = minusinput_tsvs.flatten().sort { minusinput_tsv -> minusinput_tsv.name }
                def meta_clone = sorted_metas[0].clone()
                meta_clone.id = meta_clone.exp_type + '_' + antibody + (averaged_brep ? '_bRep_avg' : '')
                [ meta_clone, sorted_scar_tsvs, sorted_input_tsvs, sorted_minusinput_tsvs ]
            }
            .combine(ch_okseq_rfd_file.map { it -> it[1] })
            .combine(ch_initiation_zones.map { it -> it[1] })
            .set { ch_part_flt_group_to_plot }

        // TODO: print for debugging
        ch_part_flt_group_to_plot
            .map { meta, scar_tsvs, input_tsvs, minusinput_tsvs, okseq_rfd_file, initiation_zones ->
                "${meta}\t${scar_tsvs}\t${input_tsvs}\t${minusinput_tsvs}\t${okseq_rfd_file}\t${initiation_zones}"
            }
            .collectFile( name: '18_ch_part_flt_group_to_plot.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_CREATE_PARTITIONS")


        //
        // MODULE: Plot all the samples of the same antibody together
        //
        PARTITION_PLOT_GROUP (
            ch_part_flt_group_to_plot,
            ch_blacklist,
            ch_chrom_sizes
        )

    }

    emit:
    tab      = PARTITION_OR_RFD_SMOOTH.out.rfd       // channel: [ val(meta), [ tab ] ]
}
