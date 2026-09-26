include { DEEPTOOLS_BAMCOVERAGE                                         } from '../../../modules/nf-core/deeptools/bamcoverage/main'
include { BEDTOOLS_MAKEWINDOWS as BEDTOOLS_MAKEWINDOWS_ENDO             } from '../../../modules/nf-core/bedtools/makewindows/main'
include { BEDTOOLS_MAKEWINDOWS as BEDTOOLS_MAKEWINDOWS_EXO              } from '../../../modules/nf-core/bedtools/makewindows/main'
include { BEDTOOLS_MAP as BEDTOOLS_MAP_ENDO                             } from '../../../modules/nf-core/bedtools/map/main'
include { BEDTOOLS_MAP as BEDTOOLS_MAP_EXO                              } from '../../../modules/nf-core/bedtools/map/main'
include { BEDGRAPH_NORMALIZE                                            } from '../../../modules/local/bedgraph_normalize/main'
include { FILE_SORT as BEDGRAPH_SORT                                    } from '../../../modules/local/file_sort/main'
include { BEDGRAPH_SIGNAL_OVER_INPUT                                    } from '../../../modules/local/bedgraph_signal_over_input/main'
include { UCSC_BEDGRAPHTOBIGWIG as UCSC_BEDGRAPHTOBIGWIG_ENDO           } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'
include { UCSC_BEDGRAPHTOBIGWIG as UCSC_BEDGRAPHTOBIGWIG_EXO            } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'
include { DEEPTOOLS_BIGWIGCOMPARE                                       } from '../../../modules/nf-core/deeptools/bigwigcompare/main'
include { DEEPTOOLS_BIGWIGAVERAGE                                       } from '../../../modules/local/deeptools/bigwigaverage/main'
include { DEEPTOOLS_MULTIBIGWIGSUMMARY                                  } from '../../../modules/nf-core/deeptools/multibigwigsummary/main'
include { DEEPTOOLS_PLOTCORRELATION                                     } from '../../../modules/nf-core/deeptools/plotcorrelation/main'
include { DEEPTOOLS_PLOTPCA                                             } from '../../../modules/nf-core/deeptools/plotpca/main'


workflow BAM_NORMALIZE_BIGWIG_DEEPTOOLS {

    take:
    ch_bam_bai              // channel: [ val(meta), [ bam ], [ bai ] ]
    ch_chrom_sizes_endo     // channel: [ bed ]
    ch_chrom_sizes_exo      // channel: [ bed ]
    genome                  // string: genome name
    spikein_genome          // string: spike-in genome name
    min_reads_for_norm
    skip_srpm               // boolean: skip the SRPM normalization step
    skip_cisrpm             // boolean: skip the CISRPM normalization step
    skip_signal_vs_input          // boolean: skip the CISRPM-SOI normalization step
    signal_vs_input_operation
    skip_bw_average
    skip_exo_bw             // boolean: skip generating bigwigs for the exogenous genome
    skip_multibigwigsummary
    skip_plotcorrelation
    skip_plotpca

    main:

    ch_multiqc_files = channel.empty()

    ch_bam_bai
        .map { meta, bam, bai ->
            def meta_clone = meta.clone()
            meta_clone.norm_factor_val = 1
            meta_clone.norm_factor_type = 'raw'
            [ meta_clone, bam, bai ]
        }
        .set { ch_bam_bai }


    if (spikein_genome) {
        // Copy exogenous total_mapped_reads meta fields to their corresponding endogenous samples
        ch_bam_bai
            .map { meta, bam, bai ->
                // samples have meta.antibody, while input controls have meta.input_control_of_antibody
                def antibody = meta.antibody ?: meta.input_control_of_antibody
                [ meta.id, antibody, meta, bam, bai ]
            }
            .branch { id, antibody, meta, bam, bai ->
                endo: meta.genome == genome
                exo: meta.genome == spikein_genome
            }
            .set { ch_bam_bai_genome }

            ch_bam_bai_genome
                .endo
                .combine(ch_bam_bai_genome.exo, by: [0,1])
                .map { id, antibody, endo_meta, endo_bam, endo_bai, exo_meta, exo_bam, exo_bai ->
                    def meta_clone = endo_meta.clone()
                    meta_clone.exo_flT1_total_mapped_reads = exo_meta.flT1_total_mapped_reads
                    meta_clone.exo_flT2_total_mapped_reads = exo_meta.flT2_total_mapped_reads ?: null
                    meta_clone.exo_flT3_total_mapped_reads = exo_meta.flT3_total_mapped_reads ?: null
                    meta_clone.exo_flTbl_total_mapped_reads = exo_meta.flTbl_total_mapped_reads ?: null
                    meta_clone.exo_dSp_total_mapped_reads = exo_meta.dSp_total_mapped_reads ?: null
                    meta_clone.exo_ref_total_mapped_reads_key = "exo_" + exo_meta.last_total_mapped_reads_key
                    meta_clone.exo_ref_total_mapped_reads_for_dSp_key = "exo_" + exo_meta.ref_total_mapped_reads_for_dSp_key
                    meta_clone.exo_ref_total_mapped_reads_for_rpm_key = "exo_" + exo_meta.ref_total_mapped_reads_for_rpm_key
                    meta_clone.exo_ref_total_mapped_reads_for_srpm_key = "exo_" + exo_meta.ref_total_mapped_reads_for_srpm_key
                    meta_clone.exo_ref_total_mapped_reads_for_cisrpm_key = "exo_" + exo_meta.ref_total_mapped_reads_for_cisrpm_key
                    [ meta_clone, endo_bam, endo_bai ]
                }
                .set { ch_bam_bai_endo }

            if (!skip_exo_bw) {
                // Mix back the exogenous BAMs to process them together with the endogenous ones
                ch_bam_bai = ch_bam_bai_endo.mix(ch_bam_bai_genome.exo)
            } else {
                ch_bam_bai = ch_bam_bai_endo
            }
    }

    // Remove empty BAMs to prevent bamCoverage errors
    ch_bam_bai =  ch_bam_bai.filter { meta, bam, bai -> meta[meta.last_total_mapped_reads_key] >= min_reads_for_norm }

    //
    // MODULE: Calculate raw coverage per bin
    //
    DEEPTOOLS_BAMCOVERAGE (
        ch_bam_bai,
        [],
        [],
        channel.value([[:], []])
    )
    ch_bdg_raw = DEEPTOOLS_BAMCOVERAGE.out.bedgraph

    // TODO: print for debugging
    ch_bdg_raw
        .map { meta, bdg ->
            "${meta}\t${bdg}"
        }
        .collectFile( name: 'ch_bdg_raw.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_NORMALIZE_BIGWIG_DEEPTOOLS" )

    // bamCoverage merges contiguous bins with the same coverage and there is no option to disable this.
    // See https://github.com/deeptools/deepTools/issues/907#issuecomment-576729674
    // Therefore, to get bins of equal size across the whole genome,
    // we have to make windows and map the bedgraph to the windows.
    // An alternative would be to use featureCounts:

    //
    // MODULE: Make windows of equal size across the whole genome
    //
    BEDTOOLS_MAKEWINDOWS_ENDO (
        ch_chrom_sizes_endo
    )
    ch_windows_endo = BEDTOOLS_MAKEWINDOWS_ENDO.out.bed

    // Create channel: [ val(meta_bdg_raw), windows, bdg_raw ]
    ch_bdg_raw
        .filter { meta, bdg ->
            meta.genome == genome
        }
        .combine(ch_windows_endo)
        .map { meta_bdg_raw, bdg_raw, meta_windows, windows ->
            [ meta_bdg_raw, windows, bdg_raw ]
        }
        .set { ch_windows_endo_bdg_raw }

    //
    // MODULE: Map the coverage bedgraph to the windows (endogenous genome)
    //
    BEDTOOLS_MAP_ENDO (
        ch_windows_endo_bdg_raw,
        ch_chrom_sizes_endo
    )
    ch_bdg_map_endo = BEDTOOLS_MAP_ENDO.out.mapped

    ch_bdg_map = ch_bdg_map_endo
    ch_windows_exo = channel.empty()
    ch_windows_exo_bdg_raw = channel.empty()
    if (spikein_genome && !skip_exo_bw) {
        BEDTOOLS_MAKEWINDOWS_EXO (
            ch_chrom_sizes_exo
        )
        ch_windows_exo = BEDTOOLS_MAKEWINDOWS_EXO.out.bed

        // Create channel: [ val(meta_bdg_raw), windows, bdg_raw ]
        ch_bdg_raw
            .filter { meta, bdg ->
                meta.genome == spikein_genome
            }
            .combine(ch_windows_exo)
            .map { meta_bdg_raw, bdg_raw, meta_windows, windows ->
                [ meta_bdg_raw, windows, bdg_raw ]
            }
            .set { ch_windows_exo_bdg_raw }

        //
        // MODULE: Map the coverage bedgraph to the windows (spike-in genome)
        //
        BEDTOOLS_MAP_EXO (
            ch_windows_exo_bdg_raw,
            ch_chrom_sizes_exo
        )
        ch_bdg_map_exo = BEDTOOLS_MAP_EXO.out.mapped

        // Merge the two channels
        ch_bdg_map = ch_bdg_map_endo.mix(ch_bdg_map_exo)
    }

    // TODO: print for debugging
    ch_bdg_map
        .map { meta, bdg ->
            "${meta}\t${bdg}"
        }
        .collectFile( name: 'ch_bdg_map.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_NORMALIZE_BIGWIG_DEEPTOOLS" )

    // RPM normalization factors
    ch_bdg_rpm = channel.empty()
    ch_bdg_map
        .filter { meta, bdg ->
            def denom = meta[meta.ref_total_mapped_reads_for_rpm_key]
            if (!denom) {
                log.warn "Skipping RPM normalization for sample '${meta.id}': total mapped reads (key '${meta.ref_total_mapped_reads_for_rpm_key}') is missing or zero."
                return false
            }
            true
        }
        .map { meta, bdg ->
            def meta_clone = meta.clone()
            meta_clone.norm_factor_val = 1e6 / meta[meta.ref_total_mapped_reads_for_rpm_key]
            meta_clone.norm_factor_type = 'rpm'
            [ meta_clone, bdg ]
        }
        .set { ch_bdg_rpm }

    // TODO: print for debugging
    ch_bdg_rpm
        .map {
            meta, bdg ->
                "${meta}\t${bdg}"
        }
        .collectFile( name: 'ch_bdg_rpm.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_NORMALIZE_BIGWIG_DEEPTOOLS" )

    // Copy and modify channel meta to add SRPM normalization factors (for ChIPs)
    ch_bdg_srpm = channel.empty()
    if (spikein_genome && !skip_srpm) {
        ch_bdg_map
            .map { meta, bdg ->
                // samples have meta.antibody, while input controls have meta.input_control_of_antibody
                def antibody_to_use = meta.antibody ?: meta.input_control_of_antibody
                [ meta.id, antibody_to_use, meta, bdg ]
            }
            .branch { id, antibody, meta, bdg ->
                endo: meta.genome == genome
                exo: meta.genome == spikein_genome
            }
            .set { ch_bdg_genome }

        // TODO: print for debugging
        ch_bdg_genome.endo
            .map { it ->
                "${it}"
            }
            .collectFile( name: 'ch_bdg_srpm_genome_endo.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_NORMALIZE_BIGWIG_DEEPTOOLS" )

        // TODO: print for debugging
        ch_bdg_genome.exo
            .map { it ->
                "${it}"
            }
            .collectFile( name: 'ch_bdg_srpm_genome_exo.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_NORMALIZE_BIGWIG_DEEPTOOLS" )

        // Combine the endo and exo BAMs if we want to produce exogenous bigWigs
        if (!skip_exo_bw) {
            ch_bdg_genome
                .endo
                .combine(ch_bdg_genome.exo, by: [0,1])
                .filter { id, antibody, endo_meta, endo_bdg, exo_meta, exo_bdg ->
                    def denom = exo_meta[exo_meta.ref_total_mapped_reads_for_srpm_key]
                    if (!denom) {
                        log.warn "Skipping SRPM normalization for sample '${endo_meta.id}': spike-in total mapped reads (key '${exo_meta.ref_total_mapped_reads_for_srpm_key}') is missing or zero."
                        return false
                    }
                    true
                }
                .map { id, antibody, endo_meta, endo_bdg, exo_meta, exo_bdg ->
                    def meta_clone = endo_meta.clone()
                    meta_clone.norm_factor_val = 1e6 / exo_meta[exo_meta.ref_total_mapped_reads_for_srpm_key]
                    meta_clone.norm_factor_type = 'srpm'
                    [ meta_clone, endo_bdg ]
                }
                .set { ch_bdg_srpm }
        // Use the exogenous totals saved in the the endogenous metadata otherwise
        } else {
            ch_bdg_genome
                .endo
                .filter { id, antibody, endo_meta, endo_bdg ->
                    def denom = endo_meta[endo_meta.exo_ref_total_mapped_reads_for_srpm_key]
                    if (!denom) {
                        log.warn "Skipping SRPM normalization for sample '${endo_meta.id}': spike-in total mapped reads (key '${endo_meta.exo_ref_total_mapped_reads_for_srpm_key}') is missing or zero."
                        return false
                    }
                    true
                }
                .map { id, antibody, endo_meta, endo_bdg ->
                    def meta_clone = endo_meta.clone()
                    meta_clone.norm_factor_val = 1e6 / endo_meta[endo_meta.exo_ref_total_mapped_reads_for_srpm_key]
                    meta_clone.norm_factor_type = 'srpm'
                    [ meta_clone, endo_bdg ]
                }
                .set { ch_bdg_srpm }
        }

    }

    // TODO: print for debugging
    ch_bdg_srpm
        .map { it ->
            "${it}"
        }
        .collectFile( name: 'ch_bdg_srpm.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_NORMALIZE_BIGWIG_DEEPTOOLS" )


    // Copy and modify channel meta to add CISRPM normalization factors
    ch_bdg_genome_type = channel.empty()
    ch_bdg_genome_ip = channel.empty()
    ch_bdg_genome_ipcontrol = channel.empty()
    ch_bdg_ip_cisrpm = channel.empty()
    ch_bdg_ipcontrol_cisrpm = channel.empty()
    ch_bdg_cisrpm = channel.empty()
    // "if (spikein_genome)" is needed, otherwise cisrpm will be attempted for
    // controls, and this will fail, since there is no ref_total_mapped_reads_for_cisrpm_key
    if (spikein_genome && !skip_cisrpm) {
        // Split BAMs by genome (endo and exo) and by type (ip and ipcontrol)
        ch_bdg_map
            .map { meta, bdg ->
                // samples have meta.antibody, while input controls have meta.input_control_of_antibody
                def antibody_to_use = meta.antibody ?: meta.input_control_of_antibody
                [ meta.id, antibody_to_use, meta, bdg ]
            }
            .branch { id, antibody, meta, bdg ->
                endo_ip: meta.genome == genome && !meta.is_input_control
                endo_ipcontrol: meta.genome == genome && meta.is_input_control
                exo_ip: meta.genome == spikein_genome && !meta.is_input_control
                exo_ipcontrol: meta.genome == spikein_genome && meta.is_input_control
            }
            .set { ch_bdg_genome_type }


        // Combine the endo and exo BAMs if we want to produce exogenous bigWigs
        if (!skip_exo_bw) {

            // Combine the endo and exo BAMs (ChIPs)
            ch_bdg_genome_type.endo_ip
                .combine(ch_bdg_genome_type.exo_ip, by: [0,1])
                .map { ip_id, ip_antibody, endo_ip_meta, endo_ip_bdg, exo_ip_meta, exo_ip_bdg ->
                        [ endo_ip_meta.input_control, endo_ip_meta.antibody, endo_ip_meta, endo_ip_bdg, exo_ip_meta, exo_ip_bdg ]
                }
                .set { ch_bdg_genome_ip }

            // Combine the endo and exo BAMs (inputs)
            ch_bdg_genome_type.endo_ipcontrol
                .combine(ch_bdg_genome_type.exo_ipcontrol, by: [0,1])
                .map { ipcontrol_id, ipcontrol_antibody, endo_ipcontrol_meta, endo_ipcontrol_bdg, exo_ipcontrol_meta, exo_ipcontrol_bdg ->
                    [ endo_ipcontrol_meta.id, endo_ipcontrol_meta.input_control_of_antibody, endo_ipcontrol_meta, endo_ipcontrol_bdg, exo_ipcontrol_meta, exo_ipcontrol_bdg ]
                }
                .set { ch_bdg_genome_ipcontrol }

            // Combine the combined ChIPs with the combined inputs
            ch_bdg_genome_ip
                .combine(ch_bdg_genome_ipcontrol, by: [0,1])
                .filter { id, antibody, endo_ip_meta, endo_ip_bdg, exo_ip_meta, exo_ip_bdg, endo_ipcontrol_meta, endo_ipcontrol_bdg, exo_ipcontrol_meta, exo_ipcontrol_bdg ->
                    def exo_ip_denom = exo_ip_meta[exo_ip_meta.ref_total_mapped_reads_for_cisrpm_key]
                    def endo_ipcontrol_denom = endo_ipcontrol_meta[endo_ipcontrol_meta.ref_total_mapped_reads_for_cisrpm_key]
                    def exo_ipcontrol_denom = exo_ipcontrol_meta[exo_ipcontrol_meta.ref_total_mapped_reads_for_cisrpm_key]
                    if (!exo_ip_denom || !endo_ipcontrol_denom || !exo_ipcontrol_denom) {
                        log.warn "Skipping CISRPM normalization for sample '${endo_ip_meta.id}': spike-in or input control total mapped reads is missing or zero."
                        return false
                    }
                    true
                }
                .map { id, antibody, endo_ip_meta, endo_ip_bdg, exo_ip_meta, exo_ip_bdg, endo_ipcontrol_meta, endo_ipcontrol_bdg, exo_ipcontrol_meta, exo_ipcontrol_bdg ->
                        def meta_clone = endo_ip_meta.clone()
                        meta_clone.norm_factor_val = (1e6 / exo_ip_meta[exo_ip_meta.ref_total_mapped_reads_for_cisrpm_key]) * (exo_ipcontrol_meta[exo_ipcontrol_meta.ref_total_mapped_reads_for_cisrpm_key] / endo_ipcontrol_meta[endo_ipcontrol_meta.ref_total_mapped_reads_for_cisrpm_key])
                        meta_clone.norm_factor_type = 'cisrpm'
                        [ meta_clone, endo_ip_bdg ]
                }
                .set { ch_bdg_ip_cisrpm }

        // Use the exogenous totals saved in the the endogenous metadata otherwise
        } else {
            ch_bdg_genome_type
                .endo_ip
                .map { ip_id, ip_antibody, endo_ip_meta, endo_ip_bdg ->
                        [ endo_ip_meta.input_control, endo_ip_meta.antibody, endo_ip_meta, endo_ip_bdg ]
                }
                .set { ch_bdg_genome_ip }

            ch_bdg_genome_type
                .endo_ipcontrol
                .map { ipcontrol_id, ipcontrol_antibody, endo_ipcontrol_meta, endo_ipcontrol_bdg ->
                    [ endo_ipcontrol_meta.id, ipcontrol_antibody, endo_ipcontrol_meta, endo_ipcontrol_bdg ]
                }
                .set { ch_bdg_genome_ipcontrol }

            ch_bdg_genome_ip
                .combine(ch_bdg_genome_ipcontrol, by: [0,1])
                .filter { ipcontrol_id, ip_antibody, endo_ip_meta, endo_ip_bdg, endo_ipcontrol_meta, endo_ipcontrol_bdg ->
                    def exo_ip_denom = endo_ip_meta[endo_ip_meta.exo_ref_total_mapped_reads_for_cisrpm_key]
                    def endo_ipcontrol_denom = endo_ipcontrol_meta[endo_ipcontrol_meta.ref_total_mapped_reads_for_cisrpm_key]
                    def exo_ipcontrol_denom = endo_ipcontrol_meta[endo_ipcontrol_meta.exo_ref_total_mapped_reads_for_cisrpm_key]
                    if (!exo_ip_denom || !endo_ipcontrol_denom || !exo_ipcontrol_denom) {
                        log.warn "Skipping CISRPM normalization for sample '${endo_ip_meta.id}': spike-in or input control total mapped reads is missing or zero."
                        return false
                    }
                    true
                }
                .map { ipcontrol_id, ip_antibody, endo_ip_meta, endo_ip_bdg, endo_ipcontrol_meta, endo_ipcontrol_bdg ->
                        def meta_clone = endo_ip_meta.clone()
                        meta_clone.norm_factor_val = (1e6 / endo_ip_meta[endo_ip_meta.exo_ref_total_mapped_reads_for_cisrpm_key]) * (endo_ipcontrol_meta[endo_ipcontrol_meta.exo_ref_total_mapped_reads_for_cisrpm_key] / endo_ipcontrol_meta[endo_ipcontrol_meta.ref_total_mapped_reads_for_cisrpm_key])
                        meta_clone.norm_factor_type = 'cisrpm'
                        [ meta_clone, endo_ip_bdg ]
                }
                .set { ch_bdg_ip_cisrpm }
        }
    }

    // Now do the missing CISRPM for the endogenous inputs
    // In this case CISRPM is the same as RPM, but we cannot
    // just copy the RPM from before, since ref_total_mapped_reads_for_cisrpm_key
    // can be different than ref_total_mapped_reads_for_rpm_key
    ch_bdg_map
        .filter { meta, bdg ->
            meta.genome == genome && meta.is_input_control
        }
        .filter { meta, bdg ->
            def denom = meta[meta.ref_total_mapped_reads_for_cisrpm_key]
            if (!denom) {
                log.warn "Skipping CISRPM normalization for sample '${meta.id}': total mapped reads (key '${meta.ref_total_mapped_reads_for_cisrpm_key}') is missing or zero."
                return false
            }
            true
        }
        .map { meta, bdg ->
            def meta_clone = meta.clone()
            meta_clone.norm_factor_val = 1e6 / meta[meta.ref_total_mapped_reads_for_cisrpm_key]
            meta_clone.norm_factor_type = 'cisrpm'
            [ meta_clone, bdg ]
        }
        .set { ch_bdg_ipcontrol_cisrpm }

    ch_bdg_cisrpm = ch_bdg_ip_cisrpm.mix(ch_bdg_ipcontrol_cisrpm)

    // TODO: print for debugging
    ch_bdg_cisrpm
        .map {
            meta, bdg ->
                "${meta}\t${bdg}"
        }
        .collectFile( name: 'ch_bdg_cisrpm.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_NORMALIZE_BIGWIG_DEEPTOOLS" )

    // Mix all to-be-normalized bedgraphs in one channel
    ch_bdg_rpm
        .mix(ch_bdg_srpm)
        .mix(ch_bdg_cisrpm)
        .set { ch_bdg_norm }

    // TODO: print for debugging
    ch_bdg_norm
        .map {
            meta, bdg ->
                "${meta}\t${bdg}"
        }
        .collectFile( name: 'ch_bdg_norm.txt', newLine: true, sort: false, storeDir: "${params.outdir}/.debug/BAM_NORMALIZE_BIGWIG_DEEPTOOLS" )

    //
    // MODULE: Normalize the bedgraphs
    //
    BEDGRAPH_NORMALIZE (
        ch_bdg_norm,
        'bedGraph'
    )
    ch_bdg_map_norm = ch_bdg_map.mix(BEDGRAPH_NORMALIZE.out.normalized)

    //
    // MODULE: Sort the bedgraph so that it works with ucsc_bedgraphtobigwig
    //
    BEDGRAPH_SORT (
        ch_bdg_map_norm,
        'bedGraph'
    )
    ch_bdg_all = BEDGRAPH_SORT.out.sorted


    if (!skip_signal_vs_input && signal_vs_input_operation == 'soi') {
        ch_bdg_all
            .branch { meta, bdg ->
                ips_with_ipcontrol: meta.input_control
                    return [ meta.input_control, meta.antibody, meta.genome, meta.norm_factor_type, meta, bdg ]
                ipcontrols: !meta.input_control && meta.is_input_control
                    return [ meta.id, meta.input_control_of_antibody, meta.genome, meta.norm_factor_type, meta, bdg ]
            }
            .set { ch_bdg_per_type }

        ch_bdg_per_type
            .ips_with_ipcontrol
            .combine(ch_bdg_per_type.ipcontrols, by: [0, 1, 2, 3])
            .map { ipcontrol_id, ip_antibody, ip_genome, ip_norm_factor_type, ip_meta, ip_bdg, ipcontrol_meta, ipcontrol_bdg ->
                def meta_clone = ip_meta.clone()
                    meta_clone.signal_vs_input = true
                    meta_clone.signal_vs_input_operation = signal_vs_input_operation
                    [ meta_clone, ip_bdg, ipcontrol_bdg ]
            }
            .set { ch_bdg_ip_control_soi }

        //
        // MODULE: Calculate signal over input
        //
        BEDGRAPH_SIGNAL_OVER_INPUT (
            ch_bdg_ip_control_soi
        )
        ch_bdg_all = BEDGRAPH_SIGNAL_OVER_INPUT.out.bedgraph.mix(ch_bdg_all)

    }

    //
    // MODULE: Convert bedgraph to bigwig
    //
    UCSC_BEDGRAPHTOBIGWIG_ENDO (
        ch_bdg_all.filter { it -> it[0].genome == genome },
        ch_chrom_sizes_endo.map { it -> it[1] }
    )

    //
    // MODULE: Convert bedgraph to bigwig for the exogenous genome
    //
    UCSC_BEDGRAPHTOBIGWIG_EXO (
        ch_bdg_all.filter { it -> it[0].genome == spikein_genome },
        ch_chrom_sizes_exo.map { it -> it[1] }
    )

    // Mix the endogenous and exogenous bigwigs back together
    UCSC_BEDGRAPHTOBIGWIG_ENDO
        .out
        .bigwig
        .mix(UCSC_BEDGRAPHTOBIGWIG_EXO.out.bigwig)
        .set { ch_bigwig }


    ch_bw_compare = channel.empty()
    if (!skip_signal_vs_input && signal_vs_input_operation != 'soi') {

        ch_bigwig
            .branch { meta, bw ->
                ips_with_ipcontrol: meta.input_control
                    return [ meta.input_control, meta.antibody, meta.genome, meta.norm_factor_type, meta, bw ]
                ipcontrols: !meta.input_control && meta.is_input_control
                    return [ meta.id, meta.input_control_of_antibody, meta.genome, meta.norm_factor_type, meta, bw ]
            }
            .set { ch_bigwig_per_type }

        ch_bigwig_per_type
            .ips_with_ipcontrol
            .combine(ch_bigwig_per_type.ipcontrols, by: [0, 1, 2, 3])
            .map { ipcontrol_id, ip_antibody, ip_genome, ip_norm_factor_type, ip_meta, ip_bw, ipcontrol_meta, ipcontrol_bw ->
                def meta_clone = ip_meta.clone()
                    meta_clone.signal_vs_input = true
                    meta_clone.signal_vs_input_operation = signal_vs_input_operation
                    [ meta_clone, ip_bw, ipcontrol_bw ]
            }
            .set { ch_bigwig_ip_control_compare }

        //
        // MODULE: Compute ratio, log2ratio, etc., between IPs and their input controls
        //
        DEEPTOOLS_BIGWIGCOMPARE (
            ch_bigwig_ip_control_compare,
            channel.value([[:], []])
        )
        ch_bw_compare = DEEPTOOLS_BIGWIGCOMPARE.out.output

    }

    ch_bw_avg = channel.empty()
    if (!skip_bw_average) {

        // Create channel: [ val(meta), [ bRep_bigwigs ] ]
        ch_bigwig
            .mix(ch_bw_compare)
            .map { meta, bw ->
                def meta_clone = meta.clone()
                def antibody = meta.antibody ?: meta.input_control_of_antibody
                meta_clone.id = meta_clone.id - ~/_bRep_.*$/
                [ meta_clone.id, antibody, meta_clone.genome, meta_clone.norm_factor_type, meta_clone.signal_vs_input, meta_clone, bw ]
            }
            .groupTuple(by: [0, 1, 2, 3, 4])
            // remove elements where there is only one biological replicate
            .filter { id, antibody, meta_genome, norm_factor_type, signal_vs_input, metas, bws ->
                bws.size() > 1
            }
            .map { id, antibody, meta_genome, norm_factor_type, signal_vs_input, metas, bws ->
                // Sort metas and bws by id to ensure consistent order
                def sorted_metas = metas.sort { meta -> meta.brep }
                def sorted_bws = bws.sort { bw -> bw.name }
                def meta_clone = sorted_metas[0].clone()
                meta_clone.averaged_brep = true
                [ meta_clone, sorted_bws ]
            }
            .set { ch_bdg_all_brep_bw }

        //
        // MODULE: Average bigwigs across replicates
        //
        DEEPTOOLS_BIGWIGAVERAGE (
            ch_bdg_all_brep_bw,
            channel.value([[:], []])
        )
        ch_bw_avg = DEEPTOOLS_BIGWIGAVERAGE.out.bigwig
    }

    ch_bigwig
        .mix(ch_bw_compare)
        .mix(ch_bw_avg)
        .set { ch_bw_all }

    ch_bigwig_all_endo  = ch_bw_all.filter { it -> it[0].genome == genome }


    if (!skip_multibigwigsummary) {

       ch_bigwig_all_endo
            .map { meta, bw ->
                def antibody = meta.antibody ?: meta.input_control_of_antibody
                [ antibody, meta.exp_type, meta.norm_factor_type, meta.signal_vs_input_operation, meta.averaged_brep, meta.id, meta, bw ]
            }
            .groupTuple(by: [0, 1, 2, 3, 4])
            .map { antibody, exp_type, norm_factor_type, signal_vs_input_op, averaged_brep, ids, metas, bws ->
                // Sort ids, metas and bws by id to ensure consistent order in plots
                def sorted_ids = ids.sort()
                def sorted_metas = metas.sort { meta -> meta.id }
                def sorted_bws = bws.sort { it -> it.name }
                [ antibody, exp_type, norm_factor_type, signal_vs_input_op, averaged_brep, sorted_ids, sorted_metas, sorted_bws ]
            }
            .map {
                antibody, exp_type, norm_factor_type, signal_vs_input_op, averaged_brep, ids, metas, bws ->
                    def meta_new = metas[0].clone()
                    meta_new.id = exp_type + '_' +
                        (antibody ? antibody : 'no_antibody') +
                        '_' + norm_factor_type +
                        (signal_vs_input_op ? '_' + signal_vs_input_op : '') +
                        (averaged_brep ? '_' + 'bRep_avg' : '')
                    meta_new.antibody = antibody
                    meta_new.ids = ids
                    [ meta_new, bws.flatten(), ids.flatten() ]
            }
            // remove cases with only one sample
            .filter { meta, bws, ids -> bws.size() > 1 }
            .set { ch_bw_for_summary }

        //
        // MODULE: multiBigwigSummary
        //
        DEEPTOOLS_MULTIBIGWIGSUMMARY (
            ch_bw_for_summary,
            channel.value([[:], []])
        )

        if (!skip_plotcorrelation) {
            //
            // MODULE: plotCorrelation
            //
            DEEPTOOLS_PLOTCORRELATION (
                DEEPTOOLS_MULTIBIGWIGSUMMARY.out.matrix,
                "spearman",
                "heatmap"
            )
            ch_multiqc_files = ch_multiqc_files.mix(DEEPTOOLS_PLOTCORRELATION.out.matrix.collect { it -> it[1] })
        }

        if (!skip_plotpca) {
            //
            // MODULE: plotPCA
            //
            DEEPTOOLS_PLOTPCA (
                DEEPTOOLS_MULTIBIGWIGSUMMARY.out.matrix
            )
            ch_multiqc_files = ch_multiqc_files.mix(DEEPTOOLS_PLOTPCA.out.tab.collect { it -> it[1] })
        }

    }


    emit:

    multiqc_files    = ch_multiqc_files                         // channel: [ path(mqc_files) ]
    bedgraph_endo    = ch_bdg_all.filter { it -> it[0].genome == genome }      // channel: [ val(meta), [ bedgraph ] ]
    bigwig_endo      = UCSC_BEDGRAPHTOBIGWIG_ENDO.out.bigwig    // channel: [ val(meta), [ bigwig ] ]
    bigwig_exo       = UCSC_BEDGRAPHTOBIGWIG_EXO.out.bigwig    // channel: [ val(meta), [ bigwig ] ]
    bigwig_cmp_endo  = ch_bw_compare.filter { it -> it[0].genome == genome }      // channel: [ val(meta), [ bigwig ] ]
    bigwig_cmp_exo   = ch_bw_compare.filter { it -> it[0].genome == spikein_genome }      // channel: [ val(meta), [ bigwig ] ]
    bigwig_avg_endo  = ch_bw_avg.filter { it -> it[0].genome == genome }      // channel: [ val(meta), [ bigwig ] ]
    bigwig_avg_exo   = ch_bw_avg.filter { it -> it[0].genome == spikein_genome }      // channel: [ val(meta), [ bigwig ] ]
    bigwig_all_endo  = ch_bigwig_all_endo      // channel: [ val(meta), [ bigwig ] ]
    bigwig_all_exo   = ch_bw_all.filter { it -> it[0].genome == spikein_genome }      // channel: [ val(meta), [ bigwig ]
}
