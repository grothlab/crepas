include { DEEPTOOLS_BAMCOVERAGE as DEEPTOOLS_BAMCOVERAGE_BINS }       from '../../../modules/nf-core/deeptools/bamcoverage/main'
include { BEDTOOLS_MAKEWINDOWS as BEDTOOLS_MAKEWINDOWS_ENDO }        from '../../../modules/nf-core/bedtools/makewindows/main'
include { BEDTOOLS_MAKEWINDOWS as BEDTOOLS_MAKEWINDOWS_EXO }        from '../../../modules/nf-core/bedtools/makewindows/main'
include { BEDTOOLS_MAP as BEDTOOLS_MAP_ENDO }                from '../../../modules/nf-core/bedtools/map/main'
include { BEDTOOLS_MAP as BEDTOOLS_MAP_EXO }                from '../../../modules/nf-core/bedtools/map/main'
include { BEDGRAPH_NORMALIZE }          from '../../../modules/local/bedgraph_normalize/main'
include { BEDGRAPH_SIGNAL_OVER_INPUT }  from '../../../modules/local/bedgraph_signal_over_input/main'
include { FILE_SORT as BEDGRAPH_SORT } from '../../../modules/local/file_sort/main'
include { UCSC_BEDGRAPHTOBIGWIG as UCSC_BEDGRAPHTOBIGWIG_ENDO }   from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'
include { UCSC_BEDGRAPHTOBIGWIG as UCSC_BEDGRAPHTOBIGWIG_EXO }   from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'
include { DEEPTOOLS_BAMCOVERAGE as DEEPTOOLS_BAMCOVERAGE_BINSIZE1 }     from '../../../modules/nf-core/deeptools/bamcoverage/main'


workflow BAM_NORMALIZE_BIGWIG_DEEPTOOLS {

    take:
    ch_bam_bai               // channel: [ val(meta), [ bam ], [ bai ] ]
    ch_chrom_sizes_endo          // channel: [ bed ]
    ch_chrom_sizes_exo          // channel: [ bed ]
    coverage_bin_size      // int: size of the coverage bin in bp
    genome                  // string: genome name
    spikein_genome          // string: spike-in genome name
    skip_srpm           // boolean: skip the SRPM normalization step
    skip_cisrpm          // boolean: skip the CISRPM normalization step
    skip_cisrpmsoi      // boolean: skip the CISRPM-SOI normalization step
    skip_plot_profile   // boolean: skip the plot profile step
    rpm_use_flT2_total  // string: comma-separated list of antibodies for which to use flT2_total_mapped_reads instead of flT3_total_mapped_reads for RPM normalization
    srpm_use_flT2_total // string: comma-separated list of antibodies for which to use flT2_total_mapped_reads instead of flT3_total_mapped_reads for SRPM normalization
    cisrpm_use_flT2_total // string: comma-separated list of antibodies for which to use flT2_total_mapped_reads instead of flT3_total_mapped_reads for CISRPM normalization

    main:

    ch_versions = Channel.empty()


    ch_bam_bai
        .map { meta, bam, bai ->
            def meta_clone = meta.clone()
            meta_clone.norm_factor_val = 1
            meta_clone.norm_factor_type = 'raw'
            [ meta_clone, bam, bai ]
        }
        .set { ch_bam_bai }

    //
    // MODULE: Calculate raw coverage per bin
    //
    DEEPTOOLS_BAMCOVERAGE_BINS (
        ch_bam_bai,
        [],
        []
    )
    ch_bdg_raw = DEEPTOOLS_BAMCOVERAGE_BINS.out.bedgraph
    ch_versions = ch_versions.mix(DEEPTOOLS_BAMCOVERAGE_BINS.out.versions.first())

    // bamCoverage merges contiguous bins with the same coverage and there is no option to disable this.
    // See https://github.com/deeptools/deepTools/issues/907#issuecomment-576729674
    // Therefore, to get bins of equal size across the whole genome,
    // we have to make windows and map the bedgraph to the windows:

    //
    // MODULE: Make windows of equal size across the whole genome
    //
    BEDTOOLS_MAKEWINDOWS_ENDO (
        ch_chrom_sizes_endo
    )
    ch_windows_endo = BEDTOOLS_MAKEWINDOWS_ENDO.out.bed
    ch_versions = ch_versions.mix(BEDTOOLS_MAKEWINDOWS_ENDO.out.versions.first())

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
    ch_versions = ch_versions.mix(BEDTOOLS_MAP_ENDO.out.versions.first())

    ch_windows_exo = Channel.empty()
    ch_windows_exo_bdg_raw = Channel.empty()
    if (spikein_genome) {
        BEDTOOLS_MAKEWINDOWS_EXO (
            ch_chrom_sizes_exo
        )
        ch_windows_exo = BEDTOOLS_MAKEWINDOWS_EXO.out.bed
        ch_versions = ch_versions.mix(BEDTOOLS_MAKEWINDOWS_EXO.out.versions.first())

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
        ch_versions = ch_versions.mix(BEDTOOLS_MAP_EXO.out.versions.first())

        // Merge the two channels
        ch_bdg_map = ch_bdg_map_endo.mix(ch_bdg_map_exo)
    }

    // TODO: print for debugging
    ch_bdg_map
        .map { meta, bdg ->
            "${meta}\t${bdg}"
        }
        .collectFile( name: 'ch_bdg_map.txt', newLine: true, sort: false, storeDir: "${params.outdir}" )

    // Copy and modify channel meta to add RPM normalization factors
    ch_bdg_rpm = Channel.empty()
    ch_bdg_map
        .map { meta, bdg ->
            def meta_clone = meta.clone()
            // if meta.antibody is in the list of antibodies to use flT2_total_mapped_reads or flT1_total_mapped_reads for RPM normalization, otherwise use flT3_total_mapped_reads
            if (rpm_use_flT2_total && meta.antibody in rpm_use_flT2_total.split(',').collect { it.trim() }) {
                if (meta_clone.flT2_total_mapped_reads) {
                    meta_clone.norm_factor_val = 1e6 / meta_clone.flT2_total_mapped_reads
                    meta_clone.norm_factor_val_used = 'flT2_total_mapped_reads'
                } else {
                    // Samples without spike-in wouldn't have flT2_total_mapped_reads, so we use flT1_total_mapped_reads instead
                    meta_clone.norm_factor_val = 1e6 / meta_clone.flT1_total_mapped_reads
                    meta_clone.norm_factor_val_used = 'flT1_total_mapped_reads'
                }
            } else {
                meta_clone.norm_factor_val = 1e6 / meta_clone.flT3_total_mapped_reads
                meta_clone.norm_factor_val_used = 'flT3_total_mapped_reads'
            }
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
        .collectFile( name: 'ch_bdg_rpm.txt', newLine: true, sort: false, storeDir: "${params.outdir}" )

    // Copy and modify channel meta to add SRPM normalization factors
    ch_bdg_srpm = Channel.empty()
    if (!skip_srpm) {
        ch_bdg_map
            .map { meta, bdg ->
                    [ meta.id, meta, bdg ]
            }
            .branch { id, meta, bdg ->
                endo: meta.genome == genome
                exo: meta.genome == spikein_genome
            }
            .set { ch_bdg_genome }

        ch_bdg_genome.endo
            .combine(ch_bdg_genome.exo, by: 0)
            .map { id, endo_meta, endo_bdg, exo_meta, exo_bdg ->
                def meta_clone = endo_meta.clone()
                if (srpm_use_flT2_total && meta_clone.antibody in srpm_use_flT2_total.split(',').collect { it.trim() }) {
                    meta_clone.norm_factor_val = 1e6 / exo_meta.flT2_total_mapped_reads
                    meta_clone.norm_factor_val_used = 'flT2_total_mapped_reads'
                } else {
                    meta_clone.norm_factor_val = 1e6 / exo_meta.flT3_total_mapped_reads
                    meta_clone.norm_factor_val_used = 'flT3_total_mapped_reads'
                }
                meta_clone.norm_factor_type = 'srpm'
                [ meta_clone, endo_bdg ]
            }
            .set { ch_bdg_srpm }
    }

    // TODO: print for debugging
    ch_bdg_srpm
        .map {
            meta, bdg ->
                "${meta}\t${bdg}"
        }
        .collectFile( name: 'ch_bdg_srpm.txt', newLine: true, sort: false, storeDir: "${params.outdir}" )

    
    // Copy and modify channel meta to add CISRPM normalization factors
    ch_bdg_genome_type = Channel.empty()
    ch_bdg_genome_ip = Channel.empty()
    ch_bdg_genome_control = Channel.empty()
    ch_bdg_ip_cisrpm = Channel.empty()
    ch_bdg_control_cisrpm = Channel.empty()
    ch_bdg_cisrpm = Channel.empty()
    if (!skip_cisrpm) {

        // Split BAMs by genome (endo and exo) and by type (ip and control)
        ch_bdg_map
            .map { meta, bdg ->
                [ meta.id, meta, bdg ]
            }
            .branch { id, meta, bdg ->
                endo_ip: meta.genome == genome && !meta.is_control
                endo_control: meta.genome == genome && meta.is_control
                exo_ip: meta.genome == spikein_genome && !meta.is_control
                exo_control: meta.genome == spikein_genome && meta.is_control
            }
            .set { ch_bdg_genome_type }

        // Combine the endo and exo BAMs (ChIPs)
        ch_bdg_genome_type.endo_ip
            .combine(ch_bdg_genome_type.exo_ip, by: 0)
            .map { ip_id, endo_ip_meta, endo_ip_bdg, exo_ip_meta, exo_ip_bdg ->
                    [ endo_ip_meta.control, endo_ip_meta, endo_ip_bdg, exo_ip_meta, exo_ip_bdg ]
            }
            .set { ch_bdg_genome_ip }

        // TODO: print for debugging
        ch_bdg_genome_ip
            .map {
                control_id, endo_ip_meta, endo_ip_bdg, exo_ip_meta, exo_ip_bdg ->
                    "${control_id}\t${endo_ip_meta}\t${endo_ip_bdg}\t${exo_ip_meta}\t${exo_ip_bdg}"
            }
            .collectFile( name: 'ch_bdg_genome_ip.txt', newLine: true, sort: false, storeDir: "${params.outdir}" )

        // Combine the endo and exo BAMs (inputs)
        ch_bdg_genome_type.endo_control
            .combine(ch_bdg_genome_type.exo_control, by: 0)
            .map { control_id, endo_control_meta, endo_control_bdg, exo_control_meta, exo_control_bdg ->
                [ endo_control_meta.id, endo_control_meta, endo_control_bdg, exo_control_meta, exo_control_bdg ]
            }
            .set { ch_bdg_genome_control }
        
        // TODO: print for debugging
        ch_bdg_genome_control
            .map {
                control_id, endo_control_meta, endo_control_bdg, exo_control_meta, exo_control_bdg ->
                    "${control_id}\t${endo_control_meta}\t${endo_control_bdg}\t${exo_control_meta}\t${exo_control_bdg}"
            }
            .collectFile( name: 'ch_bdg_genome_control.txt', newLine: true, sort: false, storeDir: "${params.outdir}" )

        // Combine the combined ChIPs with the combined inputs
        ch_bdg_genome_ip
            .combine(ch_bdg_genome_control, by: 0)
            .map { id, endo_ip_meta, endo_ip_bdg, exo_ip_meta, exo_ip_bdg, endo_control_meta, endo_control_bdg, exo_control_meta, exo_control_bdg ->
                    def meta_clone = endo_ip_meta.clone()
                    if (cisrpm_use_flT2_total && meta_clone.antibody in cisrpm_use_flT2_total.split(',').collect { it.trim() }) {
                        meta_clone.norm_factor_val = (1e6 / exo_ip_meta.flT2_total_mapped_reads) * (exo_control_meta.flT2_total_mapped_reads / endo_control_meta.flT2_total_mapped_reads)
                        meta_clone.norm_factor_val_used = 'flT2_total_mapped_reads'
                    } else {
                        meta_clone.norm_factor_val = (1e6 / exo_ip_meta.flT3_total_mapped_reads) * (exo_control_meta.flT3_total_mapped_reads / endo_control_meta.flT3_total_mapped_reads)
                        meta_clone.norm_factor_val_used = 'flT3_total_mapped_reads'
                    }
                    meta_clone.norm_factor_type = 'cisrpm'
                    [ meta_clone, endo_ip_bdg ]
            }
            .set { ch_bdg_ip_cisrpm }

        // TODO: print for debugging
        ch_bdg_ip_cisrpm
            .map {
                meta, bdg ->
                    "${meta}\t${bdg}"
            }
            .collectFile( name: 'ch_bdg_ip_cisrpm.txt', newLine: true, sort: false, storeDir: "${params.outdir}" )
        
        // Now do the missing CISRPM for the endogenous inputs
        // In this case CISRPM is the same as RPM, but we cannot
        // just copy the RPM from before, since cisrpm_use_flT2_total
        // can be different than rpm_use_flT2_total
        ch_bdg_map
            .filter { meta, bdg ->
                meta.genome == genome && meta.is_control
            }
            .map { meta, bdg ->
                def meta_clone = meta.clone()
                if (cisrpm_use_flT2_total && meta.antibody in cisrpm_use_flT2_total.split(',').collect { it.trim() }) {
                    meta_clone.norm_factor_val = 1e6 / meta_clone.flT2_total_mapped_reads
                    meta_clone.norm_factor_val_used = 'flT2_total_mapped_reads'
                } else {
                    meta_clone.norm_factor_val = 1e6 / meta_clone.flT3_total_mapped_reads
                    meta_clone.norm_factor_val_used = 'flT3_total_mapped_reads'
                }
                meta_clone.norm_factor_type = 'cisrpm'
                [ meta_clone, bdg ]
            }
            .set { ch_bdg_control_cisrpm }


        // TODO: print for debugging
        ch_bdg_control_cisrpm
            .map {
                meta, bdg ->
                    "${meta}\t${bdg}"
            }
            .collectFile( name: 'ch_bdg_control_cisrpm.txt', newLine: true, sort: false, storeDir: "${params.outdir}" )

        ch_bdg_cisrpm = ch_bdg_ip_cisrpm.mix(ch_bdg_control_cisrpm)
    }

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
        .collectFile( name: 'ch_bdg_norm.txt', newLine: true, sort: false, storeDir: "${params.outdir}" )

    //
    // MODULE: Normalize the bedgraphs
    //
    BEDGRAPH_NORMALIZE (
        ch_bdg_norm
    )
    ch_versions = ch_versions.mix(BEDGRAPH_NORMALIZE.out.versions.first())

    ch_bdg_map_norm = ch_bdg_map.mix(BEDGRAPH_NORMALIZE.out.bedgraph)

    // Create channel: [ val(meta), [ ip_bdg ], [ control_bdg ] ]
    ch_bdg_ip_control_cisrpm = Channel.empty()
    if (!skip_cisrpmsoi) {
        ch_bdg_map_norm
            .filter { meta, bdg ->
                meta.norm_factor_type == 'cisrpm'
            }
            .branch { meta, bdg ->
                ips_with_control: meta.control
                    return [ meta.control, meta, bdg ]
                    // Cannot calculate CISRPM-SOI for ChIPs without inputs
                    // ips_without_control: !meta.control && !meta.is_control
                    //     return [ meta, bdg ]
                controls: !meta.control && meta.is_control
                    return [ meta.id, bdg ]
            }
            .set { ch_bdg_ip_control_cisrpm }

        ch_bdg_ip_control_cisrpm
            .ips_with_control
            .combine(ch_bdg_ip_control_cisrpm.controls, by: 0)
            .map { control_id, ip_meta, ip_bdg, control_bdg ->
                def meta_clone = ip_meta.clone()
                    meta_clone.signal_over_input = true
                    [ meta_clone, ip_bdg, control_bdg ]
            }
            .set { ch_bdg_ip_control_cisrpm }

        //
        // MODULE: Calculate CIRSPM signal over input (ChIP over input)
        //
        BEDGRAPH_SIGNAL_OVER_INPUT (
            ch_bdg_ip_control_cisrpm
        )
        ch_bdg_all = BEDGRAPH_SIGNAL_OVER_INPUT.out.bedgraph.mix(ch_bdg_map_norm)
        ch_versions = ch_versions.mix(BEDGRAPH_SIGNAL_OVER_INPUT.out.versions.first())

    }

    //
    // MODULE: Sort the bedgraph so that it works with ucsc_bedgraphtobigwig
    //
    BEDGRAPH_SORT (
        ch_bdg_all,
        'bedgraph'
    )
    ch_bdg_all = BEDGRAPH_SORT.out.sorted
    ch_versions = ch_versions.mix(BEDGRAPH_SORT.out.versions.first())

  
    // Remove empty bedgraphs (mostly from SOI) to avoid bdg to bw conversion errors
    ch_bdg_all
        .filter {
            meta, bdg ->
                bdg.size() > 0
        }
        .set { ch_bdg_all }

    // TODO: print for debugging
    ch_bdg_all
        .map {
            meta, bdg ->
                "${meta}\t${bdg}"
        }
        .collectFile( name: 'ch_bdg_all.txt', newLine: true, sort: false, storeDir: "${params.outdir}" )


    //
    // MODULE: Convert bedgraph to bigwig
    //
    UCSC_BEDGRAPHTOBIGWIG_ENDO (
        ch_bdg_all.filter { it -> it[0].genome == genome },
        ch_chrom_sizes_endo.map { it[1] }
    )
    ch_versions = ch_versions.mix(UCSC_BEDGRAPHTOBIGWIG_ENDO.out.versions.first())
   
   ch_bw_exo = Channel.empty()
   if (spikein_genome) {
        UCSC_BEDGRAPHTOBIGWIG_EXO (
            ch_bdg_all.filter { it -> it[0].genome == spikein_genome },
            ch_chrom_sizes_exo.map { it[1] }
        )
        ch_bw_exo = UCSC_BEDGRAPHTOBIGWIG_EXO.out.bigwig
        ch_versions = ch_versions.mix(UCSC_BEDGRAPHTOBIGWIG_EXO.out.versions.first())
    }


    // if coverage_bin_size is not 1, then we need to generate bw with that binsize for computeMatrix
    ch_binsize1 = Channel.empty()
    if (coverage_bin_size != 1 && !skip_plot_profile) {

        ch_bam_bai
            // we only want these bw for the endo genome
            .filter { meta, bam, bai ->
                meta.genome == genome
            }
            .map { meta, bam, bai ->
                def meta_clone = meta.clone()
                // if meta.antibody is in the list of antibodies to use flT2_total_mapped_reads or flT1_total_mapped_reads for RPM normalization, otherwise use flT3_total_mapped_reads
                if (rpm_use_flT2_total && meta.antibody in rpm_use_flT2_total.split(',').collect { it.trim() }) {
                    if (meta_clone.flT2_total_mapped_reads) {
                        meta_clone.norm_factor_val = 1e6 / meta_clone.flT2_total_mapped_reads
                        meta_clone.norm_factor_val_used = 'flT2_total_mapped_reads'
                    } else {
                        // Samples without spike-in wouldn't have flT2_total_mapped_reads, so we use flT1_total_mapped_reads instead
                        meta_clone.norm_factor_val = 1e6 / meta_clone.flT1_total_mapped_reads
                        meta_clone.norm_factor_val_used = 'flT1_total_mapped_reads'
                    }
                } else {
                    meta_clone.norm_factor_val = 1e6 / meta_clone.flT3_total_mapped_reads
                    meta_clone.norm_factor_val_used = 'flT3_total_mapped_reads'
                }
                meta_clone.norm_factor_type = 'rpm'
                    [ meta_clone, bam, bai ]
            }
            .set { ch_binsize1 }
        
    
        DEEPTOOLS_BAMCOVERAGE_BINSIZE1 (
            ch_binsize1,
            [],
            []
        )
        ch_binsize1 = DEEPTOOLS_BAMCOVERAGE_BINSIZE1.out.bigwig
        ch_versions = ch_versions.mix(DEEPTOOLS_BAMCOVERAGE_BINSIZE1.out.versions.first())

    }

    emit:
    bigwig_endo_rpm  = UCSC_BEDGRAPHTOBIGWIG_ENDO.out.bigwig.filter { it -> it[0].norm_factor_type == 'rpm' } // channel: [ val(meta), [ bigwig ] ]
    bigwig_endo      = UCSC_BEDGRAPHTOBIGWIG_ENDO.out.bigwig        // channel: [ val(meta), [ bigwig ] ]
    bigwig_exo       = ch_bw_exo                                    // channel: [ val(meta), [ bigwig ] ]
    bigwig_binsize1  = ch_binsize1                                  // channel: [ val(meta), [ bigwig ] ]

    versions      = ch_versions                                     // channel: [ versions.yml ]
}

