process PARTITION_PLOT {
    tag "$meta.id"
    label 'process_medium_memory'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-genomicalignments_bioconductor-genomicfeatures_r-argparse_r-dplyr_pruned:aad6cf5716386302' :
        'community.wave.seqera.io/library/bioconductor-genomicalignments_bioconductor-genomicfeatures_r-argparse_r-dplyr_pruned:ff63fd989740e4c5' }"

    input:
    tuple val(meta), path(partition), path(strandedinput), path(scarminusinput), path(okazaki)
    tuple val(meta2), path(blacklist)
    tuple val(meta3), path(initiation_zones)

    output:
    path "*.scatter_plots.pdf",             optional:true, emit: scatter_pdf
    path "*.partition_plots_raw.pdf",                      emit: partition_raw_pdf
    path "*.partition_plots_smoothed.pdf",                 emit: partition_smoothed_pdf
    path "versions.yml",                                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args               = task.ext.args ?: ''
    def prefix             = task.ext.prefix ?: "${meta.id}"
    def scar_partition_arg = partition ? "--scar_partition_file $partition" : ''
    def scarminusinput_arg = scarminusinput ? "--scarminusinput_partition_file $scarminusinput" : ''
    def strandedinput_arg  = strandedinput ? "--strandedinput_partition_file $strandedinput" : ''
    def okazaki_arg        = okazaki ? "--okazaki_file $okazaki" : ''
    def iz_arg             = initiation_zones ? "--initiation_zones $initiation_zones" : ''
    def blacklist_arg      = blacklist ? "--blacklist $blacklist" : ''

    """
    partition_or_rfd_plot.R \\
        ${scar_partition_arg} \\
        ${scarminusinput_arg} \\
        ${strandedinput_arg} \\
        ${okazaki_arg} \\
        ${iz_arg} \\
        ${blacklist_arg} \\
        --prefix ${prefix} \\
        --outdir ./ \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch  ${prefix}.scatter_plots.pdf
    touch  ${prefix}.partition_plots_raw.pdf
    touch  ${prefix}.partition_plots_smoothed.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
