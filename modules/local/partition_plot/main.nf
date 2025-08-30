process PARTITION_PLOT {
    tag "$meta.id"
    label 'process_medium_memory'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/99/99b2a7149c1943265175ef013ca8c247c9847b14a0d2e802d0cbdda4a58458a5/data' :
        'community.wave.seqera.io/library/bioconductor-genomicalignments_bioconductor-genomicfeatures_r-argparse_r-ggpmisc_pruned:2c7c689513df97b4' }"

    input:
    tuple val(meta), path(partition), path(strandedinput), path(scarminusinput)
    tuple val(meta2), path(blacklist)
    tuple val(meta3), path(okseq_rfd_file)
    tuple val(meta4), path(initiation_zones)
    tuple val(meta5), path(chrom_sizes)

    output:
    path "*.partition_plots_raw.pdf",       emit: partition_raw_pdf
    path "*.partition_plots_raw.png",       emit: partition_raw_png
    path "*.partition_plots_smoothed.pdf",  emit: partition_smoothed_pdf
    path "*.partition_plots_smoothed.png",  emit: partition_smoothed_png
    path "*.scatter_plots.pdf",             emit: scatter_pdf, optional:true
    path "*.scatter_plots.png",             emit: scatter_png, optional:true
    path "versions.yml",                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args               = task.ext.args ?: ''
    def prefix             = task.ext.prefix ?: "${meta.id}"
    def scar_partition_arg = partition ? "--scar_partition_file $partition" : ''
    def scarminusinput_arg = scarminusinput ? "--scarminusinput_partition_file $scarminusinput" : ''
    def strandedinput_arg  = strandedinput ? "--strandedinput_partition_file $strandedinput" : ''
    def okazaki_arg        = okseq_rfd_file ? "--okazaki_file $okseq_rfd_file" : ''
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
        --chrom_sizes ${chrom_sizes} \\
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
    touch  ${prefix}.scatter_plots.png
    touch  ${prefix}.partition_plots_raw.pdf
    touch  ${prefix}.partition_plots_raw.png
    touch  ${prefix}.partition_plots_smoothed.pdf
    touch  ${prefix}.partition_plots_smoothed.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
